#include <stdio.h>
#include <stdlib.h>

#include <class.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#include "configs.h"
#include "utils.h"
#include "data.h"
#include "cosmology.h"
#include "numerics.h"
#include "power.h"

struct class_interface_s
{//{{{
    char *class_ini;
    char *class_pre;

    struct precision *pr;
    struct background *ba;
    struct thermo *th;
    struct primordial *pm;
    struct perturbs *pt;
    struct nonlinear *nl;
    struct transfers *tr;
    struct spectra *sp;
    struct lensing *le;
    struct output *op;

    ErrorMsg errmsg;
};//}}}

typedef struct class_interface_s cls;

static
double Pk_linear(all_data *d, double k)
// k is logk if LOGK is defined
{//{{{
    #ifdef LOGK
    k = exp(k);
    #endif
    double out;
    cls *_c = (cls *)d->cls;
    nonlinear_pk_at_k_and_z(_c->ba, _c->pm, _c->nl, pk_linear,
                            k, 0.0, _c->nl->index_pk_total, &out, NULL);
    return out;
}//}}}

static
double power_kernel(all_data *d, double k)
// kernel = k^2 P(k) / 2\pi^2
{//{{{
    #ifdef LOGK
    return exp(3.0*k)*Pk_linear(d, k)/2.0/M_PI/M_PI;
    #else
    return k*k*Pk_linear(d, k)/2.0/M_PI/M_PI;
    #endif
}//}}}

typedef struct
{//{{{
    all_data *d;
    gsl_function F;
}//}}}
power_integrand_params;

static
double power_integrand(double k, void *params)
// this is the gsl_function
{//{{{
    power_integrand_params *p = (power_integrand_params *)params;
    return power_kernel(p->d, k) * GSL_FN_EVAL(&(p->F), k);
}//}}}

static
double power_integral(power_integrand_params *p)
{//{{{
    gsl_function integrand;
    integrand.function = &power_integrand;
    integrand.params = p;
    double res, err;
    #ifdef LOGK
    double kmin = log(PKINTEGR_KMIN);
    double kmax = log(PKINTEGR_KMAX);
    #else
    double kmin = PKINTEGR_KMIN;
    double kmax = PKINTEGR_KMAX;
    #endif

    gsl_integration_workspace *ws = gsl_integration_workspace_alloc(PKINTEGR_LIMIT);
    gsl_integration_qag(&integrand, kmin, kmax,
                        PKINTEGR_EPSABS, PKINTEGR_EPSREL,
                        PKINTEGR_LIMIT,  PKINTEGR_KEY,
                        ws, &res, &err);
    gsl_integration_workspace_free(ws);

    return res;
}//}}}

static
double tophat_Wsq(double k, void *params)
{//{{{
    #ifdef LOGK
    k = exp(k);
    #endif
    double R = *(double *)params;
    double x = k * R;
    return gsl_pow_2(3.0*(sin(x) - x*cos(x))/gsl_pow_3(x));
}//}}}

static
double tophat_Wsqprime(double k, void *params)
{//{{{
    #ifdef LOGK
    k = exp(k);
    #endif
    double R = *(double *)params;
    double x = k * R;
    return k * 6.0 * (sin(x) - x*cos(x))/gsl_pow_3(x)
                   * (3.0*sin(x)/gsl_pow_2(x) - 9.0*(sin(x) - x*cos(x))/gsl_pow_4(x));
}//}}}

static
double _ssq(all_data *d, double M, double *dssq)
// return sigma^2(M), write d sigma^2 / dlogM into return value
{//{{{
    double R = cbrt(3.0*M/4.0/M_PI/d->c->rho_m_0);

    power_integrand_params p;
    p.d = d;
    p.F.params = &R;
    // cmpute derivative d sigma^2 / dlogM
    p.F.function = tophat_Wsqprime;
    *dssq = (R/3.0) * power_integral(&p);
    // compute sigma^2(M)
    p.F.function = tophat_Wsq;
    return power_integral(&p);
}//}}}

static
void create_ssq(all_data *d)
{//{{{
    printf("\tcreate_ssq\n");
    d->pwr->ssq = (double **)malloc(d->n->gr->NM * sizeof(double*));
    for (int M_index=0; M_index<d->n->gr->NM; M_index++)
    {
        d->pwr->ssq[M_index] = (double *)malloc(2 * sizeof(double));
        d->pwr->ssq[M_index][0] = _ssq(d, d->n->gr->Mgrid[M_index],
                                       d->pwr->ssq[M_index]+1);
    }
}//}}}

static
double autocorr_kernel(double k, void *params)
{//{{{
    #ifdef LOGK
    return M_PI * exp(-k);
    #else
    return M_PI * 1.0/k;
    #endif
}//}}}

static
double _autocorr(all_data *d)
{//{{{
    power_integrand_params p;
    p.d = d;
    p.F.params = NULL;
    p.F.function = autocorr_kernel;
    return power_integral(&p);
}//}}}

static
void create_autocorr(all_data *d)
{//{{{
    printf("\tcreate_autocorr\n");
    d->pwr->autocorr = _autocorr(d);
}//}}}

void create_corr(all_data *d)
// computes the z=0 matter correlation function
{//{{{
    if (d->pwr->created_corr) { return; }
    printf("\tcreate_corr_interp\n");
    double rmax = 1.1 * d->n->gr->phimax * d->c->comoving[d->n->gr->Nz-1];
    gsl_dht *t = gsl_dht_new(CORRINTERP_N, 0, rmax);
    double *Pk = (double *)malloc(CORRINTERP_N * sizeof(double));
    double *r = (double *)malloc((CORRINTERP_N+1) * sizeof(double));
    double *zeta = (double *)malloc((CORRINTERP_N+1) * sizeof(double));
    r[0] = 0.0;
    zeta[0] = d->pwr->autocorr;

    double hankel_norm = gsl_pow_2(gsl_dht_k_sample(t,0)/gsl_dht_x_sample(t,0));
    for (int ii=0; ii<CORRINTERP_N; ii++)
    {
        #ifdef LOGK
        double k = log(gsl_dht_k_sample(t, ii));
        #else
        double k = gsl_dht_k_sample(t, ii);
        #endif
        Pk[ii] = Pk_linear(d, k);
        r[ii+1] = gsl_dht_x_sample(t, ii);
    }
    gsl_dht_apply(t, Pk, zeta+1);
    gsl_dht_free(t);
    free(Pk);

    // divide by 2\pi and fix Hankel normalization
    for (int ii=1; ii<=CORRINTERP_N; ii++)
    {
        zeta[ii] *= 0.5 * M_1_PI * hankel_norm;
    }

    d->pwr->corr_interp = gsl_spline_alloc(gsl_interp_cspline, CORRINTERP_N+1);
    d->pwr->Ncorr_accel = d->Ncores;
    d->pwr->corr_accel = (gsl_interp_accel **)malloc(d->pwr->Ncorr_accel
                                                     * sizeof(gsl_interp_accel *));
    for (int ii=0; ii<d->pwr->Ncorr_accel; ii++)
    {
        d->pwr->corr_accel[ii] = gsl_interp_accel_alloc();
    }
    gsl_spline_init(d->pwr->corr_interp, r, zeta, CORRINTERP_N+1);
    free(r);
    free(zeta);

    d->pwr->created_corr = 1;
}//}}}

double corr(all_data *d, int z_index, double phi)
// FIXME this is called in multithreaded mode, need one accelerator for each core
{//{{{
    double r = d->c->comoving[z_index] * phi;
    return d->c->Dsq[z_index]
           * gsl_spline_eval(d->pwr->corr_interp, r,
                             d->pwr->corr_accel[this_core()]);
}//}}}

void init_power(all_data *d)
{//{{{
    printf("In power.h -> init_power :\n");
    create_ssq(d);
    create_autocorr(d);
}//}}}

