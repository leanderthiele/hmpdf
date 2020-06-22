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

int
null_power(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

    d->pwr->inited_power = 0;
    d->pwr->ssq = NULL;
    d->pwr->created_corr = 0;
    d->pwr->corr_interp = NULL;
    d->pwr->corr_accel = NULL;

    CHECKERR
    return hmpdf_status;
}//}}}

int
reset_power(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

    if (d->pwr->ssq != NULL)
    {
        for (int M_index=0; M_index<d->n->NM; M_index++)
        {
            if (d->pwr->ssq[M_index] != NULL)
            {
                free(d->pwr->ssq[M_index]);
            }
        }
        free(d->pwr->ssq);
    }
    if (d->pwr->corr_interp != NULL) { gsl_spline_free(d->pwr->corr_interp); }
    if (d->pwr->corr_accel != NULL)
    {
        for (int ii=0; ii<d->pwr->Ncorr_accel; ii++)
        {
            gsl_interp_accel_free(d->pwr->corr_accel[ii]);
        }
        free(d->pwr->corr_accel);
    }

    CHECKERR
    return hmpdf_status;
}//}}}

int
Pk_linear(hmpdf_obj *d, double k, double *out)
// k is logk if LOGK is defined
{//{{{
    int hmpdf_status = 0;

    #ifdef LOGK
    k = exp(k);
    #endif

    struct background *ba = (struct background *)d->cls->ba;
    struct primordial *pm = (struct primordial *)d->cls->pm;
    struct nonlinear *nl = (struct nonlinear *)d->cls->nl;

    if (nonlinear_pk_at_k_and_z(ba, pm, nl,
                                pk_linear, k, 0.0,
                                nl->index_pk_total,
                                out, NULL) == _FAILURE_)
    {
        double lnk = log(k);
        if (lnk > nl->ln_k[nl->k_size-1])
        // check if k falls out of bounds,
        //    maximum k needs to be set large enough that linear power is
        //    essentially zero for larger wavenumbers
        {
            *out = 0.0;
        }
        else
        {
            fprintf(stderr, "CLASS error %s\n", nl->error_message);
            fflush(stderr);
            ERRLOC;
            hmpdf_status |= 1;
        }
    }

    CHECKERR
    return hmpdf_status;
}//}}}

static int
power_kernel(hmpdf_obj *d, double k, double *out)
// kernel = k^2 P(k) / 2\pi^2
{//{{{
    int hmpdf_status = 0;

    SAFEHMPDF(Pk_linear(d, k, out));
    #ifdef LOGK
    *out *= exp(3.0*k)/2.0/M_PI/M_PI;
    #else
    *out *= k*k/2.0/M_PI/M_PI;
    #endif

    CHECKERR
    return hmpdf_status;
}//}}}

typedef struct
{//{{{
    int hmpdf_status;
    hmpdf_obj *d;
    gsl_function F;
}//}}}
power_integrand_params;

static double
power_integrand(double k, void *params)
// this is the gsl_function
{//{{{
    power_integrand_params *p = (power_integrand_params *)params;
    double temp;
    p->hmpdf_status = power_kernel(p->d, k, &temp);
    double out = GSL_FN_EVAL(&(p->F), k);
    return temp * out;
}//}}}

static int
power_integral(hmpdf_obj *d, power_integrand_params *p, double *out)
{//{{{
    int hmpdf_status = 0;
    p->hmpdf_status = 0;

    struct nonlinear *nl = (struct nonlinear *)d->cls->nl;

    gsl_function integrand;
    integrand.function = &power_integrand;
    integrand.params = p;
    double err;
    #ifdef LOGK
    double kmin = log(PKINTEGR_KMIN);
    double kmax = nl->ln_k[nl->k_size-1];
    #else
    double kmin = PKINTEGR_KMIN;
    double kmax = exp(nl->ln_k[nl->k_size-1]);
    #endif

    SAFEALLOC(gsl_integration_workspace *, ws,
              gsl_integration_workspace_alloc(PKINTEGR_LIMIT))
    SAFEGSL(gsl_integration_qag(&integrand, kmin, kmax,
                                PKINTEGR_EPSABS, PKINTEGR_EPSREL,
                                PKINTEGR_LIMIT,  PKINTEGR_KEY,
                                ws, out, &err))
    gsl_integration_workspace_free(ws);

    hmpdf_status |= p->hmpdf_status;

    CHECKERR
    return hmpdf_status;
}//}}}

static double
tophat_Wsq(double k, void *params)
{//{{{
    #ifdef LOGK
    k = exp(k);
    #endif
    double R = *(double *)params;
    double x = k * R;
    return gsl_pow_2(3.0*(sin(x) - x*cos(x))/gsl_pow_3(x));
}//}}}

static double
tophat_Wsqprime(double k, void *params)
{//{{{
    #ifdef LOGK
    k = exp(k);
    #endif
    double R = *(double *)params;
    double x = k * R;
    return k * 6.0 * (sin(x) - x*cos(x))/gsl_pow_3(x)
                   * (3.0*sin(x)/gsl_pow_2(x) - 9.0*(sin(x) - x*cos(x))/gsl_pow_4(x));
}//}}}

static int 
_ssq(hmpdf_obj *d, double M, double *ssq, double *dssq)
// return sigma^2(M), write d sigma^2 / dlogM into return value
{//{{{
    int hmpdf_status = 0;

    double R = cbrt(3.0*M/4.0/M_PI/d->c->rho_m_0);

    power_integrand_params p;
    p.d = d;
    p.F.params = &R;
    // cmpute derivative d sigma^2 / dlogM
    p.F.function = tophat_Wsqprime;
    SAFEHMPDF(power_integral(d, &p, dssq))
    *dssq *= (R/3.0);
    // compute sigma^2(M)
    p.F.function = tophat_Wsq;
    SAFEHMPDF(power_integral(d, &p, ssq))

    CHECKERR
    return hmpdf_status;
}//}}}

static int
create_ssq(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

    fprintf(stdout, "\tcreate_ssq\n");
    fflush(stdout);

    SAFEALLOC(, d->pwr->ssq, malloc(d->n->NM * sizeof(double*)))
    for (int M_index=0; M_index<d->n->NM; M_index++)
    {
        SAFEALLOC(, d->pwr->ssq[M_index], malloc(2 * sizeof(double)))
        SAFEHMPDF(_ssq(d, d->n->Mgrid[M_index],
                       d->pwr->ssq[M_index]+0,
                       d->pwr->ssq[M_index]+1))
    }

    CHECKERR
    return hmpdf_status;
}//}}}

static double
autocorr_kernel(double k, void *params)
{//{{{
    #ifdef LOGK
    return M_PI * exp(-k);
    #else
    return M_PI * 1.0/k;
    #endif
}//}}}

static int
_autocorr(hmpdf_obj *d, double *out)
{//{{{
    int hmpdf_status = 0;

    power_integrand_params p;
    p.d = d;
    p.F.params = NULL;
    p.F.function = autocorr_kernel;
    SAFEHMPDF(power_integral(d, &p, out))

    CHECKERR
    return hmpdf_status;
}//}}}

static int
create_autocorr(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

    fprintf(stdout, "\tcreate_autocorr\n");
    fflush(stdout);
    SAFEHMPDF(_autocorr(d, &(d->pwr->autocorr)))

    CHECKERR
    return hmpdf_status;
}//}}}

int
create_corr(hmpdf_obj *d)
// computes the z=0 matter correlation function
{//{{{
    int hmpdf_status = 0;

    if (d->pwr->created_corr) { return hmpdf_status; }

    fprintf(stdout, "\tcreate_corr_interp\n");
    fflush(stdout);
    
    double rmax = 1.1 * d->n->phimax * d->c->comoving[d->n->Nz-1];
    gsl_dht *t = gsl_dht_new(CORRINTERP_N, 0, rmax);
    SAFEALLOC(double *, Pk,   malloc(CORRINTERP_N     * sizeof(double)))
    SAFEALLOC(double *, r,    malloc((CORRINTERP_N+1) * sizeof(double)))
    SAFEALLOC(double *, zeta, malloc((CORRINTERP_N+1) * sizeof(double)))
    r[0] = 0.0;
    zeta[0] = d->pwr->autocorr;

    double hankel_norm = gsl_pow_2(gsl_dht_k_sample(t, 0)
                                   / gsl_dht_x_sample(t, 0));
    for (int ii=0; ii<CORRINTERP_N; ii++)
    {
        #ifdef LOGK
        double k = log(gsl_dht_k_sample(t, ii));
        #else
        double k = gsl_dht_k_sample(t, ii);
        #endif
        SAFEHMPDF(Pk_linear(d, k, Pk+ii))
        r[ii+1] = gsl_dht_x_sample(t, ii);
    }
    SAFEGSL(gsl_dht_apply(t, Pk, zeta+1))
    gsl_dht_free(t);
    free(Pk);

    // divide by 2\pi and fix Hankel normalization
    for (int ii=1; ii<=CORRINTERP_N; ii++)
    {
        zeta[ii] *= 0.5 * M_1_PI * hankel_norm;
    }

    SAFEALLOC(, d->pwr->corr_interp,
              gsl_spline_alloc(gsl_interp_cspline, CORRINTERP_N+1))
    d->pwr->Ncorr_accel = d->Ncores;
    SAFEALLOC(, d->pwr->corr_accel,
              malloc(d->pwr->Ncorr_accel * sizeof(gsl_interp_accel *)))
    for (int ii=0; ii<d->pwr->Ncorr_accel; ii++)
    {
        SAFEALLOC(, d->pwr->corr_accel[ii], gsl_interp_accel_alloc())
    }
    SAFEGSL(gsl_spline_init(d->pwr->corr_interp, r,
                            zeta, CORRINTERP_N+1))
    free(r);
    free(zeta);

    d->pwr->created_corr = 1;

    CHECKERR
    return hmpdf_status;
}//}}}

int
corr(hmpdf_obj *d, int z_index, double phi, double *out)
{//{{{
    int hmpdf_status = 0;

    double r = d->c->comoving[z_index] * phi;
    SAFEGSL(gsl_spline_eval_e(d->pwr->corr_interp, r,
                              d->pwr->corr_accel[this_core()],
                              out))
    *out *= d->c->Dsq[z_index];
    
    CHECKERR
    return hmpdf_status;
}//}}}

int
init_power(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

    fprintf(stdout, "In power.h -> init_power :\n");
    fflush(stdout);
    SAFEHMPDF(create_ssq(d))
    SAFEHMPDF(create_autocorr(d))

    CHECKERR
    return hmpdf_status;
}//}}}

