#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_dht.h>

#include "configs.h"
#include "utils.h"
#include "object.h"
#include "numerics.h"
#include "filter.h"
#include "noise.h"

#include "hmpdf.h"

int
null_noise(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->ns->toepl = NULL;

    d->ns->created_zeta_interp = 0;
    d->ns->zeta_interp = NULL;
    d->ns->zeta_accel = NULL;

    ENDFCT
}//}}}

int
reset_noise(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\treset_noise\n");

    if (d->ns->toepl != NULL) { free(d->ns->toepl); }
    if (d->ns->zeta_interp != NULL) { gsl_spline_free(d->ns->zeta_interp); }
    if (d->ns->zeta_accel != NULL)
    {
        for (int ii=0; ii<d->ns->Nzeta_accel; ii++)
        {
            gsl_interp_accel_free(d->ns->zeta_accel[ii]);
        }
        free(d->ns->zeta_accel);
    }

    ENDFCT
}//}}}

static int
create_noisy_grids(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\tcreate_noisy_grids\n");

    // can in principle make this a user setting
    d->ns->len_kernel = d->n->Nsignal;
    // construct the new signal grid
    d->n->Nsignal_noisy = d->n->Nsignal+2*d->ns->len_kernel;
    SAFEALLOC(d->n->signalgrid_noisy, malloc(d->n->Nsignal_noisy
                                             * sizeof(double)));
    double extra_signal = (double)(d->ns->len_kernel)/(double)(d->n->Nsignal-1)
                          *(d->n->signalmax - d->n->signalmin);
    double smin = d->n->signalmin - extra_signal;
    double smax = d->n->signalmax + extra_signal;
    SAFEHMPDF(linspace(d->n->Nsignal_noisy, smin, smax, d->n->signalgrid_noisy));

    ENDFCT
}//}}}

int
create_noise_zeta_interp(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->ns->created_zeta_interp) { return 0; }

    gsl_dht *t;
    SAFEALLOC(t, gsl_dht_new(NOISE_ZETAINTERP_N, 0, d->n->phimax));
    double *ell;
    double *Nell;
    double *phi;
    double *zeta;
    SAFEALLOC(ell,  malloc(NOISE_ZETAINTERP_N * sizeof(double)));
    SAFEALLOC(Nell, malloc(NOISE_ZETAINTERP_N * sizeof(double)));
    SAFEALLOC(phi,  malloc((NOISE_ZETAINTERP_N+1) * sizeof(double)));
    SAFEALLOC(zeta, malloc((NOISE_ZETAINTERP_N+1) * sizeof(double)));
    phi[0] = 0.0;
    zeta[0] = d->ns->sigmasq;

    double hankel_norm = gsl_pow_2(gsl_dht_k_sample(t, 0)
                                   / gsl_dht_x_sample(t, 0));

    // fill the integrand
    for (int ii=0; ii<NOISE_ZETAINTERP_N; ii++)
    {
        ell[ii] = gsl_dht_k_sample(t, ii);
        Nell[ii] = d->ns->noise_pwr(ell[ii], d->ns->noise_pwr_params);
        phi[ii+1] = gsl_dht_x_sample(t, ii);
    }

    // add the filter functions
    SAFEHMPDF(apply_filters(d, NOISE_ZETAINTERP_N, ell, Nell, Nell,
                            1, filter_ps, NULL));

    // perform the hankel transform
    SAFEGSL(gsl_dht_apply(t, Nell, zeta+1));
    free(ell);
    free(Nell);

    // fix the normalization
    for (int ii=1; ii<=NOISE_ZETAINTERP_N; ii++)
    {
        zeta[ii] *= 0.5 * M_1_PI * hankel_norm;
    }

    // interpolate
    SAFEALLOC(d->ns->zeta_interp,
              gsl_spline_alloc(gsl_interp_cspline, NOISE_ZETAINTERP_N+1));
    d->ns->Nzeta_accel = d->Ncores;
    SAFEALLOC(d->ns->zeta_accel,
              malloc(d->ns->Nzeta_accel * sizeof(gsl_interp_accel *)));
    for (int ii=0; ii<d->ns->Nzeta_accel; ii++)
    {
        SAFEALLOC(d->ns->zeta_accel[ii], gsl_interp_accel_alloc());
    }
    SAFEGSL(gsl_spline_init(d->ns->zeta_interp, phi, zeta, NOISE_ZETAINTERP_N+1));

    free(phi);
    free(zeta);

    d->ns->created_zeta_interp = 1;
    
    ENDFCT
}//}}}

typedef struct//{{{
{
    int status;
    hmpdf_obj *d;
}//}}}
sigmasq_integrand_params;

static double
sigmasq_integrand(double ell, void *params)
// if LOGELL is defined, ell = log(ell)
{//{{{
    sigmasq_integrand_params *p = (sigmasq_integrand_params *)params;
    #ifdef LOGELL
    ell = exp(ell);
    #endif

    // evaluate the noise power spectrum
    double temp = ell * p->d->ns->noise_pwr(ell, p->d->ns->noise_pwr_params);

    // take care of Jacobian if necessary
    #ifdef LOGELL
    temp *= ell;
    #endif

    // apply the other filters
    p->status = apply_filters(p->d, 1, &ell, &temp, &temp, 1, filter_ps, NULL);

    return temp;
}//}}}

static int
create_noise_sigmasq(hmpdf_obj *d)
{//{{{
    STARTFCT
    
    sigmasq_integrand_params p;
    p.status = 0;
    p.d = d;

    gsl_function integrand;
    integrand.function = &sigmasq_integrand;
    integrand.params = &p;
    double err;
    #ifdef LOGELL
    double ellmin = log(NOISE_ELLMIN);
    double ellmax = log(NOISE_ELLMAX);
    #else
    double ellmin = NOISE_ELLMIN;
    double ellmax = NOISE_ELLMAX;
    #endif

    gsl_integration_workspace *ws;
    SAFEALLOC(ws, gsl_integration_workspace_alloc(NOISE_LIMIT));
    SAFEGSL(gsl_integration_qag(&integrand, ellmin, ellmax,
                                NOISE_EPSABS*(d->n->signalgrid[1]-d->n->signalgrid[0]),
                                NOISE_EPSREL, NOISE_LIMIT, NOISE_KEY,
                                ws, &(d->ns->sigmasq), &err));
    gsl_integration_workspace_free(ws);

    HMPDFCHECK(p.status, "error encountered during integration.");

    // correct normalization
    d->ns->sigmasq *= 0.5 * M_1_PI;

    ENDFCT
}//}}}

static int
create_toepl(hmpdf_obj *d)
// creates Toeplitz matrix nrows=Nsignal, ncols=Nsignal_noisy
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\tcreate_toepl\n");

    SAFEALLOC(d->ns->toepl, malloc(d->n->Nsignal * d->n->Nsignal_noisy
                                   * sizeof(double)));
    zero_real(d->n->Nsignal*d->n->Nsignal_noisy, d->ns->toepl);
    // dimensionless version of the variance
    double var = d->ns->sigmasq
                 / gsl_pow_2(d->n->signalgrid[1] - d->n->signalgrid[0]);
    // fill the first row
    for (int ii= -(int)d->ns->len_kernel; ii<=(int)d->ns->len_kernel; ii++)
    {
        d->ns->toepl[ii+d->ns->len_kernel] = exp(-0.5 * gsl_pow_2((double)(ii)) / var)
                                             /sqrt(2.0 * M_PI * var);
    }
    // fill the remaining rows
    for (long ii=1; ii<d->n->Nsignal; ii++)
    {
        memcpy(d->ns->toepl + ii*(d->n->Nsignal_noisy+1), d->ns->toepl,
               (2*d->ns->len_kernel+1) * sizeof(double));
    }

    ENDFCT
}//}}}

int
noise_vect(hmpdf_obj *d, double *in, double *out)
// assumes len(in) = Nsignal, len(out) = Nsignal_noisy
{//{{{
    STARTFCT

    HMPDFCHECK(d->ns->toepl == NULL, "Toeplitz matrix not computed.");

    cblas_dgemv(CblasRowMajor, CblasTrans/*toepl matrix needs to be transposed*/,
                d->n->Nsignal/*rows*/, d->n->Nsignal_noisy/*cols*/,
                1.0/*alpha*/, d->ns->toepl/*matrix A*/, d->n->Nsignal_noisy/*lda*/,
                in/*input vector X*/, 1/*stride of X*/, 0.0/*beta*/,
                out/*output vector Y*/, 1/*stride of Y*/);
    
    ENDFCT
}//}}}

// TODO
int
noise_matr(hmpdf_obj *d, double *in, double *out)
// assumes [in] = Nsignal*Nsignal, [out] = Nsignal_noisy*Nsignal_noisy
{//{{{
    STARTFCT

    HMPDFCHECK(d->ns->toepl == NULL, "Toeplitz matrix not computed.");

    // no aliasing allowed, so we need intermediate storage
    double *temp;
    SAFEALLOC(temp, malloc(d->n->Nsignal * d->n->Nsignal_noisy
                                     * sizeof(double)));

    // multiply from the left with Toeplitz matrix
    cblas_dgemm(CblasRowMajor, CblasTrans/*toepl matrix needs to be transposed*/,
                CblasNoTrans/*the right matrix is not transposed*/,
                d->n->Nsignal_noisy/*rows of left matrix & output*/,
                d->n->Nsignal/*cols of right matrix & output*/,
                d->n->Nsignal/*cols of left matrix, rows of right matrix*/,
                1.0/*alpha*/, d->ns->toepl/*left matrix*/,
                d->n->Nsignal_noisy/*lda*/, in/*right matrix*/,
                d->n->Nsignal/*ldb*/, 0.0/*beta*/, temp/*output matrix*/,
                d->n->Nsignal/*ldc*/);

    // multiply from the right with Toeplitz matrix (transposed)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                d->n->Nsignal_noisy, d->n->Nsignal_noisy, d->n->Nsignal,
                1.0, temp, d->n->Nsignal, d->ns->toepl, d->n->Nsignal_noisy,
                0.0, out, d->n->Nsignal_noisy);

    free(temp);

    ENDFCT
}//}}}

int
init_noise(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->ns->noise_pwr != NULL)
    {
        HMPDFPRINT(1, "init_noise\n");

        SAFEHMPDF(create_noise_sigmasq(d));
        SAFEHMPDF(create_noisy_grids(d));
        SAFEHMPDF(create_toepl(d));
    }

    ENDFCT
}//}}}
