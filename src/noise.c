#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include <fftw3.h>

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

    d->ns->created_noise_zeta_interp = 0;
    d->ns->zeta_interp = NULL;
    d->ns->zeta_accel = NULL;

    d->ns->conv_buffer_real = NULL;
    d->ns->pconv_r2c = NULL;
    d->ns->pconv_c2r = NULL;

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
        for (int ii=0; ii<d->Ncores; ii++)
        {
            gsl_interp_accel_free(d->ns->zeta_accel[ii]);
        }
        free(d->ns->zeta_accel);
    }
    if (d->ns->conv_buffer_real != NULL)
    {
        for (int ii=0; ii<d->Ncores; ii++)
        {
            if (d->ns->conv_buffer_real[ii] != NULL)
            {
                fftw_free(d->ns->conv_buffer_real[ii]);
            }
        }
        free(d->ns->conv_buffer_real);
    }
    if (d->ns->pconv_r2c != NULL)
    {
        for (int ii=0; ii<d->Ncores; ii++)
        {
            if (d->ns->pconv_r2c[ii] != NULL)
            {
                fftw_destroy_plan(*(d->ns->pconv_r2c[ii]));
                free(d->ns->pconv_r2c[ii]);
            }
        }
        free(d->ns->pconv_r2c);
    }
    if (d->ns->pconv_c2r != NULL)
    {
        for (int ii=0; ii<d->Ncores; ii++)
        {
            if (d->ns->pconv_c2r[ii] != NULL)
            {
                fftw_destroy_plan(*(d->ns->pconv_c2r[ii]));
                free(d->ns->pconv_c2r[ii]);
            }
        }
        free(d->ns->pconv_c2r);
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
    SAFEALLOC(d->n->lambdagrid_noisy, malloc((d->n->Nsignal_noisy/2+1)
                                             * sizeof(double)));
    SAFEHMPDF(linspace(d->n->Nsignal_noisy/2+1,
                       0.0,
                       M_PI / (d->n->signalgrid_noisy[1]-d->n->signalgrid_noisy[0]),
                       d->n->lambdagrid_noisy));

    ENDFCT
}//}}}

int
create_noise_zeta_interp(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->ns->created_noise_zeta_interp) { return 0; }

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
    SAFEALLOC(d->ns->zeta_accel,
              malloc(d->Ncores * sizeof(gsl_interp_accel *)));
    for (int ii=0; ii<d->Ncores; ii++)
    {
        SAFEALLOC(d->ns->zeta_accel[ii], gsl_interp_accel_alloc());
    }
    SAFEGSL(gsl_spline_init(d->ns->zeta_interp, phi, zeta, NOISE_ZETAINTERP_N+1));

    // TODO debugging
    // savetxt("noise_corrfunc.dat", NOISE_ZETAINTERP_N+1, 2, phi, zeta);

    free(phi);
    free(zeta);

    d->ns->created_noise_zeta_interp = 1;
    
    ENDFCT
}//}}}

typedef struct
{//{{{
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
    double out = ell * p->d->ns->noise_pwr(ell, p->d->ns->noise_pwr_params);

    // take care of Jacobian if necessary
    #ifdef LOGELL
    out *= ell;
    #endif

    // apply the other filters
    p->status = apply_filters(p->d, 1, &ell, &out, &out, 1, filter_ps, NULL);

    return out;
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
                                NOISE_EPSABS
                                    * gsl_pow_2(d->n->signalgrid[1]-d->n->signalgrid[0]),
                                NOISE_EPSREL, NOISE_LIMIT, NOISE_KEY,
                                ws, &(d->ns->sigmasq), &err));
    gsl_integration_workspace_free(ws);

    HMPDFCHECK(p.status, "error encountered during integration.");

    // correct normalization
    d->ns->sigmasq *= 0.5 * M_1_PI;

    // TODO debugging
    printf("sigma = %.8e\n", sqrt(d->ns->sigmasq));

    ENDFCT
}//}}}

static inline int
safe_exp_div(double exponent, double divisor, double *out)
// computes exp(-exponent)/divisor in an underflow-safe way
{//{{{
    STARTFCT

    if (exponent > - GSL_MAX(0.0, log(divisor))
                   - log((double)FLT_RADIX) * (double)DBL_MIN_EXP
                   - 2.0)
    {
        *out = 0.0;
    }
    else
    {
        *out = exp(-exponent)/divisor;
    }

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
    for (long ii= -d->ns->len_kernel; ii<=d->ns->len_kernel; ii++)
    {
        SAFEHMPDF(safe_exp_div(0.5 * gsl_pow_2((double)(ii)) / var,
                               sqrt(2.0 * M_PI * var),
                               d->ns->toepl + ii + d->ns->len_kernel));
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
create_noise_matr_conv(hmpdf_obj *d, int Nbuffers)
{//{{{
    STARTFCT

    HMPDFCHECK(Nbuffers > d->Ncores,
               "too many buffers requested.");

    SAFEHMPDF(create_noise_zeta_interp(d));

    if (d->ns->conv_buffer_real == NULL)
    {
        SAFEALLOC(d->ns->conv_buffer_real,
                  malloc(d->Ncores * sizeof(double *)));
        SAFEALLOC(d->ns->conv_buffer_comp,
                  malloc(d->Ncores * sizeof(double complex *)));
        for (int ii=0; ii<d->Ncores; ii++)
        {
            d->ns->conv_buffer_real[ii] = NULL;
        }
    }

    if (d->ns->pconv_r2c == NULL)
    {
        SAFEALLOC(d->ns->pconv_r2c,
                  malloc(d->Ncores * sizeof(fftw_plan *)));
        SAFEALLOC(d->ns->pconv_c2r,
                  malloc(d->Ncores * sizeof(fftw_plan *)));
        for (int ii=0; ii<d->Ncores; ii++)
        {
            d->ns->pconv_r2c[ii] = NULL;
            d->ns->pconv_c2r[ii] = NULL;
        }
    }

    for (int ii=0; ii<Nbuffers; ii++)
    {
        if (d->ns->conv_buffer_real[ii] == NULL)
        {
            SAFEALLOC(d->ns->conv_buffer_real[ii],
                      fftw_malloc((d->n->Nsignal_noisy+2)
                                  * d->n->Nsignal_noisy
                                  * sizeof(double)));
            d->ns->conv_buffer_comp[ii]
                = (double complex *)d->ns->conv_buffer_real[ii];
        }

        if (d->ns->pconv_r2c[ii] == NULL)
        {
            SAFEALLOC(d->ns->pconv_r2c[ii], malloc(sizeof(fftw_plan)));
            *(d->ns->pconv_r2c[ii])
                = fftw_plan_dft_r2c_2d(d->n->Nsignal_noisy, d->n->Nsignal_noisy,
                                       d->ns->conv_buffer_real[ii],
                                       d->ns->conv_buffer_comp[ii],
                                       FFTW_MEASURE);
            SAFEALLOC(d->ns->pconv_c2r[ii], malloc(sizeof(fftw_plan)));
            *(d->ns->pconv_c2r[ii])
                = fftw_plan_dft_c2r_2d(d->n->Nsignal_noisy, d->n->Nsignal_noisy,
                                       d->ns->conv_buffer_comp[ii],
                                       d->ns->conv_buffer_real[ii],
                                       FFTW_MEASURE);
        }
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

static int
multiply_w_gaussian2d(hmpdf_obj *d, double phi, double complex *A)
// note : this fct includes the FFT normalization
{//{{{
    STARTFCT

    // evaluate pixel-pixel noise correlation function
    double zeta;
    SAFEGSL(gsl_spline_eval_e(d->ns->zeta_interp, phi,
                              d->ns->zeta_accel[this_core()],
                              &zeta));

    // multiply with the Fourier space noise kernel
    for (long ii=0; ii<d->n->Nsignal_noisy; ii++)
    // loop over long direction (rows)
    {
        // same logic as in the "redundant" function,
        //      note that the FT of a real Gaussian is real valued
        //      so conjugation not required
        double lambda1;
        if (ii<=d->n->Nsignal_noisy/2)
        {
            lambda1 = d->n->lambdagrid_noisy[ii];
        }
        else
        {
            lambda1 = d->n->lambdagrid_noisy[d->n->Nsignal_noisy-ii];
        }

        for (long jj=0; jj<d->n->Nsignal_noisy/2+1; jj++)
        // loop over short direction (cols)
        {
            double lambda2 = d->n->lambdagrid_noisy[jj];

            // evaluate the Gaussian kernel in Fourier space
            //     and normalize properly
            double temp = 0.0;
            SAFEHMPDF(safe_exp_div(0.5 * d->ns->sigmasq * (gsl_pow_2(lambda1)+gsl_pow_2(lambda2))
                                   + zeta * lambda1 * lambda2,
                                   gsl_pow_2((double)d->n->Nsignal_noisy),
                                   &temp));
            A[ii*(d->n->Nsignal_noisy/2+1) + jj]
                *= temp;
        }
    }

    ENDFCT
}//}}}

int
noise_matr(hmpdf_obj *d, double *in, double *out, int is_buffered, double phi)
// assumes [in] = Nsignal*Nsignal if !is_buffered, else (Nsignal+2)*Nsignal,
//         [out] = Nsignal_noisy*Nsignal_noisy or NULL
// CAUTION: this function is not thread safe!
//              (we need too much buffer space for that to make sense)
//          Need to consider this in any OMP environment it is called from
{//{{{
    STARTFCT

    HMPDFCHECK(d->ns->conv_buffer_real[this_core()] == NULL,
               "trying to access uninitialized buffer.");

    zero_real((d->n->Nsignal_noisy+2)*d->n->Nsignal_noisy,
              d->ns->conv_buffer_real[this_core()]);
    // copy into the buffer, which is zero padded by d->ns->len_kernel in each direction
    for (long ii=0; ii<d->n->Nsignal; ii++)
    {
        memcpy(d->ns->conv_buffer_real[this_core()]
                   + (d->ns->len_kernel+ii) * (d->n->Nsignal_noisy+2) // go downwards
                   + d->ns->len_kernel, // go right
               in
                   + ii * ((is_buffered) ? (d->n->Nsignal+2) : d->n->Nsignal),
               d->n->Nsignal * sizeof(double));
    }

    HMPDFCHECK(d->ns->pconv_r2c[this_core()] == NULL,
               "trying to use uninitialized plan.");
    
    // perform the forward FFT
    fftw_execute(*(d->ns->pconv_r2c[this_core()]));

    // apply the filter in Fourier space
    SAFEHMPDF(multiply_w_gaussian2d(d, phi, d->ns->conv_buffer_comp[this_core()]));

    HMPDFCHECK(d->ns->pconv_c2r[this_core()] == NULL,
               "trying to use uninitialized plan.");

    // transform back to real space
    fftw_execute(*(d->ns->pconv_c2r[this_core()]));

    if (out != NULL)
    {
        // copy into the output array (which is not buffered)
        //    and normalize properly
        for (long ii=0; ii<d->n->Nsignal_noisy; ii++)
        {
            memcpy(out
                       + ii * d->n->Nsignal_noisy,
                   d->ns->conv_buffer_real[this_core()]
                       + ii * (d->n->Nsignal_noisy+2),
                   d->n->Nsignal_noisy * sizeof(double));
        }
    }

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
