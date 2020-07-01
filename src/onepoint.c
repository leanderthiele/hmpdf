#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <gsl/gsl_math.h>

#include "utils.h"
#include "configs.h"
#include "object.h"
#include "profiles.h"
#include "noise.h"
#include "numerics.h"
#include "onepoint.h"

#include "hmpdf.h"

int
null_onepoint(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->op->created_op = 0;
    d->op->created_noisy_op = 0;
    d->op->PDFu = NULL;
    d->op->PDFc = NULL;
    d->op->PDFu_noisy = NULL;
    d->op->PDFc_noisy = NULL;

    ENDFCT
}//}}}

int
reset_onepoint(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\treset_onepoint\n")

    if (d->op->PDFu != NULL) { free(d->op->PDFu); }
    if (d->op->PDFc != NULL) { free(d->op->PDFc); }
    if (d->op->PDFu_noisy != NULL) { free(d->op->PDFu_noisy); }
    if (d->op->PDFc_noisy != NULL) { free(d->op->PDFc_noisy); }

    ENDFCT
}//}}}

static int
op_Mint_invertible(hmpdf_obj *d, int z_index, double *au, double *ac)
// performs the mass integrals, for the invertible profiles
// integrals are _added_ to au, ac
{//{{{
    STARTFCT

    SAFEALLOC(double *, temp, malloc(d->n->Nsignal * sizeof(double)))

    for (int M_index=0; M_index<d->n->NM; M_index++)
    {
        if (d->p->is_not_monotonic[z_index][M_index])
        {
            continue;
        }
        SAFEHMPDF(dtsq_of_s(d, z_index, M_index, temp))
        double n = d->h->hmf[z_index][M_index];
        double b = d->h->bias[z_index][M_index];
        for (int ii=0; ii<d->n->Nsignal; ii++)
        {
            au[ii] -= temp[ii] * M_PI * n
                      * d->n->Mweights[M_index];
            ac[ii] -= temp[ii] * M_PI * n * b
                      * d->n->Mweights[M_index];
        }
    }

    free(temp);

    ENDFCT
}//}}}

static int
not_inv_integral(hmpdf_obj *d, double *s, int lambda_index, complex *out)
// TODO can potentially speed this up a bit
{//{{{
    STARTFCT

    SAFEALLOC(complex *, integr, malloc(d->p->prtilde_Ntheta * sizeof(complex)))
    // fill the integrand
    for (int ii=0; ii<d->p->prtilde_Ntheta; ii++)
    {
        integr[ii] = d->p->prtilde_thetagrid[ii]
                     * cexp(- _Complex_I * s[ii] * d->n->lambdagrid[lambda_index]);
    }
    // perform the integration
    *out = integr_comp(d->p->prtilde_Ntheta, 1.0/(double)(d->p->prtilde_Ntheta-1),
                       1, integr);
    free(integr);

    ENDFCT
}//}}}

static int
lambda_loop(hmpdf_obj *d, int z_index, int M_index, complex *out)
{//{{{
    STARTFCT

    SAFEALLOC(double *, s, malloc(d->p->prtilde_Ntheta * sizeof(double)))
    // evaluate the interpolated signal profile
    SAFEHMPDF(s_of_t(d, z_index, M_index, d->p->prtilde_Ntheta,
                     d->p->prtilde_thetagrid, s))

    // loop over lambda values
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(d->Ncores)
    #endif
    for (int ii=0; ii<d->n->Nsignal/2+1; ii++)
    {
        if (hmpdf_status) { continue; }
        SAFEHMPDF_NORETURN(not_inv_integral(d, s, ii, out+ii))
        if (hmpdf_status) { continue; }
        out[ii] *= 2.0 * M_PI * gsl_pow_2(d->p->profiles[z_index][M_index][0]);
    }

    free(s);

    ENDFCT
}//}}}

static int
op_Mint_notinvertible(hmpdf_obj *d, int z_index, complex *au, complex *ac)
// performs the mass integration over the non-invertible profiles
// integrals are _added_ to au, ac
{//{{{
    STARTFCT

    SAFEALLOC(complex *, temp, malloc((d->n->Nsignal/2+1) * sizeof(complex)))
    for (int M_index=0; M_index<d->n->NM; M_index++)
    {
        if (!(d->p->is_not_monotonic[z_index][M_index]))
        {
            continue;
        }
        SAFEHMPDF(lambda_loop(d, z_index, M_index, temp))
        double n = d->h->hmf[z_index][M_index];
        double b = d->h->bias[z_index][M_index];
        for (int ii=0; ii<d->n->Nsignal/2+1; ii++)
        {
            au[ii] += temp[ii] * n
                      * d->n->Mweights[M_index];
            ac[ii] += temp[ii] * n * b
                      * d->n->Mweights[M_index];
        }
    }

    free(temp);

    ENDFCT
}//}}}

static int
op_zint(hmpdf_obj *d, complex *pu_comp, complex *pc_comp) // p is the exponent in P(lambda)
{//{{{
    STARTFCT

    SAFEALLOC(double *, ac_real, fftw_malloc((d->n->Nsignal+2) * sizeof(double)))
    complex *ac_comp = (complex *)ac_real;
    SAFEALLOC(double *, au_real, fftw_malloc((d->n->Nsignal+2) * sizeof(double)))
    complex *au_comp = (complex *)au_real;
    fftw_plan plan_u = fftw_plan_dft_r2c_1d(d->n->Nsignal, au_real, au_comp, FFTW_MEASURE);
    fftw_plan plan_c = fftw_plan_dft_r2c_1d(d->n->Nsignal, ac_real, ac_comp, FFTW_MEASURE);

    // zero the output arrays
    for (int ii=0; ii<d->n->Nsignal/2+1; ii++)
    {
        pu_comp[ii] = pc_comp[ii] = 0.0;
    }

    // fill the z-integrand
    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        // zero the arrays
        for (int ii=0; ii<d->n->Nsignal/2+1; ii++)
        {
            au_comp[ii] = ac_comp[ii] = 0.0;
        }

        if (d->p->is_not_monotonic[z_index][d->n->NM] < d->n->NM)
        // there are samples for which we can do the FFT integral
        {
            SAFEHMPDF(op_Mint_invertible(d, z_index, au_real, ac_real))
            fftw_execute(plan_u);
            fftw_execute(plan_c);
        }
        if (d->p->is_not_monotonic[z_index][d->n->NM])
        // there are samples for which we need to do the slow integral
        {
            SAFEHMPDF(op_Mint_notinvertible(d, z_index, au_comp, ac_comp))
        }

        for (int ii=0; ii<d->n->Nsignal/2+1; ii++)
        {
            // subtract the zero modes, square the clustered mass integral
            complex tempu = au_comp[ii] - au_comp[0];
            complex tempc = (ac_comp[ii] - ac_comp[0])
                            * (ac_comp[ii] - ac_comp[0]);
            // multiply with the prefactors
            tempu *= gsl_pow_2(d->c->comoving[z_index])
                     / d->c->hubble[z_index];
            tempc *= 0.5 * d->c->Dsq[z_index] * d->pwr->autocorr
                     * gsl_pow_4(d->c->comoving[z_index])
                     / d->c->hubble[z_index];
            tempc += tempu;

            pu_comp[ii] += tempu * d->n->zweights[z_index];
            pc_comp[ii] += tempc * d->n->zweights[z_index];
        }
    }
    fftw_free(au_real);
    fftw_free(ac_real);
    fftw_destroy_plan(plan_u);
    fftw_destroy_plan(plan_c);

    ENDFCT
}//}}}

static int
_mean(int N, double *x, double *p, double *out)
{//{{{
    STARTFCT

    *out = 0.0;
    double norm = 0.0;
    for (int ii=0; ii<N; ii++)
    {
        *out += p[ii] * x[ii];
        norm += p[ii];
    }
    *out /= norm;

    ENDFCT
}//}}}

static int
get_mean_signal(hmpdf_obj *d)
{//{{{
    STARTFCT

    SAFEHMPDF(_mean(d->n->Nsignal, d->n->signalgrid,
                    d->op->PDFu, &(d->op->signalmeanu)))
    SAFEHMPDF(_mean(d->n->Nsignal, d->n->signalgrid,
                    d->op->PDFc, &(d->op->signalmeanc)))

    ENDFCT
}//}}}

int
create_noisy_op(hmpdf_obj *d)
// convolves the original PDF with a Gaussian kernel of width sigma = noise
{//{{{
    STARTFCT

    if (d->op->created_noisy_op) { return hmpdf_status; }

    HMPDFPRINT(2, "\tcreate_noisy_op\n")

    double *in[] = {d->op->PDFu, d->op->PDFc};
    double **out[] = {&d->op->PDFu_noisy, &d->op->PDFc_noisy};
    for (int ii=0; ii<2; ii++)
    {
        SAFEALLOC(, *out[ii], malloc(d->n->Nsignal_noisy * sizeof(double)))
        SAFEHMPDF(noise_vect(d, in[ii], *out[ii]))
    }

    d->op->created_noisy_op = 1;

    ENDFCT
}//}}}

int
create_op(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->op->created_op) { return hmpdf_status; }

    HMPDFPRINT(2, "\tcreate_op\n")
    
    SAFEALLOC(, d->op->PDFu, fftw_malloc((d->n->Nsignal + 2) * sizeof(double)))
    SAFEALLOC(, d->op->PDFc, fftw_malloc((d->n->Nsignal + 2) * sizeof(double)))
    complex *PDFu_comp = (complex *)d->op->PDFu;
    complex *PDFc_comp = (complex *)d->op->PDFc;

    fftw_plan plan_u = fftw_plan_dft_c2r_1d(d->n->Nsignal, PDFu_comp, d->op->PDFu, FFTW_ESTIMATE);
    fftw_plan plan_c = fftw_plan_dft_c2r_1d(d->n->Nsignal, PDFc_comp, d->op->PDFc, FFTW_ESTIMATE);

    // perform redshift integration
    SAFEHMPDF(op_zint(d, PDFu_comp, PDFc_comp))

    // take exponential and normalize
    for (int ii=0; ii<d->n->Nsignal/2+1; ii++)
    {
        PDFu_comp[ii] = cexp(PDFu_comp[ii])/(double)(d->n->Nsignal);
        PDFc_comp[ii] = cexp(PDFc_comp[ii])/(double)(d->n->Nsignal);
    }

    // transform back to real space
    fftw_execute(plan_u);
    fftw_execute(plan_c);
    fftw_destroy_plan(plan_u);
    fftw_destroy_plan(plan_c);

    // compute the mean of the distributions if needed later
    SAFEHMPDF(get_mean_signal(d))

    d->op->created_op = 1;

    ENDFCT
}//}}}

static int
prepare_op(hmpdf_obj *d)
// does the necessary create calls
{//{{{
    STARTFCT

    HMPDFPRINT(1, "prepare_op\n")

    if (d->f->Nfilters > 0)
    {
        SAFEHMPDF(create_conj_profiles(d))
        SAFEHMPDF(create_filtered_profiles(d))
    }
    SAFEHMPDF(create_monotonicity(d))
    SAFEHMPDF(create_op(d))
    if (d->ns->noise > 0.0)
    // include gaussian noise
    {
        SAFEHMPDF(create_noisy_op(d))
    }
    
    ENDFCT
}//}}}

int
hmpdf_get_op(hmpdf_obj *d, int Nbins, double binedges[Nbins+1], double op[Nbins], int incl_2h, int noisy)
{//{{{
    STARTFCT

    if (noisy && d->ns->noise<0.0)
    {
        HMPDFERR("noisy pdf requested but no/invalid noise level passed.")
    }

    if (not_monotonic(Nbins+1, binedges, NULL))
    {
        HMPDFERR("binedges not monotonically increasing.")
    }

    SAFEHMPDF(prepare_op(d))
    
    double _binedges[Nbins+1];
    memcpy(_binedges, binedges, (Nbins+1) * sizeof(double));
    if (d->p->stype == hmpdf_kappa)
    {
        for (int ii=0; ii<=Nbins; ii++)
        {
            _binedges[ii] += (incl_2h) ? d->op->signalmeanc : d->op->signalmeanu;
        }
    }

    SAFEHMPDF(bin_1d((noisy) ? d->n->Nsignal_noisy : d->n->Nsignal,
                     (noisy) ? d->n->signalgrid_noisy : d->n->signalgrid,
                     (noisy) ?  ((incl_2h) ? d->op->PDFc_noisy : d->op->PDFu_noisy)
                     : ((incl_2h) ? d->op->PDFc : d->op->PDFu),
                     Nbins, _binedges, op, OPINTERP_TYPE))

    ENDFCT
}//}}}

