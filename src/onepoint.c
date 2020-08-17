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
#include <gsl/gsl_integration.h>

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

    HMPDFPRINT(2, "\treset_onepoint\n");

    if (d->op->PDFu != NULL) { fftw_free(d->op->PDFu); }
    if (d->op->PDFc != NULL) { fftw_free(d->op->PDFc); }
    if (d->op->PDFu_noisy != NULL) { free(d->op->PDFu_noisy); }
    if (d->op->PDFc_noisy != NULL) { free(d->op->PDFc_noisy); }

    ENDFCT
}//}}}

int
correct_phase1d(hmpdf_obj *d, double complex *x, int stride, int sgn)
// sign is +1 for application after real -> double complex FFT
// and -1 for application before double complex -> real FFT
{//{{{
    STARTFCT

    if (d->n->Nsignal_negative > 0)
    {
        for (int ii=0; ii<d->n->Nsignal/2+1; ii++)
        {
            x[ii*stride]
                *= cexp(- (double complex)sgn * _Complex_I * d->n->signalmin * d->n->lambdagrid[ii]);
        }
    }

    ENDFCT
}//}}}

static int
op_segmentsum(hmpdf_obj *d, int z_index, int M_index, double *au, double *ac)
{//{{{
    STARTFCT

    double n = d->h->hmf[z_index][M_index];
    double b = d->h->bias[z_index][M_index];
    for (int segment=0;
         segment<d->p->segment_boundaries[z_index][M_index][0];
         segment++)
    {
        batch_t bt;
        SAFEHMPDF(inv_profile(d, z_index, M_index, segment, dtsq_of_s, &bt));
        for (int signalindex=bt.start, ii=0;
             ii < bt.len;
             signalindex += bt.incr, ii++)
        {
            au[signalindex] += fabs(bt.data[ii]) * M_PI * n
                               * d->n->Mweights[M_index];
            ac[signalindex] += fabs(bt.data[ii]) * M_PI * n * b
                               * d->n->Mweights[M_index];
        }
        delete_batch(&bt);
    }

    ENDFCT
}//}}}

static int
op_Mint(hmpdf_obj *d, int z_index, double *au, double *ac)
// performs the mass integrals, for the invertible profiles
// integrals are _added_ to au, ac
{//{{{
    STARTFCT

    for (int M_index=0; M_index<d->n->NM; M_index++)
    {
        SAFEHMPDF(op_segmentsum(d, z_index, M_index, au, ac));
    }

    ENDFCT
}//}}}

static int
op_zint(hmpdf_obj *d, double complex *pu_comp, double complex *pc_comp) // p is the exponent in P(lambda)
{//{{{
    STARTFCT

    double *ac_real;
    SAFEALLOC(ac_real, fftw_malloc((d->n->Nsignal+2) * sizeof(double)));
    double complex *ac_comp = (double complex *)ac_real;
    double *au_real;
    SAFEALLOC(au_real, fftw_malloc((d->n->Nsignal+2) * sizeof(double)));
    double complex *au_comp = (double complex *)au_real;
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

        SAFEHMPDF(op_Mint(d, z_index, au_real, ac_real));
        // perform FFTs real -> double complex
        fftw_execute(plan_u);
        fftw_execute(plan_c);
        // correct phases
        SAFEHMPDF(correct_phase1d(d, au_comp, 1, 1));
        SAFEHMPDF(correct_phase1d(d, ac_comp, 1, 1));

        for (int ii=0; ii<d->n->Nsignal/2+1; ii++)
        {
            // subtract the zero modes, square the clustered mass integral
            double complex tempu = au_comp[ii] - au_comp[0];
            double complex tempc = (ac_comp[ii] - ac_comp[0])
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
_mean(int N, const double *const x, const double *const p, double *out)
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
                    d->op->PDFu, &(d->op->signalmeanu)));
    SAFEHMPDF(_mean(d->n->Nsignal, d->n->signalgrid,
                    d->op->PDFc, &(d->op->signalmeanc)));

    ENDFCT
}//}}}

int
create_noisy_op(hmpdf_obj *d)
// convolves the original PDF with a Gaussian kernel of width sigma = noise
{//{{{
    STARTFCT

    if (d->op->created_noisy_op) { return hmpdf_status; }

    HMPDFPRINT(2, "\tcreate_noisy_op\n");

    double *in[] = {d->op->PDFu, d->op->PDFc};
    double **out[] = {&d->op->PDFu_noisy, &d->op->PDFc_noisy};
    for (int ii=0; ii<2; ii++)
    {
        SAFEALLOC(*out[ii], malloc(d->n->Nsignal_noisy * sizeof(double)));
        SAFEHMPDF(noise_vect(d, in[ii], *out[ii]));
    }

    d->op->created_noisy_op = 1;

    ENDFCT
}//}}}

int
create_op(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->op->created_op) { return hmpdf_status; }

    HMPDFPRINT(2, "\tcreate_op\n");
    
    SAFEALLOC(d->op->PDFu, fftw_malloc((d->n->Nsignal + 2) * sizeof(double)));
    SAFEALLOC(d->op->PDFc, fftw_malloc((d->n->Nsignal + 2) * sizeof(double)));
    double complex *PDFu_comp = (double complex *)d->op->PDFu;
    double complex *PDFc_comp = (double complex *)d->op->PDFc;

    fftw_plan plan_u = fftw_plan_dft_c2r_1d(d->n->Nsignal, PDFu_comp, d->op->PDFu, FFTW_ESTIMATE);
    fftw_plan plan_c = fftw_plan_dft_c2r_1d(d->n->Nsignal, PDFc_comp, d->op->PDFc, FFTW_ESTIMATE);

    // perform redshift integration
    SAFEHMPDF(op_zint(d, PDFu_comp, PDFc_comp));

    // take exponential and normalize
    for (int ii=0; ii<d->n->Nsignal/2+1; ii++)
    {
        PDFu_comp[ii] = cexp(PDFu_comp[ii])/(double)(d->n->Nsignal);
        PDFc_comp[ii] = cexp(PDFc_comp[ii])/(double)(d->n->Nsignal);
    }

    // correct phases
    SAFEHMPDF(correct_phase1d(d, PDFu_comp, 1, -1));
    SAFEHMPDF(correct_phase1d(d, PDFc_comp, 1, -1));
    // transform back to real space
    fftw_execute(plan_u);
    fftw_execute(plan_c);
    fftw_destroy_plan(plan_u);
    fftw_destroy_plan(plan_c);

    // compute the mean of the distributions
    SAFEHMPDF(get_mean_signal(d));

    d->op->created_op = 1;

    ENDFCT
}//}}}

static int
prepare_op(hmpdf_obj *d)
// does the necessary create calls
{//{{{
    STARTFCT

    HMPDFPRINT(1, "prepare_op\n");

    if (d->f->Nfilters > 0)
    {
        SAFEHMPDF(create_conj_profiles(d));
        SAFEHMPDF(create_filtered_profiles(d));
    }
    SAFEHMPDF(create_segments(d));
    SAFEHMPDF(create_op(d));
    if (d->ns->noise > 0.0)
    // include gaussian noise
    {
        SAFEHMPDF(create_noisy_op(d));
    }

    ENDFCT
}//}}}

int
pdf_check_user_input(hmpdf_obj *d, int Nbins, double binedges[Nbins+1], int noisy)
{//{{{
    STARTFCT

    HMPDFCHECK(noisy && d->ns->noise<0.0,
               "noisy output requested but no/invalid noise level passed.");

    HMPDFCHECK(not_monotonic(Nbins+1, binedges, 1),
               "binedges not monotonically increasing.");

    ENDFCT
}//}}}

int
hmpdf_get_op(hmpdf_obj *d, int Nbins, double binedges[Nbins+1], double op[Nbins], int incl_2h, int noisy)
{//{{{
    STARTFCT

    SAFEHMPDF(pdf_check_user_input(d, Nbins, binedges, noisy));

    SAFEHMPDF(prepare_op(d));
    
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
                     (noisy) ? ((incl_2h) ? d->op->PDFc_noisy : d->op->PDFu_noisy)
                     : ((incl_2h) ? d->op->PDFc : d->op->PDFu),
                     Nbins, _binedges, op, OPINTERP_TYPE));

    ENDFCT
}//}}}

