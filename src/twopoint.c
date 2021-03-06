#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include <gsl/gsl_math.h>

#include "utils.h"
#include "configs.h"
#include "object.h"
#include "power.h"
#include "profiles.h"
#include "onepoint.h"

#include "hmpdf.h"

int
null_twopoint(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->tp->created_phi_indep = 0;
    d->tp->dtsq = NULL;
    d->tp->t = NULL;
    d->tp->ac = NULL;
    d->tp->au = NULL;
    d->tp->ws = NULL;
    d->tp->last_phi = -1.0;
    d->tp->pdf = NULL;
    d->tp->pdf_noisy = NULL;

    ENDFCT
}//}}}

int
reset_twopoint(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\treset_twopoint\n");

    if (d->tp->dtsq != NULL)
    {
        for (int z_index=0; z_index<d->n->Nz; z_index++)
        {
            if (d->tp->dtsq[z_index] != NULL)
            {
                for (int M_index=0; M_index<d->n->NM; M_index++)
                {
                    if (d->tp->dtsq[z_index][M_index] != NULL)
                    {
                        for (int segment=0;
                             segment<d->p->segment_boundaries[z_index][M_index][0];
                             segment++)
                        {
                            delete_batch(d->tp->dtsq[z_index][M_index]+segment);
                        }
                        free(d->tp->dtsq[z_index][M_index]);
                    }
                }
                free(d->tp->dtsq[z_index]);
            }
        }
        free(d->tp->dtsq);
    }
    if (d->tp->t != NULL)
    {
        for (int z_index=0; z_index<d->n->Nz; z_index++)
        {
            if (d->tp->t[z_index] != NULL)
            {
                for (int M_index=0; M_index<d->n->NM; M_index++)
                {
                    if (d->tp->t[z_index][M_index] != NULL)
                    {
                        for (int segment=0;
                             segment<d->p->segment_boundaries[z_index][M_index][0];
                             segment++)
                        {
                            delete_batch(d->tp->t[z_index][M_index]+segment);
                        }
                        free(d->tp->t[z_index][M_index]);
                    }
                }
                free(d->tp->t[z_index]);
            }
        }
        free(d->tp->t);
    }
    if (d->tp->ac != NULL)
    {
        for (int z_index=0; z_index<d->n->Nz; z_index++)
        {
            if (d->tp->ac[z_index] != NULL)
            {
                free(d->tp->ac[z_index]);
            }
        }
        free(d->tp->ac);
    }
    if (d->tp->au != NULL) { fftw_free(d->tp->au); }
    if (d->tp->ws != NULL)
    {
        fftw_free(d->tp->ws->pdf_real);
        free(d->tp->ws->bc);
        fftw_free(d->tp->ws->tempc_real);
        fftw_destroy_plan(d->tp->ws->pu_r2c);
        fftw_destroy_plan(d->tp->ws->pc_r2c);
        fftw_destroy_plan(d->tp->ws->ppdf_c2r);
        free(d->tp->ws);
    }
    if (d->tp->pdf != NULL) { free(d->tp->pdf); }
    if (d->tp->pdf_noisy != NULL) { free(d->tp->pdf_noisy); }

    ENDFCT
}//}}}

static int
correct_phase2d(hmpdf_obj *d, double complex *x, int sgn)
{//{{{
    STARTFCT

    if (d->n->Nsignal_negative > 0)
    {
        for (long ii=0; ii<d->n->Nsignal; ii++)
        // loop over long direction (rows)
        {
            double lambda1 = WAVENR(d->n->Nsignal, d->n->lambdagrid, ii);

            for (long jj=0; jj<d->n->Nsignal/2+1; jj++)
            // loop over short direction (cols)
            {
                double lambda2 = WAVENR(d->n->Nsignal, d->n->lambdagrid, jj);

                x[ii*(d->n->Nsignal/2+1)+jj]
                    *= cexp(- (double complex)sgn * _Complex_I
                              * d->n->signalmin * (lambda1 + lambda2));
            }
        }
    }

    ENDFCT
}//}}}

int
create_phi_indep(hmpdf_obj *d)
// computes tp->dtsq, tp->t, tp->ac
{//{{{
    STARTFCT

    if (d->tp->created_phi_indep) { return 0; }

    HMPDFPRINT(2, "\tcreate_phi_indep\n");
    
    SAFEALLOC(d->tp->dtsq, malloc(d->n->Nz * sizeof(batch_t **)));
    SETARRNULL(d->tp->dtsq, d->n->Nz);
    SAFEALLOC(d->tp->t,    malloc(d->n->Nz * sizeof(batch_t **)));
    SETARRNULL(d->tp->t,    d->n->Nz);
    SAFEALLOC(d->tp->ac,   malloc(d->n->Nz * sizeof(double complex *)));
    SETARRNULL(d->tp->ac,   d->n->Nz);
    SAFEALLOC(d->tp->au,   fftw_malloc((d->n->Nsignal/2+1) * sizeof(double complex)));
    double *au_real = (double *)(d->tp->au);

    double *tempc_real;
    SAFEALLOC(tempc_real, fftw_malloc((d->n->Nsignal+2)*sizeof(double)));
    double complex *tempc_comp = (double complex *)tempc_real;
    fftw_plan plan_c = fftw_plan_dft_r2c_1d(d->n->Nsignal, tempc_real,
                                            tempc_comp, FFTW_MEASURE);
    fftw_plan plan_u = fftw_plan_dft_r2c_1d(d->n->Nsignal, au_real,
                                            d->tp->au, FFTW_MEASURE);

    // zero the unclustered array
    zero_comp(d->n->Nsignal/2+1, d->tp->au);

    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        SAFEALLOC(d->tp->dtsq[z_index], malloc(d->n->NM * sizeof(batch_t *)));
        SETARRNULL(d->tp->dtsq[z_index], d->n->NM);
        SAFEALLOC(d->tp->t[z_index],    malloc(d->n->NM * sizeof(batch_t *)));
        SETARRNULL(d->tp->t[z_index],    d->n->NM);
        SAFEALLOC(d->tp->ac[z_index],   malloc((d->n->Nsignal/2+1) * sizeof(double complex)));

        // zero the FFT array
        zero_real(d->n->Nsignal+2, tempc_real);

        // integrate clustered contribution over mass;
        // and fill the theta(signal) dtheta^2/dsignal grids
        for (int M_index=0; M_index<d->n->NM; M_index++)
        {
            double n = d->h->hmf[z_index][M_index];
            double b = d->h->bias[z_index][M_index];

            SAFEALLOC(d->tp->dtsq[z_index][M_index],
                      malloc(d->p->segment_boundaries[z_index][M_index][0]
                             * sizeof(batch_t)));
            SAFEALLOC(d->tp->t[z_index][M_index],
                      malloc(d->p->segment_boundaries[z_index][M_index][0]
                             * sizeof(batch_t)));

            for (int segment=0;
                 segment<d->p->segment_boundaries[z_index][M_index][0];
                 segment++)
            {
                SAFEHMPDF(inv_profile(d, z_index, M_index, segment,
                                      dtsq_of_s, d->tp->dtsq[z_index][M_index]+segment));
                SAFEHMPDF(inv_profile(d, z_index, M_index, segment,
                                      t_of_s, d->tp->t[z_index][M_index]+segment));

                // sanity check
                HMPDFCHECK(not_monotonic(d->tp->t[z_index][M_index][segment].len,
                                         d->tp->t[z_index][M_index][segment].data,
                                         -1),
                           "theta values not monotonically decreasing in "
                           "z = %d, M = %d, segment = %d",
                           z_index, M_index, segment);

                for (long signalindex=d->tp->t[z_index][M_index][segment].start, ii=0;
                     ii < d->tp->t[z_index][M_index][segment].len;
                     signalindex += d->tp->t[z_index][M_index][segment].incr, ii++)
                {
                    au_real[signalindex] += M_PI * n
                                            * d->tp->dtsq[z_index][M_index][segment].data[ii]
                                            * d->n->Mweights[M_index] * d->n->zweights[z_index]
                                            * gsl_pow_2(d->c->comoving[z_index]) / d->c->hubble[z_index];
                    tempc_real[signalindex] += M_PI * n * b
                                               * d->tp->dtsq[z_index][M_index][segment].data[ii]
                                               * d->n->Mweights[M_index];
                }
            }
        }

        // perform FFT for clustered term
        fftw_execute(plan_c);
        // correct phases
        SAFEHMPDF(correct_phase1d(d, tempc_comp, 1));
        // write into the output array, subtracting the zero mode
        for (int ii=0; ii<(int)d->n->Nsignal/2+1; ii++)
        {
            d->tp->ac[z_index][ii] = tempc_comp[ii] - tempc_comp[0];
        }
    }
    // perform FFT for unclustered term
    fftw_execute(plan_u);
    // correct phases
    SAFEHMPDF(correct_phase1d(d, d->tp->au, 1));
    // subtract zero mode of unclustered contribution
    for (int ii=(int)d->n->Nsignal/2; ii>=0; ii--)
    {
        d->tp->au[ii] -= d->tp->au[0];
    }
    fftw_free(tempc_real);
    fftw_destroy_plan(plan_c);
    fftw_destroy_plan(plan_u);

    d->tp->created_phi_indep = 1;

    ENDFCT
}//}}}

static inline double
triang_A(double a, double b, double c)
// inverse triangle area by Heron's formula
{//{{{
    double s = 0.5 * (a + b + c);
    return pow(s * (s-a) * (s-b) * (s-c), -0.5);
}//}}}

static int
tp_segmentsum(hmpdf_obj *d, int z_index, int M_index, double phi, twopoint_workspace *ws)
{//{{{
    STARTFCT

    // read HMF and bias
    double n = d->h->hmf[z_index][M_index];
    double b = d->h->bias[z_index][M_index];

    for (int segment1 = 0;
         segment1 < d->p->segment_boundaries[z_index][M_index][0];
         segment1++)
    {
        for (int segment2 = 0;
             segment2 < d->p->segment_boundaries[z_index][M_index][0];
             segment2++)
        {
            // loop such that the theta values are always monotonically decreasing
            // so that we know when to break
            for (long signalindex1 = d->tp->t[z_index][M_index][segment1].start, ii=0;
                 ii < d->tp->t[z_index][M_index][segment1].len;
                 signalindex1 += d->tp->t[z_index][M_index][segment1].incr, ii++)
            // loop over the direction that is Nsignal long
            {
                double t1 = d->tp->t[z_index][M_index][segment1].data[ii];
                // t1 is monotonically decreasing with ii

                // check if no triangle can be formed anymore, since t1 only decreases
                if (phi >= t1 + d->p->profiles[z_index][M_index][0]) { break; }

                for (long signalindex2 = d->tp->t[z_index][M_index][segment2].start, jj=0;
                     (jj < d->tp->t[z_index][M_index][segment2].len)
                      && (signalindex2 <= signalindex1);
                     // compute only half of the matrix, because it's symmetric
                     signalindex2 += d->tp->t[z_index][M_index][segment2].incr, jj++)
                // loop over the direction that is Nsignal+2 long
                {
                    double t2 = d->tp->t[z_index][M_index][segment2].data[jj];
                    // t2 is monotonically decreasing with jj

                    // check if we can form a triangle
                    if (t1 >= t2 + phi) { break; }    // t2 only decreases
                    if (t2 >= phi + t1) { continue; } // no triangle possible here, but perhaps later
                    if (phi >= t1 + t2) { break; }    // t1,t2 only decrease

                    double Delta = triang_A(phi, t1, t2);
                    double temp = 0.25 * n * Delta
                                  * d->tp->dtsq[z_index][M_index][segment1].data[ii]
                                  * d->tp->dtsq[z_index][M_index][segment2].data[jj]
                                  * d->n->Mweights[M_index];

                    // add to clustered term
                    ws->tempc_real[signalindex1*(d->n->Nsignal+2)+signalindex2]
                        += temp * b;

                    // add to unclustered term
                    ws->pdf_real[signalindex1*(d->n->Nsignal+2)+signalindex2]
                        += temp * gsl_pow_2(d->c->comoving[z_index])
                           / d->c->hubble[z_index] * d->n->zweights[z_index];
                }
            }
        }
    }

    ENDFCT
}//}}}

static int
tp_Mint(hmpdf_obj *d, int z_index, double phi, twopoint_workspace *ws)
// adds to pdf_real, with the required zweight * Mweight, including the unclustered 1pt PDF contributions
// creates new tempc_real (nulls first) --> tempc created with fftw_malloc
// pdf_real, tempc_real are not symmetrized!
{//{{{
    STARTFCT

    // zero tempc
    zero_real(d->n->Nsignal*(d->n->Nsignal+2), ws->tempc_real);

    for (int M_index=0; M_index<d->n->NM; M_index++)
    {
        SAFEHMPDF(tp_segmentsum(d, z_index, M_index, phi, ws));
    }

    ENDFCT
}//}}}

static inline double complex
redundant(long N, double complex *a, long ii)
// extends the vector a[N/2+1] to N elements through conjugation
// N is assumed even
{//{{{
    if (ii <= N/2)
    {
        return a[ii];
    }
    else
    {
        return conj(a[N-ii]);
    }
}//}}}

static inline int
clustered_term(hmpdf_obj *d, int z_index, double corr_phi_2, double corr_phi,
               long i1/*long direction*/, long i2/*short direction*/,
               double complex *b12, double complex *out)
// computes 1/2 * (alpha1^2 + alpha2^2) * zeta(0)
//          + alpha1 * alpha2 * zeta(phi)
//          + 1/2 * beta12^2 * zeta(0)
//          + beta12 * (alpha1 + alpha2) * zeta(phi/2)
// takes care of the zero modes in b12 (a1, a2 are already zeroed)
{//{{{
    STARTFCT

    double complex a1 = redundant(d->n->Nsignal, d->tp->ac[z_index], i1);
    double complex a2 = d->tp->ac[z_index][i2];
    double complex b = b12[i1*(d->n->Nsignal/2+1)+i2]
                       - b12[i1*(d->n->Nsignal/2+1)]
                       - b12[i2] + b12[0];

    *out = 0.5 * (a1*a1 + a2*a2 + b*b) * d->c->Dsq[z_index] * d->pwr->autocorr
           + a1*a2 * corr_phi
           + b*(a1 + a2) * corr_phi_2;

    ENDFCT
}//}}}

static int
symmetrize(hmpdf_obj *d, double *A)
// fills the lower triangular part of A[Nsignal+2,Nsignal]
//     with the upper triangular part
{//{{{
    STARTFCT

    for (long ii=0; ii<d->n->Nsignal; ii++)
    {
        for (long jj=0; jj<ii; jj++)
        {
            A[jj*(d->n->Nsignal+2)+ii]
                = A[ii*(d->n->Nsignal+2)+jj];
        }
    }

    ENDFCT
}//}}}

static int
tp_zint(hmpdf_obj *d, double phi, twopoint_workspace *ws)
// z-integral of the unclustered terms, without FFT 
// z-integral of the clustered terms, including FFT (of course)
{//{{{
    STARTFCT

    // zero the integrals
    zero_real(d->n->Nsignal*(d->n->Nsignal+2),   ws->pdf_real);
    zero_comp(d->n->Nsignal*(d->n->Nsignal/2+1), ws->bc);

    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        // perform mass integration
        SAFEHMPDF(tp_Mint(d, z_index, phi, ws));

        // symmetrize the clustered beta matrix
        SAFEHMPDF(symmetrize(d, ws->tempc_real));

        // perform the FFT on the clustered part tempc_real -> tempc_comp
        fftw_execute(ws->pc_r2c);
        // correct phases
        SAFEHMPDF(correct_phase2d(d, ws->tempc_comp, 1));

        // compute the correlation function interpolator
        double corr_phi_2, corr_phi;
        SAFEHMPDF(corr(d, z_index, 0.5*phi, &corr_phi_2));
        SAFEHMPDF(corr(d, z_index, phi, &corr_phi));

        // add to the clustered output
        for (long ii=0; ii<d->n->Nsignal; ii++)
        // loop over the long direction
        {
            for (long jj=0; jj<d->n->Nsignal/2+1; jj++)
            // loop over the short direction
            {
                double complex clterm;
                SAFEHMPDF(clustered_term(d, z_index, corr_phi_2, corr_phi,
                                         ii, jj, ws->tempc_comp, &clterm));
                double complex temp = clterm
                                      * gsl_pow_4(d->c->comoving[z_index])
                                      / d->c->hubble[z_index];
                ws->bc[ii*(d->n->Nsignal/2+1)+jj] += temp * d->n->zweights[z_index];
            }
        }
    }

    ENDFCT
}//}}}

int
create_tp(hmpdf_obj *d, double phi, twopoint_workspace *ws)
// computes one 2pt PDF
{//{{{
    STARTFCT

    // perform the redshift integration
    SAFEHMPDF(tp_zint(d, phi, ws));

    // symmetrize the unclustered part
    SAFEHMPDF(symmetrize(d, ws->pdf_real));
    
    // perform the FFT on the unclustered part pdf_real -> pdf_comp
    fftw_execute(ws->pu_r2c);
    // correct phases
    SAFEHMPDF(correct_phase2d(d, ws->pdf_comp, 1));

    // add the clustering contribution,
    // subtract the zero modes in the unclustered part,
    // take exponential,
    // and normalize properly
    // loop backwards so we don't have to store the zero modes elsewhere
    for (long ii=d->n->Nsignal-1; ii>=0; ii--)
    {
        for (long jj=d->n->Nsignal/2; jj>=0; jj--)
        {
            ws->pdf_comp[ii*(d->n->Nsignal/2+1)+jj]
                = cexp(ws->pdf_comp[ii*(d->n->Nsignal/2+1)+jj]
                       - ws->pdf_comp[ii*(d->n->Nsignal/2+1)]
                       - ws->pdf_comp[jj] + ws->pdf_comp[0]
                       + d->tp->au[jj] + redundant(d->n->Nsignal, d->tp->au, ii)
                       + ws->bc[ii*(d->n->Nsignal/2+1)+jj])
                  / gsl_pow_2((double)(d->n->Nsignal));
        }
    }

    // correct phases
    SAFEHMPDF(correct_phase2d(d, ws->pdf_comp, -1));
    // perform backward FFT pdf_comp -> pdf_real
    fftw_execute(ws->ppdf_c2r);

    ENDFCT
}//}}}

static int
create_noisy_tp(hmpdf_obj *d, double phi)
{//{{{
    STARTFCT

    if (d->tp->pdf_noisy == NULL)
    {
        SAFEALLOC(d->tp->pdf_noisy,
                  malloc(d->n->Nsignal_noisy
                         * d->n->Nsignal_noisy
                         * sizeof(double)));
    }

    SAFEHMPDF(noise_matr(d, d->tp->pdf, d->tp->pdf_noisy,
                         0/*not buffered*/, phi));

    ENDFCT
}//}}}

// convenience macro which cleans up if memory allocation fails,
// since we do not consider this a critical error.
#define NEWTPWS_SAFEALLOC(var, expr)                                    \
    do {                                                                \
        var = expr;                                                     \
        if (UNLIKELY(!(var)))                                           \
        {                                                               \
            if (ws->pdf_real != NULL) { fftw_free(ws->pdf_real); }      \
            if (ws->bc != NULL) { free(ws->bc); }                       \
            if (ws->tempc_real != NULL) { fftw_free(ws->tempc_real); }  \
            free(*out);                                                 \
            return 1;                                                   \
        }                                                               \
    } while (0)

int
new_tp_ws(long N, twopoint_workspace **out)
{//{{{
    STARTFCT

    SAFEALLOC(*out, malloc(sizeof(twopoint_workspace)));

    twopoint_workspace *ws = *out; // for convenience

    // initialize to NULL so we can free reliably in case an alloc fails
    ws->pdf_real = NULL;
    ws->bc = NULL;
    ws->tempc_real = NULL;

    // do the allocs first so we don't have to worry
    // about whether fftw_plans have already been computed and need
    // to be destroyed if an alloc fails
    NEWTPWS_SAFEALLOC(ws->pdf_real, fftw_malloc(N * (N+2) * sizeof(double)));
    ws->pdf_comp = (double complex *)(ws->pdf_real);

    NEWTPWS_SAFEALLOC(ws->bc, malloc(N * (N/2+1) * sizeof(double complex)));

    NEWTPWS_SAFEALLOC(ws->tempc_real, fftw_malloc(N * (N+2) * sizeof(double)));
    ws->tempc_comp = (double complex *)(ws->tempc_real);

    ws->pu_r2c = fftw_plan_dft_r2c_2d(N, N, ws->pdf_real,
                                      ws->pdf_comp, PU_R2C_MODE);
    ws->ppdf_c2r = fftw_plan_dft_c2r_2d(N, N, ws->pdf_comp,
                                        ws->pdf_real, PPDF_C2R_MODE);
    ws->pc_r2c = fftw_plan_dft_r2c_2d(N, N, ws->tempc_real,
                                      ws->tempc_comp, PC_R2C_MODE);

    ENDFCT
}//}}}

#undef NEWTPWS_SAFEALLOC

static int
prepare_tp(hmpdf_obj *d, double phi)
{//{{{
    STARTFCT

    HMPDFPRINT(1, "prepare_tp\n");

    // run necessary code from other modules
    SAFEHMPDF(create_corr(d));
    if (d->f->Nfilters > 0)
    {
        SAFEHMPDF(create_conj_profiles(d));
        SAFEHMPDF(create_filtered_profiles(d));
    }
    SAFEHMPDF(create_segments(d));
    SAFEHMPDF(create_phi_indep(d));
    SAFEHMPDF(create_op(d));

    if (d->tp->ws == NULL)
    {
        SAFEHMPDF(new_tp_ws(d->n->Nsignal, &(d->tp->ws)));
        HMPDFCHECK(d->tp->ws==NULL, "OOM.");
    }

    // convert from arcmin to radians
    phi *= RADPERARCMIN;
    SAFEHMPDF(create_tp(d, phi, d->tp->ws));
    
    // copy PDF into contiguous array (tp->pdf_real has padding from the FFTs)
    if (d->tp->pdf == NULL)
    {
        SAFEALLOC(d->tp->pdf, malloc(d->n->Nsignal*d->n->Nsignal*sizeof(double)));
    }
    for (long ii=0; ii<d->n->Nsignal; ii++)
    {
        memcpy(d->tp->pdf+ii*d->n->Nsignal,
               d->tp->ws->pdf_real+ii*(d->n->Nsignal+2),
               d->n->Nsignal * sizeof(double));
    }

    if (d->ns->have_noise)
    {
        SAFEHMPDF(create_noise_matr_conv(d, 1/*need only one buffer*/));
        SAFEHMPDF(create_noisy_tp(d, phi));
    }

    ENDFCT
}//}}}

int
hmpdf_get_tp(hmpdf_obj *d, double phi, int Nbins, double binedges[Nbins+1], double tp[Nbins*Nbins], int noisy)
{//{{{
    STARTFCT

    CHECKINIT;

    SAFEHMPDF(pdf_check_user_input(d, Nbins, binedges, noisy));

    // perform computation if necessary
    if (fabs(1.0 - d->tp->last_phi/phi) > TP_PHI_EQ_TOL)
    {
        SAFEHMPDF(prepare_tp(d, phi));
    }
    d->tp->last_phi = phi;

    double _binedges[Nbins+1];
    SAFEHMPDF(pdf_adjust_binedges(d, Nbins, binedges, _binedges, d->op->signalmeanc));

    HMPDFPRINT(3, "\t\tbinning the twopoint pdf\n");
    SAFEHMPDF(bin_2d((noisy) ? d->n->Nsignal_noisy : d->n->Nsignal,
                     (noisy) ? d->n->signalgrid_noisy : d->n->signalgrid,
                     (noisy) ? d->tp->pdf_noisy : d->tp->pdf,
                     TPINTEGR_N, Nbins, _binedges, tp, TPINTERP_TYPE));

    ENDFCT
}//}}}

