#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include <gsl/gsl_math.h>

#include "utils.h"
#include "configs.h"
#include "data.h"
#include "power.h"
#include "profiles.h"
#include "onepoint.h"

#include "hmpdf.h"

int
null_twopoint(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

    d->tp->created_phi_indep = 0;
    d->tp->dtsq = NULL;
    d->tp->t = NULL;
    d->tp->ac = NULL;
    d->tp->au = NULL;
    d->tp->ws = NULL;
    d->tp->last_phi = -1.0;
    d->tp->pdf = NULL;
    d->tp->pdf_noisy = NULL;

    CHECKERR
    return hmpdf_status;
}//}}}

int
reset_twopoint(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

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

    CHECKERR
    return hmpdf_status;
}//}}}

int
create_phi_indep(hmpdf_obj *d)
// computes tp->dtsq, tp->t, tp->ac
{//{{{
    int hmpdf_status = 0;

    if (d->tp->created_phi_indep) { return hmpdf_status; }

    fprintf(stdout, "\tcreate_phi_indep\n");
    fflush(stdout);
    
    SAFEALLOC(, d->tp->dtsq, malloc(d->n->Nz * sizeof(double **)))
    SAFEALLOC(, d->tp->t,    malloc(d->n->Nz * sizeof(double **)))
    SAFEALLOC(, d->tp->ac, malloc(d->n->Nz * sizeof(complex *)))
    SAFEALLOC(, d->tp->au, fftw_malloc((d->n->Nsignal/2+1) * sizeof(complex)))
    double *au_real = (double *)(d->tp->au);

    SAFEALLOC(double *, tempc_real, fftw_malloc((d->n->Nsignal+2)*sizeof(double)))
    complex *tempc_comp = (complex *)tempc_real;
    fftw_plan plan_c = fftw_plan_dft_r2c_1d(d->n->Nsignal, tempc_real,
                                            tempc_comp, FFTW_MEASURE);
    fftw_plan plan_u = fftw_plan_dft_r2c_1d(d->n->Nsignal, au_real,
                                            d->tp->au, FFTW_MEASURE);

    for (int ii=0; ii<d->n->Nsignal/2+1; ii++)
    {
        d->tp->au[ii] = 0.0;
    }

    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        SAFEALLOC(, d->tp->dtsq[z_index], malloc(d->n->NM * sizeof(double *)))
        SAFEALLOC(, d->tp->t[z_index],    malloc(d->n->NM * sizeof(double *)))
        SAFEALLOC(, d->tp->ac[z_index],   malloc((d->n->Nsignal/2+1) * sizeof(complex)))
        // null the FFT array
        for (int ii=0; ii<d->n->Nsignal+2; ii++)
        {
            tempc_real[ii] = 0.0;
        }
        // integrate clustered contribution over mass;
        // and fill the theta(signal) dtheta^2/dsignal grids
        for (int M_index=0; M_index<d->n->NM; M_index++)
        {
            SAFEALLOC(, d->tp->dtsq[z_index][M_index],
                      malloc(d->n->Nsignal * sizeof(double)))
            SAFEALLOC(, d->tp->t[z_index][M_index],
                      malloc(d->n->Nsignal * sizeof(double)))

            SAFEHMPDF(dtsq_of_s(d, z_index, M_index, d->tp->dtsq[z_index][M_index]))
            SAFEHMPDF(t_of_s(d, z_index, M_index, d->tp->t[z_index][M_index]))
            
            double n = d->h->hmf[z_index][M_index];
            double b = d->h->bias[z_index][M_index];
            for (int ii=0; ii<d->n->Nsignal; ii++)
            {
                au_real[ii] -= M_PI * n * d->tp->dtsq[z_index][M_index][ii]
                               * d->n->Mweights[M_index] * d->n->zweights[z_index]
                               * gsl_pow_2(d->c->comoving[z_index]) / d->c->hubble[z_index];
                tempc_real[ii] -= M_PI * n * b
                                  * d->tp->dtsq[z_index][M_index][ii]
                                  * d->n->Mweights[M_index];
            }
        }

        // perform FFT
        fftw_execute(plan_c);
        // write into the output array, subtracting the zero mode
        for (int ii=0; ii<d->n->Nsignal/2+1; ii++)
        {
            d->tp->ac[z_index][ii] = tempc_comp[ii] - tempc_comp[0];
        }
    }
    fftw_execute(plan_u);
    // subtract zero mode of unclustered contribution
    for (int ii=d->n->Nsignal/2; ii>=0; ii--)
    {
        d->tp->au[ii] -= d->tp->au[0];
    }
    fftw_free(tempc_real);
    fftw_destroy_plan(plan_c);
    fftw_destroy_plan(plan_u);

    d->tp->created_phi_indep = 1;

    CHECKERR
    return hmpdf_status;
}//}}}

static double
triang_A(double a, double b, double c)
// inverse triangle area by Heron's formula
{//{{{
    double s = 0.5 * (a + b + c);
    return pow(s * (s-a) * (s-b) * (s-c), -0.5);
}//}}}

static int
tp_Mint(hmpdf_obj *d, int z_index, double phi, twopoint_workspace *ws)
// adds to pdf_real, with the required zweight * Mweight, including the unclustered 1pt PDF contributions
// creates new tempc_real (nulls first) --> tempc created with fftw_malloc
// pdf_real, tempc_real are not symmetrized!
{//{{{
    int hmpdf_status = 0;

    // null tempc
    for (int ii=0; ii<d->n->Nsignal*(d->n->Nsignal+2); ii++)
    {
        ws->tempc_real[ii] = 0.0;
    }

    for (int M_index=0; M_index<d->n->NM; M_index++)
    {
        double n = d->h->hmf[z_index][M_index];
        double b = d->h->bias[z_index][M_index];
        for (int ii=0; ii<d->n->Nsignal; ii++)
        // loop over the direction that is Nsignal long
        {
            double t1 = d->tp->t[z_index][M_index][ii]; // t1 is monotonically decreasing with ii
            // check if no triangle can be formed anymore, since t1 only increases
            if (phi >= t1 + d->p->profiles[z_index][M_index][0]) { break; }
            for (int jj=0; jj<=ii; jj++)
            // loop over the direction that is Nsignal+2 long
            {
                double t2 = d->tp->t[z_index][M_index][jj]; // t2 is monotonically decreasing with jj
                // check if we can form a triangle
                // we don't need to check for t1 > t2 + phi, since t2 >= t1 by construction
                if (t2 >= phi + t1) { continue; }
                if (phi >= t1 + t2) { break; }
                //double min_diff = GSL_MIN(phi+t1-t2, t1+t2-phi);
                double Delta = triang_A(phi, t1, t2);
                double temp = 0.25 * n * Delta
                              * d->tp->dtsq[z_index][M_index][ii]
                              * d->tp->dtsq[z_index][M_index][jj]
                              * d->n->Mweights[M_index];
                
                ws->tempc_real[ii*(d->n->Nsignal+2)+jj] += temp * b;
                ws->pdf_real[ii*(d->n->Nsignal+2)+jj]
                    += temp * gsl_pow_2(d->c->comoving[z_index])
                       / d->c->hubble[z_index] * d->n->zweights[z_index];
            }
        }
        // add the 1pt PDF unclustered contributions
        /*
        for (int ii=0; ii<d->n->Nsignal; ii++)
        {
            ws->pdf_real[ii*(d->n->Nsignal+2)]
                -=  d->tp->dtsq[z_index][M_index][ii]
                    * M_PI * n
                    * gsl_pow_2(d->c->comoving[z_index]) / d->c->hubble[z_index]
                    * d->n->Mweights[M_index] * d->n->zweights[z_index]
                    * ((ii==0) ? 2.0 : 1.0);
        }
        */
    }

    CHECKERR
    return hmpdf_status;
}//}}}

static complex
redundant(int N, complex *a, int ii)
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

static int
clustered_term(hmpdf_obj *d, int z_index, double phi,
                       int i1/*long direction*/, int i2/*short direction*/,
                       complex *b12, complex *out)
// computes 1/2 * (alpha1^2 + alpha2^2) * zeta(0)
//          + alpha1 * alpha2 * zeta(phi)
//          + 1/2 * beta12^2 * zeta(0)
//          + beta12 * (alpha1 + alpha2) * zeta(phi/2)
// takes care of the zero modes in b12 (a1, a2 are already zeroed)
{//{{{
    int hmpdf_status = 0;

    *out = 0.0; // to avoid maybe-uninitialized

    complex a1 = redundant(d->n->Nsignal, d->tp->ac[z_index], i1);
    complex a2 = d->tp->ac[z_index][i2];
    complex b = b12[i1*(d->n->Nsignal/2+1)+i2]
                - b12[i1*(d->n->Nsignal/2+1)]
                - b12[i2] + b12[0];
    double corr1, corr2;
    SAFEHMPDF(corr(d, z_index, phi, &corr1))
    SAFEHMPDF(corr(d, z_index, 0.5*phi, &corr2))
    *out = 0.5 * (a1*a1 + a2*a2 + b*b) * d->c->Dsq[z_index] * d->pwr->autocorr
           + a1*a2 * corr1
           + b*(a1 + a2) * corr2;

    CHECKERR
    return hmpdf_status;
}//}}}

static int
tp_zint(hmpdf_obj *d, double phi, twopoint_workspace *ws)
// z-integral of the unclustered terms, without FFT 
// z-integral of the clustered terms, including FFT (of course)
{//{{{
    int hmpdf_status = 0;

    // zero the integrals
    for (int ii=0; ii<d->n->Nsignal*(d->n->Nsignal+2); ii++)
    {
        ws->pdf_real[ii] = 0.0;
    }
    for (int ii=0; ii<d->n->Nsignal*(d->n->Nsignal/2+1); ii++)
    {
        ws->bc[ii] = 0.0;
    }

    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        // perform mass integration
        SAFEHMPDF(tp_Mint(d, z_index, phi, ws))

        // symmetrize the clustered beta matrix
        for (int ii=0; ii<d->n->Nsignal; ii++)
        {
            for (int jj=0; jj<ii; jj++)
            {
                ws->tempc_real[jj*(d->n->Nsignal+2)+ii]
                    = ws->tempc_real[ii*(d->n->Nsignal+2)+jj];
            }
        }
        // perform the FFT on the clustered part tempc_real -> tempc_comp
        fftw_execute(ws->pc_r2c);

        // add to the clustered output
        for (int ii=0; ii<d->n->Nsignal; ii++)
        // loop over the long direction
        {
            for (int jj=0; jj<d->n->Nsignal/2+1; jj++)
            // loop over the short direction
            {
                complex clterm;
                SAFEHMPDF(clustered_term(d, z_index, phi, ii, jj, ws->tempc_comp, &clterm))
                complex temp = clterm
                               * gsl_pow_4(d->c->comoving[z_index])
                               / d->c->hubble[z_index];
                ws->bc[ii*(d->n->Nsignal/2+1)+jj] += temp * d->n->zweights[z_index];
            }
        }
    }

    CHECKERR
    return hmpdf_status;
}//}}}

int
create_tp(hmpdf_obj *d, double phi, twopoint_workspace *ws)
// computes one 2pt PDF
{//{{{
    int hmpdf_status = 0;

    // perform the redshift integration
    SAFEHMPDF(tp_zint(d, phi, ws))

    // symmetrize the unclustered part
    for (int ii=0; ii<d->n->Nsignal; ii++)
    {
        for (int jj=0; jj<ii; jj++)
        {
            ws->pdf_real[jj*(d->n->Nsignal+2)+ii]
                = ws->pdf_real[ii*(d->n->Nsignal+2)+jj];
        }
    }
    
    // perform the FFT on the unclustered part pdf_real -> pdf_comp
    fftw_execute(ws->pu_r2c);

    // add the clustering contribution,
    // subtract the zero modes in the unclustered part,
    // take exponential,
    // and normalize properly
    // loop backwards so we don't have to store the zero modes elsewhere
    for (int ii=d->n->Nsignal-1; ii>=0; ii--)
    {
        for (int jj=d->n->Nsignal/2; jj>=0; jj--)
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

    // perform backward FFT pdf_comp -> pdf_real
    fftw_execute(ws->ppdf_c2r);

    CHECKERR
    return hmpdf_status;
}//}}}

int
create_noisy_tp(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

    if (d->tp->pdf_noisy == NULL)
    {
        SAFEALLOC(, d->tp->pdf_noisy,
                  malloc(d->n->Nsignal_noisy
                         * d->n->Nsignal_noisy
                         * sizeof(double)))
    }

    SAFEHMPDF(noise_matr(d, d->tp->pdf, d->tp->pdf_noisy))

    CHECKERR
    return hmpdf_status;
}//}}}

int
new_tp_ws(int N, twopoint_workspace **out)
{//{{{
    int hmpdf_status = 0;
    SAFEALLOC(, *out, malloc(sizeof(twopoint_workspace)))
    twopoint_workspace *ws = *out; // for convenience

    SAFEALLOC_NORETURN(, ws->pdf_real, fftw_malloc(N * (N+2) * sizeof(double)))
    if (hmpdf_status)
    {//{{{
        free(*out);
        return hmpdf_status;
    }//}}}
    ws->pdf_comp = (complex *)(ws->pdf_real);
    ws->pu_r2c = fftw_plan_dft_r2c_2d(N, N, ws->pdf_real,
                                      ws->pdf_comp, PU_R2C_MODE);
    ws->ppdf_c2r = fftw_plan_dft_c2r_2d(N, N, ws->pdf_comp,
                                        ws->pdf_real, PPDF_C2R_MODE);
    
    SAFEALLOC_NORETURN(, ws->bc, malloc(N * (N/2+1) * sizeof(complex)))
    if (hmpdf_status)
    {//{{{
        fftw_destroy_plan(ws->ppdf_c2r);
        fftw_destroy_plan(ws->pu_r2c);
        fftw_free(ws->pdf_real);
        free(*out);
        return hmpdf_status;
    }//}}}

    SAFEALLOC_NORETURN(, ws->tempc_real, fftw_malloc(N * (N+2) * sizeof(double)))
    if (hmpdf_status)
    {//{{{
        free(ws->bc);
        fftw_destroy_plan(ws->ppdf_c2r);
        fftw_destroy_plan(ws->pu_r2c);
        fftw_free(ws->pdf_real);
        free(*out);
        return hmpdf_status;
    }//}}}
    ws->tempc_comp = (complex *)(ws->tempc_real);
    ws->pc_r2c = fftw_plan_dft_r2c_2d(N, N, ws->tempc_real, ws->tempc_comp, PC_R2C_MODE);

    CHECKERR
    return hmpdf_status;
}//}}}

static int
prepare_tp(hmpdf_obj *d, double phi)
{//{{{
    int hmpdf_status = 0;

    fprintf(stdout, "In twopoint.h -> prepare_tp :\n");
    fflush(stdout);

    if (!(d->n->monotonize))
    {
        fprintf(stderr, "Error: Twopoint only possible if monotonize=1.\n");
        fflush(stderr);
        ERRLOC
        hmpdf_status |= 1;
        return hmpdf_status;
    }

    // run necessary code from other modules
    SAFEHMPDF(create_corr(d))
    if (d->f->Nfilters > 0)
    {
        SAFEHMPDF(create_conj_profiles(d))
        SAFEHMPDF(create_filtered_profiles(d))
    }
    SAFEHMPDF(create_breakpoints_or_monotonize(d))
    SAFEHMPDF(create_phi_indep(d))
    SAFEHMPDF(create_op(d))

    if (d->tp->ws == NULL)
    {
        SAFEHMPDF(new_tp_ws(d->n->Nsignal, &(d->tp->ws)))
        if (d->tp->ws == NULL)
        {
            fprintf(stderr, "Error: Out of memory in twopoint.h -> prepare_tp.\n");
            fflush(stderr);
            ERRLOC
            hmpdf_status |= 1;
            return hmpdf_status;
        }
    }

    SAFEHMPDF(create_tp(d, phi, d->tp->ws))
    
    // copy PDF into contiguous array (tp->pdf_real has padding from the FFTs)
    if (d->tp->pdf == NULL)
    {
        SAFEALLOC(, d->tp->pdf, malloc(d->n->Nsignal*d->n->Nsignal*sizeof(double)))
    }
    for (int ii=0; ii<d->n->Nsignal; ii++)
    {
        memcpy(d->tp->pdf+ii*d->n->Nsignal,
               d->tp->ws->pdf_real+ii*(d->n->Nsignal+2),
               d->n->Nsignal * sizeof(double));
    }

    if (d->ns->noise > 0.0)
    {
        SAFEHMPDF(create_noisy_tp(d))
    }

    CHECKERR
    return hmpdf_status;
}//}}}

int
hmpdf_get_tp(hmpdf_obj *d, double phi, int Nbins, double *binedges, double *out, int noisy)
{//{{{
    int hmpdf_status = 0;

    if (noisy && d->ns->noise<0.0)
    {
        fprintf(stderr, "Error: noisy twopoint pdf requested but no/invalid noise level passed.\n");
        fflush(stderr);
        ERRLOC
        hmpdf_status |= 1;
        return hmpdf_status;
    }

    if (not_monotonic(Nbins+1, binedges, NULL))
    {
        fprintf(stderr, "Error: binedges not monotonically increasing.\n");
        fflush(stderr);
        ERRLOC
        hmpdf_status |= 1;
        return hmpdf_status;
    }

    // perform computation if necessary
    if (fabs(1.0 - d->tp->last_phi/phi) > TP_PHI_EQ_TOL)
    {
        SAFEHMPDF(prepare_tp(d, phi))
    }
    d->tp->last_phi = phi;

    double _binedges[Nbins+1];
    memcpy(_binedges, binedges, (Nbins+1) * sizeof(double));
    // if kappa, adjust the bins
    if (d->p->stype == hmpdf_kappa)
    {
        for (int ii=0; ii<= Nbins; ii++)
        {
            _binedges[ii] += d->op->signalmeanc;
        }
    }

    fprintf(stdout, "\t\tbinning the twopoint pdf\n");
    fflush(stdout);
    SAFEHMPDF(bin_2d((noisy) ? d->n->Nsignal_noisy : d->n->Nsignal,
                     (noisy) ? d->n->signalgrid_noisy : d->n->signalgrid,
                     (noisy) ? d->tp->pdf_noisy : d->tp->pdf,
                     TPINTEGR_N, Nbins, _binedges, out, TPINTERP_TYPE))

    CHECKERR
    return hmpdf_status;
}//}}}

