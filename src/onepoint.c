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
#include "profiles.h"
#include "numerics.h"
#include "onepoint.h"

static
void op_Mint_invertible(all_data *d, int z_index, double *au, double *ac)
// performs the mass integrals, for the invertible profiles
// integrals are _added_ to au, ac
{//{{{
    double *temp = (double *)malloc(d->n->gr->Nsignal * sizeof(double));

    for (int M_index=breakpoints(d, z_index); M_index<d->n->gr->NM; M_index++)
    {
        dtsq_of_s(d, z_index, M_index, temp);
        double n = d->h->hmf[z_index*d->n->gr->NM + M_index];
        double b = d->h->bias[z_index*d->n->gr->NM + M_index];
        for (int ii=0; ii<d->n->gr->Nsignal; ii++)
        {
            au[ii] -= temp[ii] * M_PI * n
                      * d->n->gr->Mweights[M_index];
            ac[ii] -= temp[ii] * M_PI * n * b
                      * d->n->gr->Mweights[M_index];
        }
    }

    free(temp);
}//}}}

static
complex not_inv_integral(all_data *d, int z_index, int M_index, int lambda_index)
// TODO do this w/ Gauss fixed point (linear weight function)
{//{{{
    complex *integr = (complex *)malloc(d->p->prtilde_Ntheta * sizeof(complex));
    double *temp = (double *)malloc(d->p->prtilde_Ntheta * sizeof(double));
    // evaluate the interpolated signal profile
    s_of_t(d, z_index, M_index, d->p->prtilde_Ntheta, d->p->prtilde_thetagrid, temp);
    // fill the integrand
    for (int ii=0; ii<d->p->prtilde_Ntheta; ii++)
    {
        integr[ii] = d->p->prtilde_thetagrid[ii]
                     * cexp(- _Complex_I * temp[ii] * d->n->gr->lambdagrid[lambda_index]);
    }
    // TODO TESTING
    /*
    if ((z_index+1)%10==0 && (M_index+1)%10==0 && (lambda_index+1)%200==0
        && z_index>25 && M_index>35)
    {
        printf("z=%d/%d\tM=%d/%d\tlambda=%d/%d\n", z_index+1, d->n->gr->Nz,
                                                   M_index+1, d->n->gr->NM,
                                                   lambda_index+1, d->n->gr->Nsignal/2+1);
        gnuplot *gp = plot_comp(NULL, d->p->prtilde_Ntheta, d->p->prtilde_thetagrid, integr, 0);
        plot_comp(gp, d->p->prtilde_Ntheta, d->p->prtilde_thetagrid, integr, 1);
        show(gp);
    }
    */
    // perform the integration
    complex res = integr_comp(d->p->prtilde_Ntheta, 1.0/(double)(d->p->prtilde_Ntheta-1), 1, integr);
    free(integr);
    free(temp);
    return res * 2.0 * M_PI * gsl_pow_2(d->p->profiles[z_index][M_index][0]);
}//}}}

static
void lambda_loop(all_data *d, int z_index, int M_index, complex *out)
{//{{{
    for (int ii=0; ii<d->n->gr->Nsignal/2+1; ii++)
    {
        out[ii] = not_inv_integral(d, z_index, M_index, ii);
    }
}//}}}

static
void op_Mint_notinvertible(all_data *d, int z_index, complex *au, complex *ac)
// performs the mass integration over the non-invertible profiles
// integrals are _added_ to au, ac
{//{{{
    complex *temp = (complex *)malloc((d->n->gr->Nsignal/2+1) * sizeof(complex));
    for (int M_index=0; M_index<breakpoints(d, z_index); M_index++)
    {
        lambda_loop(d, z_index, M_index, temp);
        double n = d->h->hmf[z_index*d->n->gr->NM + M_index];
        double b = d->h->bias[z_index*d->n->gr->NM + M_index];
        for (int ii=0; ii<d->n->gr->Nsignal/2+1; ii++)
        {
            au[ii] += temp[ii] * n
                      * d->n->gr->Mweights[M_index];
            ac[ii] += temp[ii] * n * b
                      * d->n->gr->Mweights[M_index];
        }
    }

    free(temp);
}//}}}

static
void op_zint(all_data *d, complex *pu_comp, complex *pc_comp) // p is the exponent in P(lambda)
{//{{{
    double *ac_real = (double *)fftw_malloc((d->n->gr->Nsignal+2) * sizeof(double));
    complex *ac_comp = (complex *)ac_real;
    double *au_real = (double *)fftw_malloc((d->n->gr->Nsignal+2) * sizeof(double));
    complex *au_comp = (complex *)au_real;
    fftw_plan plan_u = fftw_plan_dft_r2c_1d(d->n->gr->Nsignal, au_real, au_comp, FFTW_MEASURE);
    fftw_plan plan_c = fftw_plan_dft_r2c_1d(d->n->gr->Nsignal, ac_real, ac_comp, FFTW_MEASURE);

    // zero the output arrays
    for (int ii=0; ii<d->n->gr->Nsignal/2+1; ii++)
    {
        pu_comp[ii] = pc_comp[ii] = 0.0;
    }

    // fill the z-integrand
    for (int z_index=0; z_index<d->n->gr->Nz; z_index++)
    {
        // zero the arrays
        for (int ii=0; ii<d->n->gr->Nsignal/2+1; ii++)
        {
            au_comp[ii] = ac_comp[ii] = 0.0;
        }

        if (breakpoints(d, z_index) < d->n->gr->NM-1)
        // there are samples for which we can do the FFT integral
        {
            op_Mint_invertible(d, z_index, au_real, ac_real);
            fftw_execute(plan_u);
            fftw_execute(plan_c);
        }
        if (breakpoints(d, z_index) > 0)
        // there are samples for which we need to do the slow integral
        {
            op_Mint_notinvertible(d, z_index, au_comp, ac_comp);
        }

        for (int ii=0; ii<d->n->gr->Nsignal/2+1; ii++)
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

            pu_comp[ii] += tempu * d->n->gr->zweights[z_index];
            pc_comp[ii] += tempc * d->n->gr->zweights[z_index];
        }
    }
    fftw_free(au_real);
    fftw_free(ac_real);
    fftw_destroy_plan(plan_u);
    fftw_destroy_plan(plan_c);
}//}}}

static
void get_mean_signal(all_data *d)
{//{{{
    d->op->signalmeanu = 0.0;
    d->op->signalmeanc = 0.0;
    double normu = 0.0;
    double normc = 0.0;
    for (int ii=0; ii<d->n->gr->Nsignal; ii++)
    {
        d->op->signalmeanu += d->n->gr->signalgrid[ii] * d->op->PDFu[ii];
        d->op->signalmeanc += d->n->gr->signalgrid[ii] * d->op->PDFc[ii];
        normu += d->op->PDFu[ii];
        normc += d->op->PDFc[ii];
    }
    d->op->signalmeanu /= normu;
    d->op->signalmeanc /= normc;
}//}}}

void create_op(all_data *d)
{//{{{
    if (d->op->created_op) { return; }
    printf("\tcreate_op\n");
    d->op->PDFu = (double *)fftw_malloc((d->n->gr->Nsignal + 2) * sizeof(double));
    d->op->PDFc = (double *)fftw_malloc((d->n->gr->Nsignal + 2) * sizeof(double));
    complex *PDFu_comp = (complex *)d->op->PDFu;
    complex *PDFc_comp = (complex *)d->op->PDFc;

    fftw_plan plan_u = fftw_plan_dft_c2r_1d(d->n->gr->Nsignal, PDFu_comp, d->op->PDFu, FFTW_ESTIMATE);
    fftw_plan plan_c = fftw_plan_dft_c2r_1d(d->n->gr->Nsignal, PDFc_comp, d->op->PDFc, FFTW_ESTIMATE);

    // perform redshift integration
    op_zint(d, PDFu_comp, PDFc_comp);

    // take exponential and normalize
    for (int ii=0; ii<d->n->gr->Nsignal/2+1; ii++)
    {
        PDFu_comp[ii] = cexp(PDFu_comp[ii])/(double)(d->n->gr->Nsignal);
        PDFc_comp[ii] = cexp(PDFc_comp[ii])/(double)(d->n->gr->Nsignal);
    }

    // transform back to real space
    fftw_execute(plan_u);
    fftw_execute(plan_c);
    fftw_destroy_plan(plan_u);
    fftw_destroy_plan(plan_c);
    
    // compute the mean of the distributions if needed later
    get_mean_signal(d);

    d->op->created_op = 1;
}//}}}

static
void prepare_op(all_data *d)
// does the necessary create calls
{//{{{
    printf("In onepoint.h -> prepare_op :\n");
    create_op(d);
}//}}}

void get_op(all_data *d, int Nbins, double *binedges, double *out, pdf_cl_uncl mode)
{//{{{
    prepare_op(d);
    
    double _binedges[Nbins+1];
    memcpy(_binedges, binedges, (Nbins+1) * sizeof(double));
    if (d->p->stype == kappa)
    {
        for (int ii=0; ii<=Nbins; ii++)
        {
            _binedges[ii] += ((mode==uncl) ? d->op->signalmeanu : d->op->signalmeanc);
        }
    }

    bin_1d(d->n->gr->Nsignal, d->n->gr->signalgrid,
           (mode==uncl) ? d->op->PDFu : d->op->PDFc,
           Nbins, _binedges, out, OPINTERP_TYPE);
}//}}}

