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
    d->op->created_rolled_op = 0;
    d->op->created_noisy_op = 0;
    d->op->PDFu = NULL;
    d->op->PDFc = NULL;
    d->op->PDFu_rolled = NULL;
    d->op->PDFc_rolled = NULL;
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
    if (d->op->PDFu_rolled != NULL) { free(d->op->PDFu_rolled); }
    if (d->op->PDFc_rolled != NULL) { free(d->op->PDFc_rolled); }
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

typedef struct
{//{{{
    double (*f)(double);
    double lambda;
    interp1d *s;
    int *hmpdf_status;
    gsl_integration_workspace *ws;
}//}}}
not_inv_integrand_s;

static double
not_inv_integrand(double t, void *params)
{//{{{
    not_inv_integrand_s *p = (not_inv_integrand_s *)params;
    double s;
    *(p->hmpdf_status) = interp1d_eval(p->s, t, &s);
    return t * p->f(s * p->lambda);
}//}}}

static int
not_inv_integral(hmpdf_obj *d, not_inv_integrand_s *params, int lambda_index, complex *out)
// TODO can potentially speed this up a bit
{//{{{
    STARTFCT

    gsl_function f;
    f.function = not_inv_integrand;
    f.params = params;
    params->lambda = d->n->lambdagrid[lambda_index];
    params->hmpdf_status = &hmpdf_status;
    double integrals[2];
    double err;
    size_t neval;
    for (int ii=0; ii<2; ii++)
    {
        params->f = (ii==0) ? &cos : &sin;
//        SAFEGSL(gsl_integration_qng(&f, 0.0, 1.0, 1e-6, 1e-2, integrals+ii, &err, &neval))
        SAFEGSL(gsl_integration_qag(&f, 0.0, 1.0, 1e-6, 1e-2, 1000, 1, params->ws, integrals+ii, &err))
    }
    *out = integrals[0] - _Complex_I * integrals[1];

    ENDFCT
}//}}}

static int
lambda_loop(hmpdf_obj *d, int z_index, int M_index, complex *out)
{//{{{
    STARTFCT

    fprintf(stdout, "%d\t%d\n", z_index, M_index);
    fflush(stdout);

    // reverse the signal profile
    SAFEALLOC(double *, temp, malloc((d->p->Ntheta+1) * sizeof(double)))
    reverse(d->p->Ntheta+1, d->p->profiles[z_index][M_index]+1, temp);
    SAFEALLOC(not_inv_integrand_s *, params,
              malloc(d->Ncores * sizeof(not_inv_integrand_s)))
    for (int ii=0; ii<d->Ncores; ii++)
    {
        SAFEHMPDF(new_interp1d(d->p->Ntheta+1, d->p->incr_tgrid, temp, temp[0], 0.0,
                               PRINTERP_TYPE, d->p->incr_tgrid_accel[ii], &(params[ii].s)))
        SAFEALLOC(, params[ii].ws, gsl_integration_workspace_alloc(1000))
    }

    // loop over lambda values
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(1) /*FIXME*/ schedule(dynamic)
    #endif
    for (int ii=0; ii<d->n->Nsignal/2+1; ii++)
    {
        if (temp[0] * d->n->lambdagrid[ii] > 1e1) { continue; }
        if (hmpdf_status) { continue; }
        SAFEHMPDF_NORETURN(not_inv_integral(d, params+this_core(), ii, out+ii))
        // FIXME
        if (temp[0] * d->n->lambdagrid[ii] > 8e0)
        {
            printf("%.8e\n", temp[0] * d->n->lambdagrid[ii]);
            int N = 10000;
            double *tgrid = malloc(N * sizeof(double));
            double *y = malloc(N * sizeof(double));
            double *s = malloc(N * sizeof(double));
            linspace(N, 0.0, 1.0, tgrid);
            double maxy = 0.0;
            for (int jj=0; jj<N; jj++)
            {
                y[jj] = not_inv_integrand(tgrid[jj], params+this_core());
                interp1d_eval(params[this_core()].s, tgrid[jj], s+jj);
                s[jj] /= temp[0];
                maxy = GSL_MAX(maxy, y[jj]);
            }
            for (int jj=0; jj<N; jj++)
            {
                y[jj] /= maxy;
            }
            gnuplot *gp = plot(NULL, N, tgrid, y);
            plot(gp, N, tgrid, s);
            
            show(gp);
        }
        out[ii] *= 2.0 * M_PI * gsl_pow_2(d->p->profiles[z_index][M_index][0]);
    }

    free(temp);
    for (int ii=0; ii<d->Ncores; ii++)
    {
        delete_interp1d(params[ii].s);
        gsl_integration_workspace_free(params[ii].ws);
    }
    free(params);

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
                    d->op->PDFu, &(d->op->signalmeanu)))
    SAFEHMPDF(_mean(d->n->Nsignal, d->n->signalgrid,
                    d->op->PDFc, &(d->op->signalmeanc)))

    ENDFCT
}//}}}

static int
roll_op(int N, double *xorig, double *xnew, double *yorig, double *ynew)
// assumes that min(xorig) = 0.0
{//{{{
    STARTFCT
    // interpolate the y-vector
    interp1d *interp;
    SAFEHMPDF(new_interp1d(N, xorig, yorig, 0.0, 0.0, ROLLOP_INTERP, NULL, &interp))
    // roll the array
    for (int ii=0; ii<N; ii++)
    {
        double x = xnew[ii];
        if (x < 0.0)
        {
            x += xorig[N-1];
        }
        SAFEHMPDF(interp1d_eval(interp, x, ynew+ii))
    }
    // destroy the interpolator
    delete_interp1d(interp);

    ENDFCT
}//}}}

int
create_noisy_op(hmpdf_obj *d)
// convolves the original PDF with a Gaussian kernel of width sigma = noise
{//{{{
    STARTFCT

    if (d->op->created_noisy_op) { return hmpdf_status; }

    HMPDFPRINT(2, "\tcreate_noisy_op\n")

    double *in[] = {d->op->PDFu_rolled, d->op->PDFc_rolled};
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

    // compute the mean of the distributions
    SAFEHMPDF(get_mean_signal(d))

    d->op->created_op = 1;

    ENDFCT
}//}}}

int
create_rolled_op(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->op->created_rolled_op) { return hmpdf_status; }

    HMPDFPRINT(2, "\tcreate_rolled_op\n")

    SAFEALLOC(, d->op->PDFu_rolled, malloc(d->n->Nsignal * sizeof(double)))
    SAFEALLOC(, d->op->PDFc_rolled, malloc(d->n->Nsignal * sizeof(double)))

    SAFEHMPDF(roll_op(d->n->Nsignal, d->n->signalgrid, d->n->user_signalgrid,
                      d->op->PDFu, d->op->PDFu_rolled))
    SAFEHMPDF(roll_op(d->n->Nsignal, d->n->signalgrid, d->n->user_signalgrid,
                      d->op->PDFc, d->op->PDFc_rolled))

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
    SAFEHMPDF(create_rolled_op(d))
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
                     (noisy) ? d->n->signalgrid_noisy : d->n->user_signalgrid,
                     (noisy) ? ((incl_2h) ? d->op->PDFc_noisy : d->op->PDFu_noisy)
                     : ((incl_2h) ? d->op->PDFc_rolled : d->op->PDFu_rolled),
                     Nbins, _binedges, op, OPINTERP_TYPE))

    ENDFCT
}//}}}

