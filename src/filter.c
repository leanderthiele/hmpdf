#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>

#include "configs.h"
#include "utils.h"
#include "data.h"
#include "filter.h"

void null_filters(all_data *d)
{//{{{
    d->f->inited_filters = 0;
    d->f->ffilters = NULL;
    d->f->z_dependent = NULL;
    d->f->quadraticpixel_interp = NULL;
    d->f->quadraticpixel_accel = NULL;
    d->f->quadraticpixel_ellmin = NULL;
    d->f->quadraticpixel_ellmax = NULL;
}//}}}

void reset_filters(all_data *d)
{//{{{
    if (d->f->z_dependent != NULL) { free(d->f->z_dependent); }
    if (d->f->ffilters != NULL) { free(d->f->ffilters); }
    if (d->f->quadraticpixel_interp != NULL)
    {
        if (d->f->quadraticpixel_interp[0] != NULL) { gsl_spline_free(d->f->quadraticpixel_interp[0]); }
        if (d->f->quadraticpixel_interp[1] != NULL) { gsl_spline_free(d->f->quadraticpixel_interp[1]); }
        free(d->f->quadraticpixel_interp);
    }
    if (d->f->quadraticpixel_accel != NULL)
    {
        for (int ii=0; ii<2; ii++)
        {
            if (d->f->quadraticpixel_accel[ii] != NULL)
            {
                for (int jj=0; jj<d->Ncores; jj++)
                {
                    gsl_interp_accel_free(d->f->quadraticpixel_accel[ii][jj]);
                }
                free(d->f->quadraticpixel_accel[ii]);
            }
        }
        free(d->f->quadraticpixel_accel);
    }
    if (d->f->quadraticpixel_ellmin != NULL) { free(d->f->quadraticpixel_ellmin); }
    if (d->f->quadraticpixel_ellmax != NULL) { free(d->f->quadraticpixel_ellmax); }
}//}}}
    
char filter_pdf_ps[][256] = {"pdf", "ps"};

static
double sinc(double x)
{//{{{
    if (x > 1e-4)
    {
        return sin(x)/x;
    }
    else
    {
        return 1.0-gsl_pow_2(x)/6.0+gsl_pow_4(x)/120.0-gsl_pow_6(x)/5040.0;
    }
}//}}}

static
double Bell_pdf(double phi, void *params)
{//{{{
    #ifdef LOGELL
    double ell = exp(*(double *)params);
    #else
    double ell = *(double *)params;
    #endif
    return sinc(0.5*ell*cos(phi)) * sinc(0.5*ell*sin(phi));
}//}}}

static
double Bell_ps(double phi, void *params)
{//{{{
    return gsl_pow_2(Bell_pdf(phi, params));
}//}}}

static
double (*Bell[])(double, void *) = {Bell_pdf, Bell_ps};

static
void _quadraticpixelinterp(all_data *d, filter_mode mode)
{//{{{
    double **Well;
    char fname[512];
    #ifdef LOGELL
    char log_lin[] = "log";
    #else
    char log_lin[] = "lin";
    #endif
    sprintf(fname, "%sell_W_%s.bin", log_lin, filter_pdf_ps[mode]);
    int Nell;
    Well = fromfile(fname, &Nell, 2);
    if (Well == NULL)
    // File not found, need to interpolate
    {
        Nell = PRWINDOW_INTERP_NELL;
        Well = (double **)malloc(2 * sizeof(double *));
        Well[0] = (double *)malloc(Nell * sizeof(double));
        Well[1] = (double *)malloc(Nell * sizeof(double));

        #ifdef LOGELL
        linspace(Nell, log(PRWINDOW_INTERP_ELLMIN), log(PRWINDOW_INTERP_ELLMAX), Well[0]);
        #else
        linspace(Nell, PRWINDOW_INTERP_ELLMIN, PRWINDOW_INTERP_ELLMAX, Well[0]);
        #endif
        gsl_function integrand;
        integrand.function = Bell[mode];
        gsl_integration_workspace *ws = gsl_integration_workspace_alloc(PRWINDOW_INTEGR_LIMIT);
        for (int ii=0; ii<Nell; ii++)
        {
            integrand.params = Well[0]+ii;
            double err;
            gsl_integration_qag(&integrand, 0.0, M_PI_4,
                                PRWINDOW_INTEGR_EPSABS, PRWINDOW_INTEGR_EPSREL,
                                PRWINDOW_INTEGR_LIMIT, PRWINDOW_INTEGR_KEY,
                                ws, Well[1]+ii, &err);
            Well[1][ii] *= 4.0 * M_1_PI;
        }
        gsl_integration_workspace_free(ws);

        tofile(fname, Nell, 2, Well[0], Well[1]);
    }

    d->f->quadraticpixel_interp[mode] = gsl_spline_alloc(gsl_interp_cspline, Nell);
    d->f->quadraticpixel_accel[mode]
        = (gsl_interp_accel **)malloc(d->Ncores * sizeof(gsl_interp_accel *));
    for (int ii=0; ii<d->Ncores; ii++)
    {
        d->f->quadraticpixel_accel[mode][ii] = gsl_interp_accel_alloc();
    }
    gsl_spline_init(d->f->quadraticpixel_interp[mode], Well[0], Well[1], Nell);

    #ifdef LOGELL
    d->f->quadraticpixel_ellmin[mode] = exp(Well[0][0]);
    d->f->quadraticpixel_ellmax[mode] = exp(Well[0][Nell-1]);
    #else
    d->f->quadraticpixel_ellmin[mode] = Well[0][0];
    d->f->quadraticpixel_ellmax[mode] = Well[0][Nell-1];
    #endif

    free(Well[0]);
    free(Well[1]);
    free(Well);
}//}}}

static
double _tophat(double x)
{//{{{
    if (x > 1e-4)
    {
        return 2.0 * gsl_sf_bessel_J1(x) / x;
    }
    else
    {
        return 1.0 - gsl_pow_2(x)/8.0 + gsl_pow_4(x)/192.0 - gsl_pow_6(x)/9216.0;
    }
}//}}}

static double _tophatsq(double x)
{//{{{
    if (x > 1e-4)
    {
        return 4.0 * gsl_pow_2(gsl_sf_bessel_J1(x) / x);
    }
    else
    {
        return 1.0 - gsl_pow_2(x)/4.0 + 5.0*gsl_pow_4(x)/192.0 - 7.0*gsl_pow_6(x)/7608.0;
    }
}//}}}

static
double filter_quadraticpixel(void *d, double ell, filter_mode m, int *discard)
// assumes called with the physical reci_theta, i.e. reci_theta = j_n_0/theta_out
{//{{{
    all_data *_d = (all_data *)d;
    // rescale ell to unit half pixel sidelenght
    ell *= 0.5 * _d->f->pixelside;
    if (ell < _d->f->quadraticpixel_ellmin[m])
    {
        return 1.0;
    }
    else if (ell > _d->f->quadraticpixel_ellmax[m])
    {
        return 0.0;
    }
    else
    {
        #ifdef LOGELL
        ell = log(ell);
        #endif
        return gsl_spline_eval(_d->f->quadraticpixel_interp[m], ell,
                               _d->f->quadraticpixel_accel[m][this_core()]);
    }
}//}}}

static
double filter_tophat(void *d, double ell, filter_mode m, int *discard)
{//{{{
    all_data *_d = (all_data *)d;
    ell *= _d->f->tophat_radius;

    switch (m)
    {
        case filter_pdf : return _tophat(ell);
        case filter_ps  : return _tophatsq(ell);
        default         : fprintf(stdout, "Unknown filter mode in filter_tophat\n");
                          fflush(stdout);
                          return 0.0;
    }
}//}}}

static
double filter_gaussian(void *d, double ell, filter_mode m, int *discard)
{//{{{
    all_data *_d = (all_data *)d;
    ell *= _d->f->gaussian_sigma;

    switch (m)
    {
        case filter_pdf : return exp(-0.5*ell*ell);
        case filter_ps  : return exp(-ell*ell);
        default         : fprintf(stdout, "Unknown filter mode in filter_gaussian\n");
                          fflush(stdout);
                          return 0.0;
    }
}//}}}

static
double filter_custom_ell(void *d, double ell, filter_mode m, int *discard)
{//{{{
    all_data *_d = (all_data *)d;
    double w =  _d->f->custom_ell(ell, _d->f->custom_ell_p);
    switch (m)
    {
        case filter_pdf : return w;
        case filter_ps  : return w*w;
        default         : fprintf(stdout, "Unknown filter_mode in filter_custom_ell.\n");
                          fflush(stdout);
                          return 0.0;
    }
}//}}}

static
double filter_custom_k(void *d, double ell, filter_mode m, int *z_index)
{//{{{
    all_data *_d = (all_data *)d;
    // figure out comoving wavenumber corresponding to ell
    double k = ell / _d->c->comoving[*z_index];
    double w =  _d->f->custom_k(k, _d->n->zgrid[*z_index], _d->f->custom_k_p);
    switch (m)
    {
        case filter_pdf : return w;
        case filter_ps  : return w*w;
        default         : fprintf(stdout, "Unknown filter_mode in filter_custom_k.\n");
                          fflush(stdout);
                          return 0.0;
    }
}//}}}

void apply_filters(all_data *d, int N, double *ell, double *in, double *out, int stride, filter_mode mode, int *z_index)
// (1) if *z_index == NULL, applies only the z-independent filters
// (2) if *z_index != NULL and mode == filter_pdf, apply all filters
// (3)                     and mode == filter_ps,  apply only the z-dependent filters
// NOTE : function is safe if in and out point to same memory
{//{{{
    if (in != out)
    {
        memcpy(out, in, N * sizeof(double));
    }
    for (int ii=0; ii<d->f->Nfilters; ii++)
    {
        // check if we need to apply this filter
        if ((z_index == NULL && d->f->z_dependent[ii]) // case (1)
            || (z_index != NULL && mode == filter_ps && !(d->f->z_dependent[ii]))) // case (3)
        {
            continue;
        }
        else
        {
            for (int jj=0; jj<N; jj++)
            {
                out[jj*stride] *= (*(d->f->ffilters[ii]))((void *)d, ell[jj], mode, z_index);
            }
        }
    }

}//}}}

void init_filters(all_data *d)
{//{{{
    if (d->f->inited_filters) { return; }
    fprintf(stdout, "In filter.h -> init_filters.\n");
    fflush(stdout);
    d->f->ffilters = (filter_fct *)malloc(10 * sizeof(filter_fct));
    d->f->z_dependent = (int *)malloc(10 * sizeof(int));
    d->f->Nfilters = 0;

    if (d->f->pixelside > 0.0)
    {//{{{
        fprintf(stdout, "\twill apply quadratic pixel filter\n");
        fflush(stdout);
        d->f->quadraticpixel_interp = (gsl_spline **)malloc(2 * sizeof(gsl_spline *));
        d->f->quadraticpixel_accel = (gsl_interp_accel ***)malloc(2 * sizeof(gsl_interp_accel **));
        d->f->quadraticpixel_ellmin = (double *)malloc(2 * sizeof(double));
        d->f->quadraticpixel_ellmax = (double *)malloc(2 * sizeof(double));
        _quadraticpixelinterp(d, filter_pdf);
        _quadraticpixelinterp(d, filter_ps);
        d->f->ffilters[d->f->Nfilters] = &filter_quadraticpixel;
        d->f->z_dependent[d->f->Nfilters] = 0;
        ++d->f->Nfilters;
    }//}}}
    if (d->f->tophat_radius > 0.0)
    {//{{{
        fprintf(stdout, "\twill apply tophat filter\n");
        fflush(stdout);
        d->f->ffilters[d->f->Nfilters] = &filter_tophat;
        d->f->z_dependent[d->f->Nfilters] = 0;
        ++d->f->Nfilters;
    }//}}}
    if (d->f->gaussian_sigma > 0.0)
    {//{{{
        fprintf(stdout, "\twill apply gaussian filter\n");
        fflush(stdout);
        d->f->ffilters[d->f->Nfilters] = &filter_gaussian;
        d->f->z_dependent[d->f->Nfilters] = 0;
        ++d->f->Nfilters;
    }//}}}
    if (d->f->custom_ell != NULL)
    {//{{{
        fprintf(stdout, "\twill apply user-supplied ell-space filter\n");
        fflush(stdout);
        d->f->ffilters[d->f->Nfilters] = &filter_custom_ell;
        d->f->z_dependent[d->f->Nfilters] = 0;
        ++d->f->Nfilters;
    }//}}}
    if (d->f->custom_k != NULL)
    {//{{{
        fprintf(stdout, "\twill apply user-supplied k-space filter\n");
        fflush(stdout);
        d->f->ffilters[d->f->Nfilters] = &filter_custom_k;
        d->f->z_dependent[d->f->Nfilters] = 1;
        ++d->f->Nfilters;
    }//}}}

    d->f->inited_filters = 1;
}//}}}

