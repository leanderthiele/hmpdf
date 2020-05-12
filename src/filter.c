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

char filter_pdf_ps[][256] = {"pdf", "ps"};

typedef struct//{{{
{
    all_data *d;
    filter_mode mode;
}//}}}
filter_params;

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
// FIXME not in agreement with python calculation!
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
    d->f->quadraticpixel_accel[mode] = gsl_interp_accel_alloc();
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
double filter_quadraticpixel(double ell, void *params)
// assumes called with the physical reci_theta, i.e. reci_theta = j_n_0/theta_out
{//{{{
    filter_params *p = (filter_params *)params;
    // rescale ell to unit half pixel sidelenght
    ell *= 0.5 * p->d->f->pixelside;
    if (ell < p->d->f->quadraticpixel_ellmin[p->mode])
    {
        return 1.0;
    }
    else if (ell > p->d->f->quadraticpixel_ellmax[p->mode])
    {
        return 0.0;
    }
    else
    {
        #ifdef LOGELL
        ell = log(ell);
        #endif
        return gsl_spline_eval(p->d->f->quadraticpixel_interp[p->mode], ell,
                               p->d->f->quadraticpixel_accel[p->mode]);
    }
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
double filter_tophat(double ell, void *params)
{//{{{
    filter_params *p = (filter_params *)params;
    ell *= p->d->f->tophat_radius;

    switch (p->mode)
    {
        case filter_pdf : return _tophat(ell);
        case filter_ps  : return _tophatsq(ell);
        default         : printf("Unknown filter mode in filter_tophat\n");
                          return 0.0;
    }
}//}}}

static
double filter_gaussian(double ell, void *params)
{
    filter_params *p = (filter_params *)params;
    ell *= p->d->f->gaussian_sigma;

    switch (p->mode)
    {
        case filter_pdf : return exp(-0.5*ell*ell);
        case filter_ps  : return exp(-ell*ell);
        default         : printf("Unknown filter mode in filter_gaussian\n");
                          return 0.0;
    }
}

void apply_filters(all_data *d, int N, double *ell, double *in, double *out, filter_mode mode)
{//{{{
    if (in != out)
    {
        memcpy(out, in, N * sizeof(double));
    }
    for (int ii=0; ii<d->f->Nfilters; ii++)
    {
        filter_params *p = (filter_params *)d->f->ffilters[ii].params;
        p->mode = mode;
        for (int jj=0; jj<N; jj++)
        {
            out[jj] *= GSL_FN_EVAL(d->f->ffilters+ii, ell[jj]);
        }
    }

}//}}}

void init_filters(all_data *d)
{//{{{
    if (d->f->inited_filters) { return; }
    printf("In filter.h -> init_filters.\n");
    d->f->ffilters = (gsl_function *)malloc(10 * sizeof(gsl_function));
    d->f->Nfilters = 0;

    if (d->f->pixelside > 0.0)
    {
        printf("\twill apply quadratic pixel filter\n");
        filter_params *p = (filter_params *)malloc(sizeof(filter_params));
        p->d = d;
        p->mode = filter_end;
    
        d->f->quadraticpixel_interp = (gsl_spline **)malloc(2 * sizeof(gsl_spline *));
        d->f->quadraticpixel_accel = (gsl_interp_accel **)malloc(2 * sizeof(gsl_interp_accel *));
        d->f->quadraticpixel_ellmin = (double *)malloc(2 * sizeof(double));
        d->f->quadraticpixel_ellmax = (double *)malloc(2 * sizeof(double));
        _quadraticpixelinterp(d, filter_pdf);
        _quadraticpixelinterp(d, filter_ps);
        d->f->ffilters[d->f->Nfilters].function = filter_quadraticpixel;
        d->f->ffilters[d->f->Nfilters].params = p;
        ++d->f->Nfilters;
    }
    if (d->f->tophat_radius > 0.0)
    {
        printf("\twill apply tophat filter\n");
        filter_params *p = (filter_params *)malloc(sizeof(filter_params));
        p->d = d;
        p->mode = filter_end;

        d->f->ffilters[d->f->Nfilters].function = filter_tophat;
        d->f->ffilters[d->f->Nfilters].params = p;
        ++d->f->Nfilters;
    }
    if (d->f->gaussian_sigma > 0.0)
    {
        printf("\twill apply gaussian filter\n");
        filter_params *p = (filter_params *)malloc(sizeof(filter_params));
        p->d = d;
        p->mode = filter_end;

        d->f->ffilters[d->f->Nfilters].function = filter_gaussian;
        d->f->ffilters[d->f->Nfilters].params = p;
        ++d->f->Nfilters;
    }

    d->f->inited_filters = 1;
}//}}}

