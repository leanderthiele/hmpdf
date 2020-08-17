#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_result.h>

#include "configs.h"
#include "utils.h"
#include "object.h"
#include "filter.h"

int
null_filters(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->f->inited_filters = 0;
    d->f->ffilters = NULL;
    d->f->z_dependent = NULL;
    d->f->quadraticpixel_interp = NULL;
    d->f->quadraticpixel_accel = NULL;
    d->f->quadraticpixel_ellmin = NULL;
    d->f->quadraticpixel_ellmax = NULL;

    ENDFCT
}//}}}

int
reset_filters(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\treset_filters\n");

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

    ENDFCT
}//}}}

static double
sinc(double x)
{//{{{
    if (x > 1e-4)
    {
        return sin(x) / x;
    }
    else
    {
        return 1.0-gsl_pow_2(x)/6.0+gsl_pow_4(x)/120.0-gsl_pow_6(x)/5040.0;
    }
}//}}}

static double
Bell_pdf(double phi, void *params)
{//{{{
    #ifdef LOGELL
    double ell = exp(*(double *)params);
    #else
    double ell = *(double *)params;
    #endif
    return sinc(0.5*ell*cos(phi)) * sinc(0.5*ell*sin(phi));
}//}}}

static double
Bell_ps(double phi, void *params)
{//{{{
    return gsl_pow_2(Bell_pdf(phi, params));
}//}}}

static double
(*Bell[])(double, void *) = {Bell_pdf, Bell_ps};

static int
_quadraticpixelinterp(hmpdf_obj *d, filter_mode mode)
{//{{{
    STARTFCT

    int Nell = PRWINDOW_INTERP_NELL;
    double *ell;
    double *Well;
    SAFEALLOC(ell,  malloc(Nell * sizeof(double)));
    SAFEALLOC(Well, malloc(Nell * sizeof(double)));

    #ifdef LOGELL
    SAFEHMPDF(linspace(Nell, log(PRWINDOW_INTERP_ELLMIN), log(PRWINDOW_INTERP_ELLMAX), ell));
    #else
    SAFEHMPDF(linspace(Nell, PRWINDOW_INTERP_ELLMIN, PRWINDOW_INTERP_ELLMAX, ell));
    #endif

    gsl_function integrand;
    integrand.function = Bell[mode];
    gsl_integration_workspace *ws;
    SAFEALLOC(ws, gsl_integration_workspace_alloc(PRWINDOW_INTEGR_LIMIT));

    for (int ii=0; ii<Nell; ii++)
    {
        integrand.params = ell+ii;
        double err;
        SAFEGSL(gsl_integration_qag(&integrand, 0.0, M_PI_4,
                                    PRWINDOW_INTEGR_EPSABS,
                                    PRWINDOW_INTEGR_EPSREL,
                                    PRWINDOW_INTEGR_LIMIT,
                                    PRWINDOW_INTEGR_KEY,
                                    ws, Well+ii, &err));
        Well[ii] *= 4.0 * M_1_PI;
    }

    gsl_integration_workspace_free(ws);

    SAFEALLOC(d->f->quadraticpixel_interp[mode],
              gsl_spline_alloc(gsl_interp_cspline, Nell));
    SAFEALLOC(d->f->quadraticpixel_accel[mode],
              malloc(d->Ncores * sizeof(gsl_interp_accel *)));
    for (int ii=0; ii<d->Ncores; ii++)
    {
        SAFEALLOC(d->f->quadraticpixel_accel[mode][ii],
                  gsl_interp_accel_alloc());
    }
    SAFEGSL(gsl_spline_init(d->f->quadraticpixel_interp[mode],
                            ell, Well, Nell));

    #ifdef LOGELL
    d->f->quadraticpixel_ellmin[mode] = exp(ell[0]);
    d->f->quadraticpixel_ellmax[mode] = exp(ell[Nell-1]);
    #else
    d->f->quadraticpixel_ellmin[mode] = ell[0];
    d->f->quadraticpixel_ellmax[mode] = ell[Nell-1];
    #endif

    free(ell);
    free(Well);

    ENDFCT
}//}}}

static int
_tophat(double x, double *out)
{//{{{
    STARTFCT

    if (x > 1e-4)
    {
        gsl_sf_result result;
        SAFEGSL(gsl_sf_bessel_J1_e(x, &result));
        *out = 2.0 * result.val / x;
    }
    else if (x > 0.0)
    {
        *out = 1.0 - gsl_pow_2(x)/8.0 + gsl_pow_4(x)/192.0 - gsl_pow_6(x)/9216.0;
    }
    else
    {
        *out = 0.0;
        HMPDFERR("argument negative");
    }

    ENDFCT
}//}}}

static int
_tophatsq(double x, double *out)
{//{{{
    STARTFCT

    if (x > 1e-4)
    {
        gsl_sf_result result;
        SAFEGSL(gsl_sf_bessel_J1_e(x, &result));
        *out = 4.0 * gsl_pow_2(result.val / x);
    }
    else if (x > 0.0)
    {
        *out = 1.0 - gsl_pow_2(x)/4.0 + 5.0*gsl_pow_4(x)/192.0 - 7.0*gsl_pow_6(x)/7608.0;
    }
    else
    {
        *out = 0.0;
        HMPDFERR("argument negative");
    }

    ENDFCT
}//}}}

static int
filter_quadraticpixel(void *d, double ell, filter_mode m, int *discard, double *out)
// assumes called with the physical reci_theta, i.e. reci_theta = j_n_0/theta_out
{//{{{
    STARTFCT

    hmpdf_obj *_d = (hmpdf_obj *)d;
    // rescale ell to unit half pixel sidelength
    ell *= 0.5 * _d->f->pixelside;
    if (ell < _d->f->quadraticpixel_ellmin[m])
    {
        *out = 1.0;
    }
    else if (ell > _d->f->quadraticpixel_ellmax[m])
    {
        *out = 0.0;
    }
    else
    {
        #ifdef LOGELL
        ell = log(ell);
        #endif
        SAFEGSL(gsl_spline_eval_e(_d->f->quadraticpixel_interp[m], ell,
                                  _d->f->quadraticpixel_accel[m][this_core()],
                                  out));
    }

    ENDFCT
}//}}}

static int
filter_tophat(void *d, double ell, filter_mode m, int *discard, double *out)
{//{{{
    STARTFCT

    hmpdf_obj *_d = (hmpdf_obj *)d;
    ell *= _d->f->tophat_radius;

    switch (m)
    {
        case filter_pdf : SAFEHMPDF(_tophat(ell, out));
                          break;
        case filter_ps  : SAFEHMPDF(_tophatsq(ell, out));
                          break;
        default         : *out = 0.0; // to avoid maybe-uninitialized 
                          HMPDFERR("unknown filter mode.");
    }

    ENDFCT
}//}}}

static int
filter_gaussian(void *d, double ell, filter_mode m, int *discard, double *out)
{//{{{
    STARTFCT

    hmpdf_obj *_d = (hmpdf_obj *)d;
    ell *= _d->f->gaussian_sigma;

    switch (m)
    {
        case filter_pdf : *out = exp(-0.5*ell*ell);
                          break;
        case filter_ps  : *out = exp(-ell*ell);
                          break;
        default         : *out = 0.0; // to avoid maybe-uninitialized
                          HMPDFERR("Unknown filter mode.");
    }

    ENDFCT
}//}}}

static int
filter_custom_ell(void *d, double ell, filter_mode m, int *discard, double *out)
{//{{{
    STARTFCT

    hmpdf_obj *_d = (hmpdf_obj *)d;
    double w =  _d->f->custom_ell(ell, _d->f->custom_ell_p);
    switch (m)
    {
        case filter_pdf : *out = w;
                          break;
        case filter_ps  : *out = w*w;
                          break;
        default         : *out = 0.0; // to avoid maybe-uninitialized
                          HMPDFERR("Unkown filter mode.");
    }

    ENDFCT
}//}}}

static int
filter_custom_k(void *d, double ell, filter_mode m, int *z_index, double *out)
{//{{{
    STARTFCT
    
    hmpdf_obj *_d = (hmpdf_obj *)d;
    // figure out comoving wavenumber corresponding to ell
    double k = ell / _d->c->comoving[*z_index];
    double w =  _d->f->custom_k(k, _d->n->zgrid[*z_index], _d->f->custom_k_p);
    switch (m)
    {
        case filter_pdf : *out = w;
                          break;
        case filter_ps  : *out = w*w;
                          break;
        default         : *out = 0.0; // to avoid maybe-uninitialized
                          HMPDFERR("Unknown filter mode.");
    }

    ENDFCT
}//}}}

int
apply_filters(hmpdf_obj *d, int N, double *ell, double *in, double *out,
              int stride, filter_mode mode, int *z_index)
// (1) if *z_index == NULL, applies only the z-independent filters
// (2) if *z_index != NULL and mode == filter_pdf, apply all filters
// (3)                     and mode == filter_ps,  apply only the z-dependent filters
// NOTE : function is safe if in and out point to same memory
{//{{{
    STARTFCT

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
                double temp;
                SAFEHMPDF((*(d->f->ffilters[ii]))((void *)d, ell[jj], mode, z_index, &temp));
                out[jj*stride] *= temp;
            }
        }
    }

    ENDFCT
}//}}}

int
init_filters(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->f->inited_filters) { return 0; }

    HMPDFPRINT(1, "init_filters\n");

    SAFEALLOC(d->f->ffilters,    malloc(10 * sizeof(filter_fct)));
    SAFEALLOC(d->f->z_dependent, malloc(10 * sizeof(int)));
    d->f->Nfilters = 0;

    if (d->f->pixelside > 0.0)
    {//{{{
        HMPDFPRINT(2, "\twill apply pixel filter\n");

        SAFEALLOC(d->f->quadraticpixel_interp,
                  malloc(2 * sizeof(gsl_spline *)));
        SAFEALLOC(d->f->quadraticpixel_accel,
                  malloc(2 * sizeof(gsl_interp_accel **)));
        SAFEALLOC(d->f->quadraticpixel_ellmin,
                  malloc(2 * sizeof(double)));
        SAFEALLOC(d->f->quadraticpixel_ellmax,
                  malloc(2 * sizeof(double)));
        SAFEHMPDF(_quadraticpixelinterp(d, filter_pdf));
        SAFEHMPDF(_quadraticpixelinterp(d, filter_ps));
        d->f->ffilters[d->f->Nfilters] = &filter_quadraticpixel;
        d->f->z_dependent[d->f->Nfilters] = 0;
        ++d->f->Nfilters;
    }//}}}
    if (d->f->tophat_radius > 0.0)
    {//{{{
        HMPDFPRINT(2, "\twill apply tophat filter\n");

        d->f->ffilters[d->f->Nfilters] = &filter_tophat;
        d->f->z_dependent[d->f->Nfilters] = 0;
        ++d->f->Nfilters;
    }//}}}
    if (d->f->gaussian_sigma > 0.0)
    {//{{{
        HMPDFPRINT(2, "\twill apply gaussian filter\n");

        d->f->ffilters[d->f->Nfilters] = &filter_gaussian;
        d->f->z_dependent[d->f->Nfilters] = 0;
        ++d->f->Nfilters;
    }//}}}
    if (d->f->custom_ell != NULL)
    {//{{{
        HMPDFPRINT(2, "\twill apply user-supplied ell-space filter\n");

        d->f->ffilters[d->f->Nfilters] = &filter_custom_ell;
        d->f->z_dependent[d->f->Nfilters] = 0;
        ++d->f->Nfilters;
    }//}}}
    if (d->f->custom_k != NULL)
    {//{{{
        HMPDFPRINT(2, "\twill apply user-supplied k-space filter\n");

        d->f->ffilters[d->f->Nfilters] = &filter_custom_k;
        d->f->z_dependent[d->f->Nfilters] = 1;
        ++d->f->Nfilters;
    }//}}}

    d->f->inited_filters = 1;

    ENDFCT
}//}}}

