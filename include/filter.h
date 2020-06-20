#ifndef FILTER_H
#define FILTER_H

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "hmpdf.h"

typedef enum//{{{
{
    filter_pdf,
    filter_ps,
    filter_end,
}//}}}
filter_mode;

typedef double (*filter_fct)(void * /*hmpdf_obj*/, double /*ell*/,
                             filter_mode /*pdf or ps*/, int * /*z_index*/);

typedef struct//{{{
{
    int inited_filters;

    int Nfilters;
    int *z_dependent;
    filter_fct *ffilters;

    double pixelside;
    gsl_spline **quadraticpixel_interp;
    gsl_interp_accel ***quadraticpixel_accel;
    double *quadraticpixel_ellmin;
    double *quadraticpixel_ellmax;

    double tophat_radius;
    double gaussian_sigma;

    hmpdf_ell_filter_f custom_ell;
    void *custom_ell_p;
    hmpdf_k_filter_f custom_k;
    void *custom_k_p;
}//}}}
filters_t;

void null_filters(hmpdf_obj *d);
void reset_filters(hmpdf_obj *d);
void apply_filters(hmpdf_obj *d, int N, double *ell,
                   double *in, double *out, int stride,
                   filter_mode mode, int *z_index);
void init_filters(hmpdf_obj *d);

#endif
