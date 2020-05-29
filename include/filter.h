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
typedef double (*filter_fct)(void * /*all_data*/, double /*ell*/, filter_mode /*pdf or ps*/, int * /*z_index*/);

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

    ell_filter custom_ell;
    void *custom_ell_p;
    k_filter custom_k;
    void *custom_k_p;
}//}}}
filters_t;

void null_filters(all_data *d);
void reset_filters(all_data *d);
void apply_filters(all_data *d, int N, double *ell,
                   double *in, double *out, int stride,
                   filter_mode mode, int *z_index);
void init_filters(all_data *d);

#endif
