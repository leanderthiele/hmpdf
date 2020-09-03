#ifndef NOISE_H
#define NOISE_H

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "utils.h"

#include "hmpdf.h"

typedef struct
{//{{{
    hmpdf_noise_pwr_f noise_pwr;
    void *noise_pwr_params;

    double sigmasq;
    int created_zeta_interp;
    gsl_spline *zeta_interp;
    int Nzeta_accel;
    gsl_interp_accel **zeta_accel;

    long len_kernel;
    double *toepl;
}//}}}
noise_t;

int null_noise(hmpdf_obj *d);
int reset_noise(hmpdf_obj *d);
int init_noise(hmpdf_obj *d);

int noise_vect(hmpdf_obj *d, double *in, double *out);
int noise_matr(hmpdf_obj *d, double *in, double *out);

#endif
