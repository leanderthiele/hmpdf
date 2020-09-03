#ifndef NOISE_H
#define NOISE_H

#include <complex.h>

#include <fftw3.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "utils.h"

#include "hmpdf.h"

typedef struct
{//{{{
    hmpdf_noise_pwr_f noise_pwr;
    void *noise_pwr_params;

    double sigmasq;

    int created_noise_zeta_interp;
    gsl_spline *zeta_interp;
    int Nzeta_accel;
    gsl_interp_accel **zeta_accel;

    double **conv_buffer_real; // [ (Nsignal_noisy+2) * Nsignal_noisy ]
    double complex **conv_buffer_comp; // not malloced
    fftw_plan **pconv_r2c; // conv_buffer_real -> conv_buffer_comp
    fftw_plan **pconv_c2r; // conv_buffer_comp -> conv_buffer_real

    long len_kernel;
    double *toepl;
}//}}}
noise_t;

int null_noise(hmpdf_obj *d);
int reset_noise(hmpdf_obj *d);
int init_noise(hmpdf_obj *d);

int noise_vect(hmpdf_obj *d, double *in, double *out);
int noise_matr(hmpdf_obj *d, double *in, double *out, int is_buffered, double phi);

#endif
