#ifndef POWER_H
#define POWER_H

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "utils.h"
#include "hmpdf.h"

typedef struct//{{{
{
    int inited_power;

    double **ssq;
    double autocorr;

    int created_corr;
    double corr_rmax;
    gsl_spline *corr_interp;
    int Ncorr_accel;
    gsl_interp_accel **corr_accel;
}//}}}
power_t;

int null_power(hmpdf_obj *d);
int reset_power(hmpdf_obj *d);
int create_corr(hmpdf_obj *d);
int Pk_linear(hmpdf_obj *d, double k/*if LOGK is defined, this is log(k)*/, double *out);
int corr(hmpdf_obj *d, int z_index, double phi, double *out);
int init_power(hmpdf_obj *d);

#endif
