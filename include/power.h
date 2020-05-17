#ifndef POWER_H
#define POWER_H

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "utils.h"
#include "hmpdf.h"

typedef struct//{{{
{
    int inited_power;

    double *k_arr;
    double *Pk_arr;
    interp1d *Pk_interp;

    double **ssq;
    double autocorr;

    int created_corr;
    gsl_spline *corr_interp;
    int Ncorr_accel;
    gsl_interp_accel **corr_accel;
}//}}}
power_t;

void null_power(all_data *d);
void reset_power(all_data *d);
void create_corr(all_data *d);
double Pk_linear(all_data *d, double k/*if LOGK is defined, this is log(k)*/);
double corr(all_data *d, int z_index, double phi);
void init_power(all_data *d);

#endif
