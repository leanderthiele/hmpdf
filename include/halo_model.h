#ifndef HALO_MODEL_H
#define HALO_MODEL_H

#include <gsl/gsl_interp.h>

#include "hmpdf.h"

typedef struct//{{{
{
    int inited_halo;

    double *Duffy08_params;
    double *Tinker10_params;

    double **hmf;
    double **bias;

    gsl_spline *c_interp;
    gsl_interp_accel **c_accel;
}//}}}
halo_model_t;

void null_halo_model(all_data *d);
void reset_halo_model(all_data *d);
double NFW_fundamental(all_data *d, int z_index, int M_index, double *rs);
double Mconv(all_data *d, int z_index, int M_index, mdef mdef_out, double *R, double *c);
void init_halo_model(all_data *d);

#endif
