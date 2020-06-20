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

void null_halo_model(hmpdf_obj *d);
void reset_halo_model(hmpdf_obj *d);
double NFW_fundamental(hmpdf_obj *d, int z_index, int M_index, double *rs);
double Mconv(hmpdf_obj *d, int z_index, int M_index, hmpdf_mdef_e mdef_out, double *R, double *c);
void init_halo_model(hmpdf_obj *d);

#endif
