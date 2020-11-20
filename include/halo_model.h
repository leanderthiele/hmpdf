#ifndef HALO_MODEL_H
#define HALO_MODEL_H

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "hmpdf.h"

typedef struct//{{{
{
    int inited_halo;

    double *Duffy08_params;
    double *Tinker10_params;

    hmpdf_massfunc_corr_f massfunc_corr;
    void *massfunc_corr_params;

    hmpdf_conc_resc_f conc_resc;
    void *conc_resc_params;

    double **hmf;
    double **bias;

    gsl_spline *c_interp;
    gsl_interp_accel **c_accel;
}//}}}
halo_model_t;

int null_halo_model(hmpdf_obj *d);
int reset_halo_model(hmpdf_obj *d);
int NFW_fundamental(hmpdf_obj *d, int z_index, int M_index,
                    double mass_resc,
                    double *rhos, double *rs);
int Mconv(hmpdf_obj *d, int z_index, int M_index, hmpdf_mdef_e mdef_out,
          double mass_resc,
          double *M, double *R, double *c);
int init_halo_model(hmpdf_obj *d);

#endif
