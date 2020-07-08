#ifndef PROFILES_H
#define PROFILES_H

#include <gsl/gsl_interp.h>
#include <gsl/gsl_dht.h>

#include "hmpdf.h"

typedef struct//{{{
{
    int inited_profiles;

    double *Battaglia12_params;

    hmpdf_signaltype_e stype;

    double rout_scale;
    int rout_def;
    int Ntheta;
    double *decr_tgrid;
    double *incr_tgrid;
    double *decr_tsqgrid;
    double *incr_tsqgrid;
    double *reci_tgrid; // reciprocal space grid

    gsl_interp_accel **incr_tgrid_accel;
    gsl_interp_accel *reci_tgrid_accel;

    double ***profiles; // each profile has as zero entry theta out and then the profile

    int created_conj_profiles;
    double ***conj_profiles; // each profile has as zero entry the rescaling such that reci_thetagrid -> ell

    int created_filtered_profiles;

    int prtilde_Ntheta;
    double *prtilde_thetagrid;

    int created_monotonicity;
    int **is_not_monotonic;

    int created_segments;
    int ***segment_boundaries;

    gsl_dht *dht_ws;
}//}}}
profiles_t;

typedef enum
{
    dtsq_of_s,
    t_of_s,
}
inv_profile_e;

int null_profiles(hmpdf_obj *d);
int reset_profiles(hmpdf_obj *d);
int init_profiles(hmpdf_obj *d);
int create_conj_profiles(hmpdf_obj *d);
int create_filtered_profiles(hmpdf_obj *d);
int create_monotonicity(hmpdf_obj *d);
int create_segments(hmpdf_obj *d);

int s_of_t(hmpdf_obj *d, int z_index, int M_index, int Nt, double *t, double *s);
int s_of_ell(hmpdf_obj *d, int z_index, int M_index, int Nell, double *ell, double *s);
int inv_profile(hmpdf_obj *d, int z_index, int M_index, int segment, int *filled, inv_profile_e mode, double *out);

#endif
