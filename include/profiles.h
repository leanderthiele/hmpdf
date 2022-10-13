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
    double ***filtered_profiles;

    int created_segments;
    int ***segment_boundaries;

    gsl_dht *dht_ws;

    hmpdf_mass_resc_f mass_resc;
    void *mass_resc_params;

    // for printing profiles to file
    int tot_profiles_N;
    char **tot_profiles_fnames;
    double *tot_profiles_where;
    int tot_profiles_Nr;
    double *tot_profiles_r;

    // populated later
    int *tot_profiles_indices;
}//}}}
profiles_t;

typedef enum
{//{{{
    dtsq_of_s,
    t_of_s,
}//}}}
inv_profile_e;

typedef struct
{//{{{
    long start; // the start index in the signal grid
    long len;   // length of this batch
    int incr;     // +-1 --> loop over signal grid such that
                  //    theta is always decreasing
    double *data; // either t_of_s or dtsq_of_s, of length len
}//}}}
batch_t;

typedef struct
{//{{{
    int Nbatches;
    batch_t *batches;
}//}}}
batch_container_t;

void delete_batch(batch_t *b);

int null_profiles(hmpdf_obj *d);
int reset_profiles(hmpdf_obj *d);
int init_profiles(hmpdf_obj *d);
int create_conj_profiles(hmpdf_obj *d);
int create_filtered_profiles(hmpdf_obj *d);
int create_segments(hmpdf_obj *d);

int s_of_t(hmpdf_obj *d, int z_index, int M_index, long Nt, double *t, double *s);
int s_of_ell(hmpdf_obj *d, int z_index, int M_index, int Nell, double *ell, double *s);
int inv_profile(hmpdf_obj *d, int z_index, int M_index, int segment,
                inv_profile_e mode, batch_t *b);

#endif
