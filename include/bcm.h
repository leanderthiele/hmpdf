#ifndef BCM_H
#define BCM_H

#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>

#include "hmpdf.h"
#include "utils.h"

typedef struct//{{{
{
    // DMO NFW profile
    double M200c,
           R200c,
           rs,
           rhos;

    // bound gas properties
    double bg_y0,
           bg_r_inn, bg_r_out,
           bg_beta_i;
    gsl_integration_workspace *bg_integr_ws;

    // central galaxy properties
    double cg_y0,
           cg_Rh;

    // re-accreted gas properties
    double rg_y0,
           rg_sigma, rg_mu;

    // ejected gas properties
    double eg_rej,
           eg_f;

    // relaxed dark matter
    double dm_f;
    double *dm_xi;
    gsl_interp_accel *dm_r_accel;
    gsl_interp *dm_xi_interp;
}//}}}
bcm_ws;

typedef struct//{{{
{
    int inited_bcm;

    // parameterization according to 2009.14225
    // { eta, M_c, beta, theta_inn, theta_out, M_inn, M_r, M1_z0_cen, }
    double *Arico20_params;

    // for the relaxed dark matter profile interpolation
    int Nradii;
    int R200c_idx; // the index closest to R200c in the array
    double rmin, rmax; // in units of R200c
    double *radii; // where xi(r) will be evaluated

    bcm_ws **ws; // one for each thread
                 // probably want to have these somewhat separate in memory to avoid false sharing
}//}}}
bcm_t;

int null_bcm(hmpdf_obj *d);
int reset_bcm(hmpdf_obj *d);
int init_bcm(hmpdf_obj *d);

// this function does the initial mallocs
int bcm_new_ws(hmpdf_obj *d, bcm_ws *ws);

// this function can re-use previously allocated bcm_halo_t
int bcm_init_ws(hmpdf_obj *d, int z_index, int M_index, double mass_resc, bcm_ws *ws);

int bcm_delete_ws(bcm_ws *ws);

int bcm_density_profile(hmpdf_obj *d, bcm_ws *ws, double r, double *out);



#endif
