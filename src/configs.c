#include <limits.h>

#include "configs.h"

//                               A       B      C
static double
def_Duffy08_conc_params[] = { 5.71, -0.087, -0.47,   // M200c
                              7.85, -0.081, -0.71,   // Mvir
                             10.14, -0.081, -1.01, 0.0, 0.0, 0.0};// M200m

//                 param = factor * (1+z)^pwr
static double
def_Tinker10_hmf_params[] = { 0.589,  0.20,   // beta
                             -0.729, -0.08,   // phi
                             -0.243,  0.27,   // eta
                              0.864, -0.01,   // gamma
                              0.368,  0.00, };// alpha

//                      param = A0 * (M/1e14)^am * (1+z)^az
static double
def_Battaglia12_tsz_params[] = { 18.1  ,  0.154  , -0.758,   // P0
                                  0.497, -0.00865,  0.731,   // xc
                                  1.0  ,  0.0    ,  0.0  ,   // alpha
                                  4.35 ,  0.0393 ,  0.415,   // beta
                                 -0.3  ,  0.0    ,  0.0  , };// gamma

struct DEFAULTS def = { .Ncores={1,1,1000}, .verbosity=0, .warn_is_err=1,
                        .class_pre="none",
                        .Npoints_z={65,10,1000}, .z_min={0.0,0.0,6.0}, .z_max={6.0,0.1,10.0},
                        .dndz=NULL, .dndz_params=NULL,
                        .Npoints_M={65,10,1000}, .M_min={1e11,1e7,1e14}, .M_max={1e16,1e13,1e19},
                        .Npoints_signal={1024UL,32UL,10000UL},
                        .min_kappa={0.0,-10.0,0.0}, .min_tsz={0.0, -1e-2, 0.0},
                        .max_kappa={1.0,0.0,10.0}, .max_tsz={2e-4,1e-6,1e-2},
                        .Npoints_theta={500,10,10000}, .rout_scale={2.0,0.1,20.0},
                        .rout_rdef={hmpdf_mdef_v,0,hmpdf_mdef_m},
                        .pixel_sidelength={-1.0,0.0,100.0},
                        .tophat_radius={-1.0,0.0,100.0},
                        .gaussian_fwhm={-1.0,0.0,100.0},
                        .custom_ell_filter=NULL, .custom_ell_filter_params=NULL,
                        .custom_k_filter=NULL, .custom_k_filter_params=NULL,
                        .massfunc_corr=NULL, .massfunc_corr_params=NULL,
                        .mass_resc=NULL, .mass_resc_params=NULL,
                        .conc_resc=NULL, .conc_resc_params=NULL,
                        .mass_cuts=NULL, .mass_cuts_params=NULL,
                        .bias_resc=NULL, .bias_resc_params=NULL,
                        .Arico20_Nz=1, .Arico20_z=NULL, .Arico20_params=NULL,
                        .profiles_N=0, .profiles_fnames=NULL, .profiles_where=NULL, .profiles_Nr=0, .profiles_r=NULL,
                        .tot_profiles_N=0, .tot_profiles_fnames=NULL, .tot_profiles_where=NULL, .tot_profiles_Nr=0, .tot_profiles_r=NULL,
                        .DM_conc_params=NULL, .bar_conc_params=NULL,
                        .Nphi={1000,50,50000}, .phimax={150.0,10.0,1000.0},
                        .pixelexactmax={20,3,50},
                        .phijitter={0.02,1e-10,1e0}, .phipwr=2.0,
                        /* the integr_types have little influence,
                         * but can be varied to assess precision of integrations */
                        .zintegr_type={hmpdf_legendre,0,hmpdf_chebyshev2},
                        .zintegr_alpha=0.0, .zintegr_beta=0.0,
                        .Mintegr_type={hmpdf_legendre,0,hmpdf_chebyshev2},
                        .Mintegr_alpha=0.0, .Mintegr_beta=0.0,
                        .Duffy08_p=def_Duffy08_conc_params,
                        .Tinker10_p=def_Tinker10_hmf_params,
                        .Battaglia12_p=def_Battaglia12_tsz_params,
                        .noise_pwr=NULL, .noise_pwr_params=NULL,
                        .fsky={-1.0,0.0,1.0}, .pxlgrid={3,1,20}, .mappoisson=1, .mapseed=INT_MAX};

// The following is only needed for more reliable interaction
//     with the python wrapper

#define CHECK_EQU_SIZE(dtype)                         \
    do {                                              \
        if (sizeof(hmpdf_configs_e) == sizeof(dtype)) \
        {                                             \
            sprintf(out, #dtype);                     \
            return 0;                                 \
        }                                             \
    } while (0)

int enum_size_for_py(char *out)
// a simple helper function to find out which
//     C data type has the same size as an enum
//     on this platform,
//     so the python wrapper can pass objects of
//     the correct size
{
    CHECK_EQU_SIZE(short);
    CHECK_EQU_SIZE(int); 
    CHECK_EQU_SIZE(long);
    CHECK_EQU_SIZE(long long);

    return 1;
}

#undef CHECK_EQU_SIZE
