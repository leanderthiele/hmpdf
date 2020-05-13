#include "configs.h"

                              //    A       B      C
double def_Duffy08_conc_params[] = { 5.71, -0.087, -0.47,   // M200c
                                     7.85, -0.081, -0.71,   // Mvir
                                    10.14, -0.081, -1.01, };// M200m

                       // param = factor x (1+x)^pwr
double def_Tinker10_hmf_params[] = { 0.589,  0.20,   // beta
                                    -0.729, -0.08,   // phi
                                    -0.243,  0.27,   // eta
                                     0.864, -0.01,   // gamma
                                     0.368,  0.00, };// alpha

                                    //  A0     alpham    alphaz
double def_Battaglia12_tsz_params[] = { 18.1  ,  0.154  , -0.758,   // P0
                                         0.497, -0.00865,  0.731,   // xc
                                         1.0  ,  0.0    ,  0.0  ,   // alpha
                                         4.35 ,  0.0393 ,  0.415,   // beta
                                        -0.3  ,  0.0    ,  0.0  , };// gamma

struct DEFAULTS def = { .Ncores=1, .class_pre="none",
                        .Npoints_z=65, .z_min=0.005, .z_max=6.0, .z_source=-1.0,
                        .Npoints_M=65, .M_min=1e11, .M_max=1e16,
                        .Npoints_signal=1024, .signal_min=0.0, .max_kappa=1.0, .max_tsz=2e-4,
                        .Npoints_theta=200, .rout_scale=2.0, .rout_rdef=mdef_v,
                        .pixel_sidelength=-1.0, .tophat_radius=-1.0, .gaussian_fwhm=-1.0,
                        .custom_ell_filter=NULL, .custom_ell_filter_params=NULL,
                        .custom_k_filter=NULL, .custom_k_filter_params=NULL,
                        .Nphi=10000, .phimax=150.0, .pixelexactmax=20,
                        .phijitter=0.02, .phipwr=2.0, .regularize_tp=0,
                        .monotonize=1,
                        /* the integr_types have little influence,
                         * but can be varied to assess precision of integrations */
                        .zintegr_type=legendre, .zintegr_alpha=0.0, .zintegr_beta=0.0,
                        .Mintegr_type=legendre, .Mintegr_alpha=0.0, .Mintegr_beta=0.0,
                        .Duffy08_p=def_Duffy08_conc_params,
                        .Tinker10_p=def_Tinker10_hmf_params,
                        .Battaglia12_p=def_Battaglia12_tsz_params, };
