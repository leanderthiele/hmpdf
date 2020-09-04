#ifndef CONFIGS_H
#define CONFIGS_H

#include "utils.h"

#include "hmpdf.h"

// defines{{{
#define LOGK // if this macro is defined, integrals over k are evaluated on log grid
#define LOGELL

#define PKINTERP_NK 10000
#define PKINTERP_KMIN 1e-10
#define PKINTERP_KMAX 1e2
#define PKINTERP_TYPE interp_cspline

#define PKINTEGR_KMIN PKINTERP_KMIN
#define PKINTEGR_LIMIT 10000
#define PKINTEGR_KEY 6
#define PKINTEGR_EPSABS 0.0
#define PKINTEGR_EPSREL 1e-6

#define MDEF_GLOBAL hmpdf_mdef_m // this should not be changed,
                                 //    since it's assumed in computation
                                 //    of the mass function

#define CINTERP_NC 10000
#define CINTERP_CMIN 0.1
#define CINTERP_CMAX 100.0
#define CINTERP_TYPE interp_cspline

#define PRINTERP_TYPE interp_steffen // ensure monotonicity is preserved

#define CORRINTERP_N 1000
#define CORRINTERP_TYPE interp_cspline

#define PRWINDOW_INTERP_ELLMIN 1e-2
#define PRWINDOW_INTERP_ELLMAX 1e3
#define PRWINDOW_INTERP_NELL   1000
#define PRWINDOW_INTERP_TYPE   interp_cspline

#define PRWINDOW_INTEGR_LIMIT  1000
#define PRWINDOW_INTEGR_EPSABS 0.0
#define PRWINDOW_INTEGR_EPSREL 1e-6
#define PRWINDOW_INTEGR_KEY    6

#define INTEGR_MINSAMPLES 5

#define PRTILDE_INTEGR_NTHETA 513

#define OPINTERP_TYPE interp_linear
#define TPINTERP_TYPE interp2d_bilinear
#define TPINTEGR_N 30 // sqrt of sample points per pixel
                      // -- verified that this gives sub-percent accuracy
#define COVINTEGR_N 100 // this is slow but we want to be sure that we pick out
                        //   possible oscillatory behaviour

#define BATTINTEGR_LIMIT 1000
#define BATTINTEGR_KEY 6
#define BATTINTEGR_EPSABS 1e-1 // in units of the signal grid spacing
#define BATTINTEGR_EPSREL 1e-4

#define TP_PHI_EQ_TOL 1e-10

#define PU_R2C_MODE FFTW_MEASURE
#define PPDF_C2R_MODE FFTW_MEASURE
#define PC_R2C_MODE FFTW_MEASURE

#define SELL_INTERP_TYPE interp_cspline
#define PS_NELL 1000
#define PS_ELLMIN 1e0
#define PS_ELLMAX 1e5
#define CELL_INTERP_TYPE interp_cspline
#define CPHI_INTERP_TYPE interp_cspline
#define PS_COVINTEGR_N 100

#define COV_STATUS_PERIOD 100

#define NOISE_ELLMIN 1e-2
#define NOISE_ELLMAX 1e12
#define NOISE_LIMIT  1000
#define NOISE_EPSABS 1e-1 // in units of the signal grid spacing
#define NOISE_EPSREL 1e-3
#define NOISE_KEY    6
#define NOISE_ZETAINTERP_N 1000
//}}}

struct DEFAULTS {int Ncores[3]; int verbosity; int warn_is_err;
                 char *class_pre;
                 int Npoints_z[3]; double z_min[3]; double z_max[3];
                 int Npoints_M[3]; double M_min[3]; double M_max[3];
                 long Npoints_signal[3]; double min_kappa[3]; double min_tsz[3]; double max_kappa[3]; double max_tsz[3];
                 int Npoints_theta[3]; double rout_scale[3]; hmpdf_mdef_e rout_rdef[3];
                 double pixel_sidelength[3]; double tophat_radius[3]; double gaussian_fwhm[3];
                 hmpdf_ell_filter_f custom_ell_filter; void *custom_ell_filter_params;
                 hmpdf_k_filter_f custom_k_filter; void *custom_k_filter_params;
                 int Nphi[3]; double phimax[3]; int pixelexactmax[3];
                 double phijitter[3]; double phipwr;
                 hmpdf_integr_mode_e zintegr_type[3]; double zintegr_alpha; double zintegr_beta;
                 hmpdf_integr_mode_e Mintegr_type[3]; double Mintegr_alpha; double Mintegr_beta;
                 double *Duffy08_p; double *Tinker10_p; double *Battaglia12_p;
                 hmpdf_noise_pwr_f noise_pwr; void *noise_pwr_params; };

struct DEFAULTS def;

#endif
