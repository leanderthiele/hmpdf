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

#define MDEF_GLOBAL mdef_m // FIXME this needs to be automatic

#define CINTERP_NC 10000
#define CINTERP_CMIN 0.1
#define CINTERP_CMAX 100.0
#define CINTERP_TYPE interp_cspline

#define PRINTERP_TYPE interp_cspline

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

#define PRTILDE_INTEGR_NTHETA 257

#define OPINTERP_TYPE interp_cspline
#define TPINTERP_TYPE interp2d_bilinear
#define TPINTEGR_N 30 // sqrt of sample points per pixel
                      // -- verified that this gives sub-percent accuracy
#define COVINTEGR_N 100 // this is slow but we want to be sure that we pick out
                        //   possible oscillatory behaviour

// TODO go to fixed point integration here
#define BATTINTEGR_EPSABS 0.0
#define BATTINTEGR_EPSREL 1e-4

// TODO think about this
#define TPREG_MAXPHI (2.0/*arcmin*/*M_PI/60.0/180.0)

#define PU_R2C_MODE FFTW_MEASURE
#define PPDF_C2R_MODE FFTW_MEASURE
#define PC_R2C_MODE FFTW_MEASURE

#define SELL_INTERP_TYPE interp_cspline
#define PS_NELL 1000
#define PS_ELLMIN 1e0
#define PS_ELLMAX 1e5
#define CELL_INTERP_TYPE interp_cspline
#define CPHI_INTERP_TYPE interp_cspline

#define COV_STATUS_PERIOD 100
//}}}

struct DEFAULTS {int Ncores; char *class_pre;
                 int Npoints_z; double z_min; double z_max; double z_source;
                 int Npoints_M; double M_min; double M_max;
                 int Npoints_signal; double signal_min; double max_kappa; double max_tsz;
                 int Npoints_theta; double rout_scale; mdef rout_rdef;
                 double pixel_sidelength; double tophat_radius; double gaussian_fwhm;
                 ell_filter custom_ell_filter; void *custom_ell_filter_params;
                 k_filter custom_k_filter; void *custom_k_filter_params;
                 int Nphi; double phimax; int pixelexactmax;
                 double phijitter; double phipwr; int regularize_tp;
                 int monotonize;
                 integr_mode zintegr_type; double zintegr_alpha; double zintegr_beta;
                 integr_mode Mintegr_type; double Mintegr_alpha; double Mintegr_beta;
                 double *Duffy08_p; double *Tinker10_p; double *Battaglia12_p; };

struct DEFAULTS def;

#endif
