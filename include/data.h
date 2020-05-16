#ifndef DATA_H
#define DATA_H

#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_dht.h>

#include "hmpdf.h"

#include "utils.h"

struct grids
{//{{{
    int Nz;
    double zmin;
    double zmax;
    double *zgrid;
    double *zweights;

    int NM;
    double Mmin;
    double Mmax;
    double *Mgrid;
    double *Mweights;

    int Nsignal;
    double signalmin;
    double signalmax;
    double *signalgrid;
    double *lambdagrid;

    // TODO move these into covariance
    int Nphi;
    int pixelexactmax;
    double phimax;
    double phijitter;
    double phipwr;
    double *phigrid;
    double *phiweights;
};//}}}

struct numerics
{//{{{
    int inited_numerics;

    int monotonize;

    // following options have no effect if
    // monotonize is set to false = 0
    integr_mode zintegr_type;
    double zintegr_alpha;
    double zintegr_beta;

    integr_mode Mintegr_type;
    double Mintegr_alpha;
    double Mintegr_beta;
    //

    struct grids *gr;

    double zsource;
};//}}}

struct cosmo
{//{{{
    int inited_cosmo;

    double *hubble;
    double *comoving;
    double *angular_diameter;
    double *Scrit;
    double *Dsq;
    // simple quantities
    double h;
    double rho_c_0;
    double rho_m_0;
    double *rho_m;
    double *rho_c;
    double *Om;
    double Om_0;
    double Ob_0;
};//}}}

struct power
{//{{{
    int inited_power;

    double *k_arr;
    double *Pk_arr;
    interp1d *Pk_interp;

    double **ssq;
    double autocorr;

    int created_corr;
    gsl_spline *corr_interp;
    int Ncorr_accel;
    gsl_interp_accel **corr_accel;
};//}}}

struct halo
{//{{{
    int inited_halo;

    double *Duffy08_params;
    double *Tinker10_params;

    double **hmf;
    double **bias;

    gsl_spline *c_interp;
    gsl_interp_accel *c_accel;
};//}}}

struct profiles
{//{{{
    int inited_profiles;

    double *Battaglia12_params;

    signaltype stype;

    double rout_scale;
    int rout_def;
    int Ntheta;
    double *decr_tgrid;
    double *incr_tgrid;
    double *decr_tsqgrid;
    double *reci_tgrid; // reciprocal space grid

    gsl_interp_accel *incr_tgrid_accel;
    gsl_interp_accel *reci_tgrid_accel;

    double ***profiles; // each profile has as zero entry theta out and then the profile

    int created_conj_profiles;
    double ***conj_profiles; // each profile has as zero entry the rescaling such that reci_thetagrid -> ell

    int created_filtered_profiles;

    int prtilde_Ntheta;
    double *prtilde_thetagrid;

    int created_breakpoints;
    int *breakpoints;

    gsl_dht *dht_ws;
};//}}}

struct filters
{//{{{
    int inited_filters;

    int Nfilters;
    int *z_dependent;
    filter_fct *ffilters;

    double pixelside;
    gsl_spline **quadraticpixel_interp;
    gsl_interp_accel **quadraticpixel_accel;
    double *quadraticpixel_ellmin;
    double *quadraticpixel_ellmax;

    double tophat_radius;
    double gaussian_sigma;

    ell_filter custom_ell;
    void *custom_ell_p;
    k_filter custom_k;
    void *custom_k_p;
};//}}}

struct onepoint
{//{{{
    int created_op;
    double *PDFu;
    double *PDFc;

    double signalmeanu;
    double signalmeanc;
};//}}}

struct twopoint_workspace_s
{//{{{
    // holds the unclustered term, bc is added in the end
    double *pdf_real; // [ Nsignal * Nsignal+2 ]
    complex *pdf_comp; // not malloced
    fftw_plan pu_r2c; // pdf_real -> pdf_comp
    fftw_plan ppdf_c2r; // pdf_comp -> pdf_real

    // holds the clustering term
    complex *bc; // [ Nsignal * Nsignal/2+1

    // holds the z-specific clustering contribution
    double *tempc_real; // [ Nsignal * Nsignal+2 ]
    complex *tempc_comp; // not malloced
    fftw_plan pc_r2c; // tempc_real -> tempc_comp
};//}}}

typedef struct twopoint_workspace_s twopoint_workspace;

struct twopoint
{//{{{
    // phi-independent quantities, to compute only once
    int created_phi_indep;
    double ***dtsq; // [ z_index, M_index, signal_index ]
    double ***t; // [ z_index, M_index, signal_index ]
    complex **ac; // [ z_index, lambda_index ]
    complex *au; // [ lambda_index ] // allocated with fftw_malloc
    
    int regularize;

    // buffer regions --> one for each core
    twopoint_workspace *ws;
};//}}}

struct powerspectrum
{//{{{
    int Nell;
    int Nell_corr;

    int created_Cell;
    double *ell;
    double *Cell_1h;
    double *Cell_2h;
    double *Cell_tot;

    int created_Cphi;
    double *phi;
    double *Cphi_1h;
    double *Cphi_2h;
    double *Cphi_tot;
};//}}}

struct covariance
{//{{{
    double *Cov;
    double *corr_diagn;

    int Nws;
    twopoint_workspace **ws;
};//}}}

struct all_data_s
{//{{{
    int Ncores;

    char *class_ini;
    char *class_pre;

    struct numerics *n;
//    struct class_interface *cls;
    void *cls;
    struct cosmo *c;
    struct power *pwr;
    struct halo *h;
    struct filters *f;
    struct profiles *p;
    struct onepoint *op;
    struct twopoint *tp;
    struct powerspectrum *ps;
    struct covariance *cov;
};//}}}

all_data *new_data(void);
void reset_data(all_data *d);
void delete_data(all_data *d);

#endif
