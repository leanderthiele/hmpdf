#ifndef HMPDF_H
#define HMPDF_H

typedef enum//{{{
{
    mdef_c,
    mdef_v,
    mdef_m,
}//}}}
mdef;

typedef enum//{{{
{
    kappa,
    tsz,
    end_signaltype,
}//}}}
signaltype;

typedef enum//{{{
{
    legendre,
    chebyshev,
    gegenbauer,
    jacobi,
    laguerre,
    hermite,
    exponential,
    rational,
    chebyshev2,
}//}}}
integr_mode;

typedef enum//{{{
{
    hmpdf_N_cores,
    hmpdf_class_pre,
    hmpdf_N_z,
    hmpdf_z_min,
    hmpdf_z_max,
    hmpdf_z_source,
    hmpdf_N_M,
    hmpdf_M_min,
    hmpdf_M_max,
    hmpdf_N_signal,
    hmpdf_signal_min,
    hmpdf_signal_max,
    hmpdf_N_theta,
    hmpdf_rout_scale,
    hmpdf_rout_rdef,
    hmpdf_pixel_side,
    hmpdf_tophat_radius,
    hmpdf_gaussian_fwhm,
    hmpdf_custom_ell_filter,
    hmpdf_custom_ell_filter_params,
    hmpdf_custom_k_filter,
    hmpdf_custom_k_filter_params,
    hmpdf_N_phi,
    hmpdf_phi_max,
    hmpdf_pixelexact_max,
    hmpdf_phi_jitter,
    hmpdf_phi_pwr,
    hmpdf_regularize_tp,
    hmpdf_monotonize,
    hmpdf_zintegr_type,
    hmpdf_zintegr_alpha,
    hmpdf_zintegr_beta,
    hmpdf_Mintegr_type,
    hmpdf_Mintegr_alpha,
    hmpdf_Mintegr_beta,
    hmpdf_Duffy08_conc_params,
    hmpdf_Tinker10_hmf_params,
    hmpdf_Battaglia12_tsz_params,
    hmpdf_noise,
    hmpdf_end_configs, // keep this last
}//}}}
configs;

// custom filters
typedef double (*ell_filter)(double /*ell 2d wavenumber*/,
                             void * /*user parameters*/);
typedef double (*k_filter)(double /*k comoving 1/Mpc*/,
                           double /*z redshift*/,
                           void * /*user parameters*/);

typedef struct all_data_s all_data;

all_data *new_data(void);

void delete_data(all_data *d);

void init(all_data *d, char *class_ini, signaltype stype, ...);

typedef enum//{{{
{
    uncl,
    cl,
}//}}}
pdf_cl_uncl;

void get_op(all_data *d, int Nbins, double *binedges, double *out, pdf_cl_uncl mode, int noisy);

void get_tp(all_data *d, double phi, int Nbins, double *binedges, double *out);

void get_cov(all_data *d, int Nbins, double *binedges, double *out, int noisy);

void get_cov_diagnostics(all_data *d, int *Nphi, double **phi, double **phiweights, double **corr_diagn);

typedef enum//{{{
{
    onehalo,
    twohalo,
    total,
}//}}}
Cell_mode;

void get_Cell(all_data *d, int Nell, double *ell, double *Cell, Cell_mode mode);

void get_Cphi(all_data *d, int Nphi, double *phi, double *Cphi, Cell_mode mode);

#endif
