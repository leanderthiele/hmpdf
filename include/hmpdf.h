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
    hmpdf_end_configs, // keep this last
}//}}}
configs;

void savetxt(char *fname, int Nlines, int Nvec, ...);
double **loadtxt(char *fname, int *Nlines, int Nvec);

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

void get_op(all_data *d, int Nbins, double *binedges, double *out, pdf_cl_uncl mode);

void get_tp(all_data *d, double phi); // TODO binning

void get_cov(all_data *d, int Nbins, double *binedges, double *out, char *fname);

#endif
