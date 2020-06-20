#ifndef COVARIANCE_H
#define COVARIANCE_H

#include "twopoint.h"
#include "hmpdf.h"

typedef struct//{{{
{
    double *Cov;
    double *Cov_noisy;
    double *corr_diagn;

    int Nws;
    int created_tp_ws;
    twopoint_workspace **ws;

    int created_phigrid;

    int created_cov;
    int created_noisy_cov;
}//}}}
covariance_t;

void null_covariance(hmpdf_obj *d);
void reset_covariance(hmpdf_obj *d);
void hmpdf_get_cov(hmpdf_obj *d, int Nbins, double *binedges, double *out, int noisy);
void hmpdf_get_cov_diagnostics(hmpdf_obj *d, int *Nphi, double **phi,
                               double **phiweights, double **corr_diagn);

#endif
