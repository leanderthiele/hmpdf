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

void null_covariance(all_data *d);
void reset_covariance(all_data *d);
void get_cov(all_data *d, int Nbins, double *binedges, double *out, int noisy);
void get_cov_diagnostics(all_data *d, int *Nphi, double *phi,
                         double *phiweights, double *corr_diagn);

#endif
