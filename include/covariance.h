#ifndef COVARIANCE_H
#define COVARIANCE_H

#include "twopoint.h"
#include "hmpdf.h"

typedef struct//{{{
{
    double *Cov;
    double *corr_diagn;

    int Nws;
    int created_tp_ws;
    twopoint_workspace **ws;

    int created_phigrid;
}//}}}
covariance_t;

void null_covariance(all_data *d);
void reset_covariance(all_data *d);
void get_cov(all_data *d, int Nbins, double *binedges, double *out, int noisy, char *fname);

#endif
