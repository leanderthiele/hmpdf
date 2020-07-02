#ifndef NUMERICS_H
#define NUMERICS_H

#include "hmpdf.h"

typedef struct//{{{
{
    int inited_numerics;

    int monotonize;

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
    double *user_signalgrid;
    double *lambdagrid;

    int Nsignal_noisy;
    double *signalgrid_noisy;

    // TODO move these into covariance
    int Nphi;
    int pixelexactmax;
    double phimax;
    double phijitter;
    double phipwr;
    double *phigrid;
    double *phiweights;

    // following options have no effect if
    // monotonize is set to false = 0
    hmpdf_integr_mode_e zintegr_type;
    double zintegr_alpha;
    double zintegr_beta;

    hmpdf_integr_mode_e Mintegr_type;
    double Mintegr_alpha;
    double Mintegr_beta;
    //

    double zsource;
}//}}}
numerics_t;

int null_numerics(hmpdf_obj *d);
int reset_numerics(hmpdf_obj *d);
double integr_real(int N, double dx, int stride, double *f);
complex integr_comp(int N, double dx, int stride, complex *f);
int init_numerics(hmpdf_obj *d);

#endif
