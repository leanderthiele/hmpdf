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

void null_numerics(hmpdf_obj *d);
void reset_numerics(hmpdf_obj *d);
double integr_real(int N, double dx, int stride, double *f);
complex integr_comp(int N, double dx, int stride, complex *f);
void gauss_fixed_point(hmpdf_integr_mode_e m, int N,
                       double a, double b, double alpha, double beta,
                       double *nodes, double *weights,
                       int neutralize_weights);
void init_numerics(hmpdf_obj *d);

#endif
