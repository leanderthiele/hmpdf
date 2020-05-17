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
    integr_mode zintegr_type;
    double zintegr_alpha;
    double zintegr_beta;

    integr_mode Mintegr_type;
    double Mintegr_alpha;
    double Mintegr_beta;
    //

    double zsource;
}//}}}
numerics_t;

void null_numerics(all_data *d);
void reset_numerics(all_data *d);
double integr_real(int N, double dx, int stride, double *f);
complex integr_comp(int N, double dx, int stride, complex *f);
void init_numerics(all_data *d);

#endif
