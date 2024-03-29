#ifndef NUMERICS_H
#define NUMERICS_H

#include <complex.h>

#include "hmpdf.h"

typedef struct//{{{
{
    int inited_numerics;

    int Nz;
    double zmin;
    double zmax;
    double *zgrid;
    double *zweights;

    hmpdf_dndz_f dndz;
    void *dndz_params;

    int NM;
    double Mmin;
    double Mmax;
    double *Mgrid;
    double *Mweights;

    hmpdf_mass_cuts_f mass_cuts;
    void *mass_cuts_params;

    long Nsignal;
    long Nsignal_negative;
    double signalmin;
    double signalmax;
    double *signalgrid;
    double *lambdagrid;

    long Nsignal_noisy;
    double *signalgrid_noisy;
    double *lambdagrid_noisy;

    int Nphi;
    int pixelexactmax;
    double phimax;
    double phijitter;
    double phipwr;
    double *phigrid;
    double *phiweights;

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
double complex integr_comp(int N, double dx, int stride, double complex *f);
int init_numerics(hmpdf_obj *d);

#endif
