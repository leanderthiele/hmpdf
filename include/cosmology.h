#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include "hmpdf.h"

typedef struct//{{{
{
    int inited_cosmo;

    double *hubble;
    double *comoving;
    double *angular_diameter;
    double *invScrit;
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
}//}}}
cosmology_t;

int null_cosmology(hmpdf_obj *d);
int reset_cosmology(hmpdf_obj *d);
int init_cosmology(hmpdf_obj *d);

#endif
