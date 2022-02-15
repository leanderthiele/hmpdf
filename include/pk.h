#ifndef PK_H
#define PK_H

#include <gsl/gsl_dht.h>

#include "hmpdf.h"

typedef struct
{
    int inited_pk;
    
    int dht_Nk;
    double *dht_kgrid;
    double *dht_rgrid; // convention is that rout==1

    gsl_dht *dht_ws;
} pk_t;

int null_pk(hmpdf_obj *d);
int reset_pk(hmpdf_obj *d);
int init_pk(hmpdf_obj *d);
int hmpdf_get_pk(hmpdf_obj *d, double z, int Npoints, double k[Npoints],
                 double Pk_1h[Npoints], double Pk_2h[Npoints]);

#endif
