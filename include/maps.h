#ifndef MAPS_H
#define MAPS_H

#include <complex.h>

#include <fftw3.h>

#include <gsl/gsl_rng.h>

#include "hmpdf.h"

typedef struct//{{{
{
    int for_fft;
    double *map; // the complete map in this work space
                 //    includes all objects handled by a specific thread
    double complex *map_comp;
    fftw_plan *p_r2c;
    fftw_plan *p_c2r;

    long bufside; // sidelength of this specific buffer
    double *pos;  // angular separation from object center
                  //     (same shape as buf)
    double *buf;  // buffer for a single object

    gsl_rng *rng;
}//}}}
map_ws;

typedef struct//{{{
{
    double area; // in physical units (rad^2)

    int pxlgrid;

    int created_mem;
    int need_ft;

    int created_ellgrid;
    double *ellgrid;

    int created_sidelengths;
    long Nside;
    long buflen;

    int created_map;
    double *map_real;
    double complex *map_comp;
    fftw_plan *p_r2c;
    fftw_plan *p_c2r;
    double mean;

    int Nws;
    int created_map_ws;
    map_ws **ws;
}//}}}
maps_t;

int null_maps(hmpdf_obj *d);
int reset_maps(hmpdf_obj *d);
int hmpdf_get_map_op(hmpdf_obj *d, int Nbins, double binedges[Nbins+1], double op[Nbins], int new_map);

#endif
