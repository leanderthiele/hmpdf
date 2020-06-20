#ifndef TWOPOINT_H
#define TWOPOINT_H

#include <fftw3.h>

#include "hmpdf.h"

typedef struct//{{{
{
    // holds the unclustered term, bc is added in the end
    double *pdf_real; // [ Nsignal * Nsignal+2 ]
    complex *pdf_comp; // not malloced
    fftw_plan pu_r2c; // pdf_real -> pdf_comp
    fftw_plan ppdf_c2r; // pdf_comp -> pdf_real

    // holds the clustering term
    complex *bc; // [ Nsignal * Nsignal/2+1

    // holds the z-specific clustering contribution
    double *tempc_real; // [ Nsignal * Nsignal+2 ]
    complex *tempc_comp; // not malloced
    fftw_plan pc_r2c; // tempc_real -> tempc_comp
}//}}}
twopoint_workspace;

typedef struct//{{{
{
    // phi-independent quantities, to compute only once
    int created_phi_indep;
    double ***dtsq; // [ z_index, M_index, signal_index ]
    double ***t; // [ z_index, M_index, signal_index ]
    complex **ac; // [ z_index, lambda_index ]
    complex *au; // [ lambda_index ] // allocated with fftw_malloc
    
    double last_phi;

    // buffer regions --> one for each core
    twopoint_workspace *ws;

    double *pdf;
    double *pdf_noisy;
}//}}}
twopoint_t;

twopoint_workspace *new_tp_ws(int N);

void null_twopoint(hmpdf_obj *d);
void reset_twopoint(hmpdf_obj *d);
void create_phi_indep(hmpdf_obj *d);
void create_tp(hmpdf_obj *d, double phi, twopoint_workspace *ws);
void hmpdf_get_tp(hmpdf_obj *d, double phi, int Nbins, double *binedges, double *out, int noisy);

#endif
