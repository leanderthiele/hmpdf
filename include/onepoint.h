#ifndef ONEPOINT_H
#define ONEPOINT_H

#include <complex.h>

#include "hmpdf.h"

typedef struct
{//{{{
    int created_op;
    double *PDFu;
    double *PDFc;

    int created_noisy_op;
    double *PDFu_noisy;
    double *PDFc_noisy;

    double signalmeanu;
    double signalmeanc;
}//}}}
onepoint_t;

int null_onepoint(hmpdf_obj *d);
int reset_onepoint(hmpdf_obj *d);
int correct_phase1d(hmpdf_obj *d, double complex *x, long stride, int sgn);
int create_op(hmpdf_obj *d);
int create_noisy_op(hmpdf_obj *d);
int pdf_adjust_binedges(hmpdf_obj *d, int Nbins, double binedges_in[Nbins+1], double binedges_out[Nbins+1], double mean);
int pdf_check_user_input(hmpdf_obj *d, int Nbins, double binedges[Nbins+1], int noisy);
int hmpdf_get_op(hmpdf_obj *d, int Nbins, double binedges[Nbins+1], double op[Nbins], int incl_2h, int noisy);

#endif
