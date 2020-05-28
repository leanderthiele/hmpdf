#ifndef ONEPOINT_H
#define ONEPOINT_H

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

void null_onepoint(all_data *d);
void reset_onepoint(all_data *d);
void create_op(all_data *d);
void get_op(all_data *d, int Nbins, double *binedges, double *out, pdf_cl_uncl mode, int noisy);

#endif
