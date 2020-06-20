#ifndef POWERSPECTRUM_H
#define POWERSPECTRUM_H

#include "hmpdf.h"

typedef struct//{{{
{
    int Nell;
    int Nell_corr;

    int created_Cell;
    double *ell;
    double *Cell_1h;
    double *Cell_2h;
    double *Cell_tot;

    int created_Cphi;
    double *phi;
    double *Cphi_1h;
    double *Cphi_2h;
    double *Cphi_tot;
}//}}}
powerspectrum_t;

void null_powerspectrum(hmpdf_obj *d);
void reset_powerspectrum(hmpdf_obj *d);
void hmpdf_get_Cell(hmpdf_obj *d, int Nell, double *ell, double *Cell, hmpdf_Cell_mode_e mode);
void hmpdf_get_Cphi(hmpdf_obj *d, int Nphi, double *phi, double *Cphi, hmpdf_Cell_mode_e mode);

#endif
