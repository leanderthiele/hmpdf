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

void null_powerspectrum(all_data *d);
void reset_powerspectrum(all_data *d);
void get_Cell(all_data *d, int Nell, double *ell, double *Cell, Cell_mode mode);
void get_Cphi(all_data *d, int Nphi, double *phi, double *Cphi, Cell_mode mode);

#endif
