#ifndef POWERSPECTRUM_H
#define POWERSPECTRUM_H

#include "hmpdf.h"

void get_Cell(all_data *d, int Nell, double *ell, double *Cell, Cell_mode mode);
void get_Cphi(all_data *d, int Nphi, double *phi, double *Cphi, Cell_mode mode);

#endif
