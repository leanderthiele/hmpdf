#ifndef NUMERICS_H
#define NUMERICS_H

#include "data.h"

double integr_real(int N, double dx, int stride, double *f);
complex integr_comp(int N, double dx, int stride, complex *f);
void init_numerics(all_data *d);

#endif
