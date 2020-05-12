#ifndef FILTER_H
#define FILTER_H

#include "data.h"

#include "hmpdf.h"

void apply_filters(all_data *d, int N, double *ell,
                   double *in, double *out,
                   filter_mode mode, int *z_index);
void init_filters(all_data *d);

#endif
