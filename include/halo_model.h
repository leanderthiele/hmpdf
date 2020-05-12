#ifndef HALO_MODEL_H
#define HALO_MODEL_H

#include "utils.h"
#include "data.h"

#include "hmpdf.h"

double NFW_fundamental(all_data *d, int z_index, int M_index, double *rs);
double Mconv(all_data *d, int z_index, int M_index, mdef mdef_out, double *R, double *c);
void init_halo_model(all_data *d);

#endif
