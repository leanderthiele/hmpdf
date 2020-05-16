#ifndef POWER_H
#define POWER_H

#include "data.h"

void create_corr(all_data *d);
double Pk_linear(all_data *d, double k/*if LOGK is defined, this is log(k)*/);
double corr(all_data *d, int z_index, double phi);
void init_power(all_data *d);

#endif
