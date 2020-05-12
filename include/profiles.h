#ifndef PROFILES_H
#define PROFILES_H

#include "data.h"

#include "hmpdf.h"

int breakpoints(all_data *d, int z_index);
// TODO remove this function
void create_conj_profiles(all_data *d);
void s_of_t(all_data *d, int z_index, int M_index, int Nt, double *t, double *s);
void s_of_l(all_data *d, int z_index, int M_index, int Nl, double *l, double *s);
void dtsq_of_s(all_data *d, int z_index, int M_index, double *dtsq);
void t_of_s(all_data *d, int z_index, int M_index, double *t);
void init_profiles(all_data *d);

#endif
