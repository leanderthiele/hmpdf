#ifndef PROFILES_H
#define PROFILES_H

#include "data.h"

#include "hmpdf.h"

void init_profiles(all_data *d);
void create_conj_profiles(all_data *d);
void create_filtered_profiles(all_data *d);
void create_breakpoints_or_monotonize(all_data *d);

void s_of_t(all_data *d, int z_index, int M_index, int Nt, double *t, double *s);
void s_of_ell(all_data *d, int z_index, int M_index, int Nell, double *ell, double *s);
void dtsq_of_s(all_data *d, int z_index, int M_index, double *dtsq);
void t_of_s(all_data *d, int z_index, int M_index, double *t);

#endif
