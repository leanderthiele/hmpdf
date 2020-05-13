#ifndef TWOPOINT_H
#define TWOPOINT_H

#include "data.h"

#include "hmpdf.h"

twopoint_workspace *new_tp_ws(int N);

void create_phi_indep(all_data *d);
void create_tp(all_data *d, double phi, twopoint_workspace *ws);
void get_tp(all_data *d, double phi, int Nbins, double *binedges, double *out);

#endif
