#ifndef ONEPOINT_H
#define ONEPOINT_H

#include "data.h"

#include "hmpdf.h"

void create_op(all_data *d);
void get_op(all_data *d, int Nbins, double *binedges, double *out, pdf_cl_uncl mode);

#endif
