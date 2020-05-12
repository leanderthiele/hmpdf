#ifndef FILTER_H
#define FILTER_H

#include "data.h"

#include "hmpdf.h"

typedef enum//{{{
{
    filter_pdf,
    filter_ps,
    filter_end,
}//}}}
filter_mode;

void apply_filters(all_data *d, int N, double *ell,
                   double *in, double *out, filter_mode mode);
void init_filters(all_data *d);
#endif
