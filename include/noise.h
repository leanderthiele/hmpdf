#ifndef NOISE_H
#define NOISE_H

#include "hmpdf.h"

typedef struct
{//{{{
    double noise;

    int len_kernel;
    double *toepl;
}//}}}
noise_t;

void null_noise(all_data *d);
void reset_noise(all_data *d);
void init_noise(all_data *d);

void noise_vect(all_data *d, double *in, double *out);
void noise_matr(all_data *d, double *in, double *out);

#endif
