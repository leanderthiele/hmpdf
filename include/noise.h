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

void null_noise(hmpdf_obj *d);
void reset_noise(hmpdf_obj *d);
void init_noise(hmpdf_obj *d);

void noise_vect(hmpdf_obj *d, double *in, double *out);
void noise_matr(hmpdf_obj *d, double *in, double *out);

#endif
