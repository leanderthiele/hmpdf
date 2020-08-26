#ifndef NOISE_H
#define NOISE_H

#include <stdlib.h>

#include "hmpdf.h"

typedef struct
{//{{{
    double noise;

    size_t len_kernel;
    double *toepl;
}//}}}
noise_t;

int null_noise(hmpdf_obj *d);
int reset_noise(hmpdf_obj *d);
int init_noise(hmpdf_obj *d);

int noise_vect(hmpdf_obj *d, double *in, double *out);
int noise_matr(hmpdf_obj *d, double *in, double *out);

#endif
