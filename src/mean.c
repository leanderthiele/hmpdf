#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include <fftw3.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "utils.h"
#include "configs.h"
#include "object.h"
#include "profiles.h"
#include "numerics.h"
#include "mean.h"

#include "hmpdf.h"



int
hmpdf_get_mean(hmpdf_obj *d, double *mean)
{
    STARTFCT

    CHECKINIT;

    *mean = 0.0;

    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        for (int M_index=0; M_index<d->n->NM; M_index++)
        {
            double temp;
            SAFEHMPDF(integrate_profile(d, z_index, M_index, &temp));

            *mean += temp * d->h->hmf[z_index][M_index]
                     * gsl_pow_2(d->c->comoving[z_index])
                     / d->c->hubble[z_index]
                     * d->n->Mweights[M_index] * d->n->zweights[z_index];
        }
    }

    ENDFCT
}
