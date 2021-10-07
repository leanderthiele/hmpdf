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
#include "cosmology.h"
#include "profiles.h"
#include "halo_model.h"
#include "numerics.h"
#include "mean.h"

#include "hmpdf.h"



static int
Arnaud05_T_of_M(hmpdf_obj *d, int z_index, int M_index, double *out)
{//{{{
    STARTFCT

    double M200c, R200c, c200c;
    SAFEHMPDF(Mconv(d, z_index, M_index, hmpdf_mdef_c, &M200c, &R200c, &c200c));

    static const double A = 5.34e14; // Msun
    static const double alpha = 1.72;
    static const double m_bias = 0.2;

    // definition of Arnaud+05: h(z) = H(z) / H0
    double hz = d->c->hubble[z_index] * SPEEDOFLIGHT / d->c->h;

    *out = 5.0 * pow( hz * (1.0-m_bias) * M200c / A, 1.0/alpha );

    ENDFCT
}//}}}


int
hmpdf_get_mean(hmpdf_obj *d, double *mean)
{//{{{
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
}//}}}


int
hmpdf_get_mean_T(hmpdf_obj *d, double *mean)
{
    STARTFCT

    CHECKINIT;

    HMPDFCHECK(d->p->stype != hmpdf_tsz,
               "hmpdf_get_mean_T only makes sense for tSZ signal.");

    double mean_yT = 0.0;
    double mean_y = 0.0;

    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        for (int M_index=0; M_index<d->n->NM; M_index++)
        {
            double temp_y;
            double temp_T;
            SAFEHMPDF(integrate_profile(d, z_index, M_index, &temp_y));
            SAFEHMPDF(Arnaud05_T_of_M(d, z_index, M_index, &temp_T));

            double weight = d->h->hmf[z_index][M_index]
                            * gsl_pow_2(d->c->comoving[z_index])
                            / d->c->hubble[z_index]
                            * d->n->Mweights[M_index] * d->n->zweights[z_index];

            mean_y += temp_y * weight;
            mean_yT += temp_y * temp_T * weight;
        }
    }

    *mean = mean_yT / mean_y;

    ENDFCT
}
