#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <class.h>

#include <gsl/gsl_math.h>

#include "utils.h"
#include "data.h"
#include "class_interface.h"
#include "cosmology.h"

void null_cosmology(all_data *d)
{//{{{
    d->c->inited_cosmo = 0;
    d->c->hubble = NULL;
    d->c->comoving = NULL;
    d->c->angular_diameter = NULL;
    d->c->Scrit = NULL;
    d->c->Dsq = NULL;
    d->c->rho_m = NULL;
    d->c->rho_c = NULL;
    d->c->Om = NULL;
}//}}}

void reset_cosmology(all_data *d)
{//{{{
    if (d->c->hubble != NULL) { free(d->c->hubble); }
    if (d->c->comoving != NULL) { free(d->c->comoving); }
    if (d->c->angular_diameter != NULL) { free(d->c->angular_diameter); }
    if (d->c->Scrit != NULL) { free(d->c->Scrit); }
    if (d->c->Dsq != NULL) { free(d->c->Dsq); }
    if (d->c->rho_m != NULL) { free(d->c->rho_m); }
    if (d->c->rho_c != NULL) { free(d->c->rho_c); }
    if (d->c->Om != NULL) { free(d->c->Om); }
}//}}}

static
void alloc_cosmo(all_data *d)
{//{{{
    printf("\talloc_cosmo\n");
    d->c->hubble = (double *)malloc(d->n->Nz * sizeof(double));
    d->c->comoving = (double *)malloc(d->n->Nz * sizeof(double));
    d->c->angular_diameter = (double *)malloc(d->n->Nz * sizeof(double));
    d->c->Dsq = (double *)malloc(d->n->Nz * sizeof(double));
    d->c->Om = (double *)malloc(d->n->Nz * sizeof(double));
    d->c->rho_m = (double *)malloc(d->n->Nz * sizeof(double));
    d->c->rho_c = (double *)malloc(d->n->Nz * sizeof(double));
    if (d->p->stype == kappa)
    {
        d->c->Scrit = (double *)malloc(d->n->Nz * sizeof(double));
    }
}//}}}

static
void fill_background(all_data *d)
{//{{{
    printf("\tfill_background\n");
    double tau; // conformal time in Mpc
    int index = 0; // some internal CLASS thing

    struct background *ba = (struct background *)d->cls->ba;

    double *pvecback = (double *)malloc(ba->bg_size * sizeof(double));
    // get z=0 numbers
    d->c->h = ba->h;
    d->c->rho_c_0 = 3.0*gsl_pow_2(SPEEDOFLIGHT)/8.0/M_PI/GNEWTON
                    * gsl_pow_2(ba->H0);
    d->c->Om_0 = ba->Omega0_m;
    d->c->rho_m_0 = d->c->Om_0 * d->c->rho_c_0;
    d->c->Ob_0 = ba->Omega0_b;

    // get background
    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        // get conformal time at given redshift
        background_tau_of_z(ba, d->n->zgrid[z_index], &tau);
        // write interpolated background quantities into pvecback
        background_at_tau(ba, tau, ba->long_info,
                          ba->inter_normal, &index, pvecback);

        d->c->hubble[z_index] = pvecback[ba->index_bg_H]; // 1/Mpc
        d->c->comoving[z_index] = pvecback[ba->index_bg_conf_distance]; // Mpc
        d->c->angular_diameter[z_index] = pvecback[ba->index_bg_ang_distance]; // Mpc, physical
        d->c->Dsq[z_index] = gsl_pow_2(pvecback[ba->index_bg_D]); // squared growth factor
        d->c->Om[z_index] = pvecback[ba->index_bg_Omega_m];
        d->c->rho_c[z_index] = 3.0*gsl_pow_2(SPEEDOFLIGHT)/8.0/M_PI/GNEWTON
                               * gsl_pow_2(d->c->hubble[z_index]);
        d->c->rho_m[z_index] = d->c->Om[z_index] * d->c->rho_c[z_index];
    }
    if (d->p->stype == kappa) // need to compute critical surface density
    {
        // find distances to source position
        background_tau_of_z(ba, d->n->zsource, &tau);
        background_at_tau(ba, tau, ba->long_info,
                          ba->inter_normal, &index, pvecback);
        double chi_s = pvecback[ba->index_bg_conf_distance];
        double dA_s = pvecback[ba->index_bg_ang_distance];
        // fill the Scrit grid
        for (int z_index=0; z_index<d->n->Nz; z_index++)
        {
            d->c->Scrit[z_index] = gsl_pow_2(SPEEDOFLIGHT)/4.0/M_PI/GNEWTON
                                   * dA_s / d->c->angular_diameter[z_index]
                                   * (1.0+d->n->zsource)
                                   / (chi_s - d->c->comoving[z_index]);
        }
    }

    free(pvecback);
}//}}}

void init_cosmology(all_data *d)
{//{{{
    printf("In cosmology.h -> init_cosmo.\n");
    alloc_cosmo(d);
    fill_background(d);
}//}}}

