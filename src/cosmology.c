#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <class.h>

#include <gsl/gsl_math.h>

#include "utils.h"
#include "data.h"
#include "cosmology.h"

struct class_interface_s
{//{{{
    char *class_ini;
    char *class_pre;

    struct precision *pr;
    struct background *ba;
    struct thermo *th;
    struct primordial *pm;
    struct perturbs *pt;
    struct nonlinear *nl;
    struct transfers *tr;
    struct spectra *sp;
    struct lensing *le;
    struct output *op;

    ErrorMsg errmsg;
};//}}}

typedef struct class_interface_s cls;

static
void alloc_cosmo(all_data *d)
{//{{{
    printf("\talloc_cosmo\n");
    d->c->hubble = (double *)malloc(d->n->gr->Nz * sizeof(double));
    d->c->comoving = (double *)malloc(d->n->gr->Nz * sizeof(double));
    d->c->angular_diameter = (double *)malloc(d->n->gr->Nz * sizeof(double));
    d->c->Dsq = (double *)malloc(d->n->gr->Nz * sizeof(double));
    d->c->Om = (double *)malloc(d->n->gr->Nz * sizeof(double));
    d->c->rho_m = (double *)malloc(d->n->gr->Nz * sizeof(double));
    d->c->rho_c = (double *)malloc(d->n->gr->Nz * sizeof(double));
    if (d->p->stype == kappa)
    {
        d->c->Scrit = (double *)malloc(d->n->gr->Nz * sizeof(double));
    }
}//}}}

static
void fill_background(all_data *d)
{//{{{
    printf("\tfill_background\n");
    double tau; // conformal time in Mpc
    int index = 0; // some internal CLASS thing
    cls *_c = (cls *)d->cls;
    double *pvecback = (double *)malloc(_c->ba->bg_size * sizeof(double));
    // get z=0 numbers
    d->c->h = _c->ba->h;
    d->c->rho_c_0 = 3.0*gsl_pow_2(SPEEDOFLIGHT)/8.0/M_PI/GNEWTON
                    * gsl_pow_2(_c->ba->H0);
    d->c->Om_0 = _c->ba->Omega0_m;
    d->c->rho_m_0 = d->c->Om_0 * d->c->rho_c_0;
    d->c->Ob_0 = _c->ba->Omega0_b;

    // get background
    for (int z_index=0; z_index<d->n->gr->Nz; z_index++)
    {
        // get conformal time at given redshift
        background_tau_of_z(_c->ba, d->n->gr->zgrid[z_index], &tau);
        // write interpolated background quantities into pvecback
        background_at_tau(_c->ba, tau, _c->ba->long_info,
                          _c->ba->inter_normal, &index, pvecback);

        d->c->hubble[z_index] = pvecback[_c->ba->index_bg_H]; // 1/Mpc
        d->c->comoving[z_index] = pvecback[_c->ba->index_bg_conf_distance]; // Mpc
        d->c->angular_diameter[z_index] = pvecback[_c->ba->index_bg_ang_distance]; // Mpc
        d->c->Dsq[z_index] = gsl_pow_2(pvecback[_c->ba->index_bg_D]); // squared growth factor
        d->c->Om[z_index] = pvecback[_c->ba->index_bg_Omega_m];
        d->c->rho_c[z_index] = 3.0*gsl_pow_2(SPEEDOFLIGHT)/8.0/M_PI/GNEWTON
                               * gsl_pow_2(d->c->hubble[z_index]);
        d->c->rho_m[z_index] = d->c->Om[z_index] * d->c->rho_c[z_index];
    }
    if (d->p->stype == kappa) // need to compute critical surface density
    {
        // find distances to source position
        background_tau_of_z(_c->ba, d->n->zsource, &tau);
        background_at_tau(_c->ba, tau, _c->ba->long_info,
                          _c->ba->inter_normal, &index, pvecback);
        double chi_s = pvecback[_c->ba->index_bg_conf_distance];
        double dA_s = pvecback[_c->ba->index_bg_ang_distance];
        // fill the Scrit grid
        for (int z_index=0; z_index<d->n->gr->Nz; z_index++)
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

