#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <class.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "utils.h"
#include "object.h"
#include "class_interface.h"
#include "configs.h"
#include "cosmology.h"

int
null_cosmology(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->c->inited_cosmo = 0;
    d->c->hubble = NULL;
    d->c->comoving = NULL;
    d->c->angular_diameter = NULL;
    d->c->invScrit = NULL;
    d->c->Dsq = NULL;
    d->c->rho_m = NULL;
    d->c->rho_c = NULL;
    d->c->Om = NULL;

    ENDFCT
}//}}}

int
reset_cosmology(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\treset_cosmology\n");

    if (d->c->hubble != NULL) { free(d->c->hubble); }
    if (d->c->comoving != NULL) { free(d->c->comoving); }
    if (d->c->angular_diameter != NULL) { free(d->c->angular_diameter); }
    if (d->c->invScrit != NULL) { free(d->c->invScrit); }
    if (d->c->Dsq != NULL) { free(d->c->Dsq); }
    if (d->c->rho_m != NULL) { free(d->c->rho_m); }
    if (d->c->rho_c != NULL) { free(d->c->rho_c); }
    if (d->c->Om != NULL) { free(d->c->Om); }

    ENDFCT
}//}}}

static int
alloc_cosmo(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\talloc_cosmo\n");

    SAFEALLOC(d->c->hubble, malloc(d->n->Nz * sizeof(double)));
    SAFEALLOC(d->c->comoving, malloc(d->n->Nz * sizeof(double)));
    SAFEALLOC(d->c->angular_diameter, malloc(d->n->Nz * sizeof(double)));
    SAFEALLOC(d->c->Dsq, malloc(d->n->Nz * sizeof(double)));
    SAFEALLOC(d->c->Om, malloc(d->n->Nz * sizeof(double)));
    SAFEALLOC(d->c->rho_m, malloc(d->n->Nz * sizeof(double)));
    SAFEALLOC(d->c->rho_c, malloc(d->n->Nz * sizeof(double)));

    if (d->p->stype == hmpdf_kappa)
    {
        SAFEALLOC(d->c->invScrit, malloc(d->n->Nz * sizeof(double)));
    }

    ENDFCT
}//}}}

typedef struct
{
    // CLASS stuff
    double *pvecback;
    struct background *ba;
    int index;

    hmpdf_dndz_f dndz;
    void *dndz_params;

    double chi_z;

    int status;
} dndz_integr_params;

static int
dndz_integr_kernel(double z, dndz_integr_params *p, double *out)
{
    STARTFCT

    double tau;
    SAFECLASS(background_tau_of_z(p->ba, z, &tau),
              p->ba->error_message);
    SAFECLASS(background_at_tau(p->ba, tau, long_info, inter_normal, &p->index, p->pvecback),
              p->ba->error_message);

    double chi_this_z = p->pvecback[p->ba->index_bg_conf_distance];

    *out = (1.0 - chi_this_z/p->chi_z) * p->dndz(z, p->dndz_params);

    ENDFCT
}

static double
dndz_integr_f(double z, void *params)
{
    dndz_integr_params *p = (dndz_integr_params *)params;

    double out;
    p->status = dndz_integr_kernel(z, p, &out);
    return out;
}

static int
fill_background(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\tfill_background\n");

    double tau; // conformal time in Mpc
    int index = 0; // some internal CLASS thing

    struct background *ba = (struct background *)d->cls->ba;

    double *pvecback;
    SAFEALLOC(pvecback, malloc(ba->bg_size * sizeof(double)));
    // get z=0 numbers
    d->c->h = ba->h;
    d->c->rho_c_0 = 3.0 * gsl_pow_2(SPEEDOFLIGHT) / 8.0 / M_PI / GNEWTON
                    * gsl_pow_2(ba->H0);
    d->c->Om_0 = ba->Omega0_m;
    d->c->rho_m_0 = d->c->Om_0 * d->c->rho_c_0;
    d->c->Ob_0 = ba->Omega0_b;

    // get background
    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        // get conformal time at given redshift
        SAFECLASS(background_tau_of_z(ba, d->n->zgrid[z_index], &tau),
                  ba->error_message);
        // write interpolated background quantities into pvecback
        SAFECLASS(background_at_tau(ba, tau, long_info,
                                    inter_normal, &index, pvecback),
                  ba->error_message);

        d->c->hubble[z_index] = pvecback[ba->index_bg_H]; // 1/Mpc
        d->c->comoving[z_index] = pvecback[ba->index_bg_conf_distance]; // Mpc
        d->c->angular_diameter[z_index] = pvecback[ba->index_bg_ang_distance]; // Mpc, physical
        d->c->Dsq[z_index] = gsl_pow_2(pvecback[ba->index_bg_D]); // squared growth factor
        d->c->Om[z_index] = pvecback[ba->index_bg_Omega_m];
        d->c->rho_c[z_index] = 3.0 * gsl_pow_2(SPEEDOFLIGHT) / 8.0 / M_PI / GNEWTON
                               * gsl_pow_2(d->c->hubble[z_index]);
        d->c->rho_m[z_index] = d->c->Om[z_index] * d->c->rho_c[z_index];
    }

    if (d->p->stype == hmpdf_kappa) // need to compute critical surface density
    {
        if (d->n->dndz != NULL)
        // non-trivial source distribution
        {
            // allocate some integration resources
            gsl_integration_workspace *ws;
            SAFEALLOC(ws, gsl_integration_workspace_alloc(DNDZ_INTEGR_LIMIT));
            gsl_function F;

            // first figure out the normalization
            double norm, err;
            F.function = d->n->dndz;
            F.params = d->n->dndz_params;
            SAFEGSL(gsl_integration_qag(&F, 0.0, d->n->zsource,
                                        0.0, DNDZ_INTEGR_EPSREL,
                                        DNDZ_INTEGR_LIMIT, DNDZ_INTEGR_KEY, ws,
                                        &norm, &err));

            // prepare our integration struct
            dndz_integr_params p = { .ba=ba, .pvecback=pvecback, .index=0,
                                     .dndz=d->n->dndz, .dndz_params=d->n->dndz_params,
                                     .status=0 };
            F.function = dndz_integr_f;
            F.params = &p;

            // now compute the critical surface densities
            for (int z_index=0; z_index<d->n->Nz; z_index++)
            {
                p.chi_z = d->c->comoving[z_index];
                double out;
                SAFEGSL(gsl_integration_qag(&F, d->n->zgrid[z_index], d->n->zsource,
                                            DNDZ_INTEGR_EPSABS*norm, DNDZ_INTEGR_EPSREL,
                                            DNDZ_INTEGR_LIMIT, DNDZ_INTEGR_KEY, ws,
                                            &out, &err));

                HMPDFCHECK(p.status, "error encountered during integration");

                d->c->invScrit[z_index] = 4.0*M_PI*GNEWTON*p.chi_z
                                          /gsl_pow_2(SPEEDOFLIGHT)/(1.0 + d->n->zgrid[z_index])
                                          * out / norm;
            }
            
            // clean up
            gsl_integration_workspace_free(ws);
        }
        else
        // Dirac delta source distribution
        {
            // find distances to source position
            SAFECLASS(background_tau_of_z(ba, d->n->zsource, &tau),
                      ba->error_message);
            SAFECLASS(background_at_tau(ba, tau, long_info,
                                        inter_normal, &index, pvecback),
                      ba->error_message);
            double chi_s = pvecback[ba->index_bg_conf_distance];
            double dA_s = pvecback[ba->index_bg_ang_distance];
            // fill the Scrit grid
            for (int z_index=0; z_index<d->n->Nz; z_index++)
            {
                d->c->invScrit[z_index] = 4.0*M_PI*GNEWTON/gsl_pow_2(SPEEDOFLIGHT)/(1.0+d->n->zsource)
                                          * (chi_s - d->c->comoving[z_index]) * d->c->angular_diameter[z_index]
                                          / dA_s;
            }
        }
    }

    free(pvecback);

    ENDFCT
}//}}}

int
init_cosmology(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(1, "init_cosmo\n");

    SAFEHMPDF(alloc_cosmo(d));
    SAFEHMPDF(fill_background(d));

    ENDFCT
}//}}}

