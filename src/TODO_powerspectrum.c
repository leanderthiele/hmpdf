#ifndef POWERSPECTRUM_H
#define POWERSPECTRUM_H

#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_sf_bessel.h>

#include "utils.h"
#include "data.h"
#include "numerics.h"
#include "cosmology.h"
#include "power.h"
#include "filter.h"

#ifdef XXX
void ps_Mint(all_data *d, int z_index, double *oneh, double *twoh)
// oneh, twoh are d->ps->Nell long
// this function nulls oneh, twoh first
{//{{{
    double *oneh_M = (double *)malloc(d->ps->Nell * sizeof(double)
                                      * ((d->n->monotonize) ? 1 : d->n->gr->NM));
    double *twoh_M = (double *)malloc(d->ps->Nell * sizeof(double)
                                      * ((d->n->monotonize) ? 1 : d->n->gr->NM));
    if (d->n->monotonize)
    {
        for (int ii=0; ii<d->ps->Nell; ii++)
        {
            oneh[ii] = twoh[ii] = 0.0;
        }
    }

    for (int M_index=0; M_index<d->n->gr->NM; M_index++)
    {
        int start = (d->n->monotonize) ? 0 : M_index*d->ps->Nell;
        // evaluate the conjugate profile interpolator
        s_of_l(d, z_index, M_index, d->ps->Nell, d->ps->ell, twoh_M+start);
        double n, b;
        n = dndlogM(d, z_index, M_index, &b);
        for (int ii=0; ii<d->ps->Nell; ii++)
        {
            oneh_M[start+ii] = n * gsl_pow_2(twoh_M[start+ii]);
            twoh_M[start+ii] *= n * b;
        }

        if (d->n->monotonize)
        {
            for (int ii=0; ii<d->ps->Nell; ii++)
            {
                oneh[ii] += oneh_M[ii] * d->n->gr->Mweights[M_index];
                twoh[ii] += twoh_M[ii] * d->n->gr->Mweights[M_index];
            }
        }
    }

    // perform mass integrations if we have to
    if (!(d->n->monotonize))
    {
        for (int ii=0; ii<d->ps->Nell; ii++)
        {
            oneh[ii] = integr_real(d->n->gr->NM, d->n->gr->dlogM, d->ps->Nell, oneh_M+ii);
            twoh[ii] = integr_real(d->n->gr->NM, d->n->gr->dlogM, d->ps->Nell, twoh_M+ii);
        }
    }

    free(oneh_M);
    free(twoh_M);
}//}}}

void ps_zint(all_data *d, double *oneh, double *twoh)
// oneh, twoh are d->ps->Nell long
// this function nulls oneh, twoh first
{//{{{
    double *oneh_z = (double *)malloc(d->ps->Nell * sizeof(double)
                                      * ((d->n->monotonize) ? 1 : d->n->gr->Nz));
    double *twoh_z = (double *)malloc(d->ps->Nell * sizeof(double)
                                      * ((d->n->monotonize) ? 1 : d->n->gr->Nz));
    if (d->n->monotonize)
    {
        for (int ii=0; ii<d->ps->Nell; ii++)
        {
            oneh[ii] = twoh[ii] = 0.0;
        }
    }

    for (int z_index=0; z_index<d->n->gr->Nz; z_index++)
    {
        int start = (d->n->monotonize) ? 0 : z_index*d->ps->Nell;
        ps_Mint(d, z_index, oneh_z+start, twoh_z+start);

        for (int ii=0; ii<d->ps->Nell; ii++)
        {
            double k = d->ps->ell[ii] / comoving(d, z_index);
            #ifdef LOGK
            k = log(k);
            #endif
            oneh_z[start+ii] *= gsl_pow_2(comoving(d, z_index)) / hubble(d, z_index);
            twoh_z[start+ii] = gsl_pow_2(comoving(d, z_index)) / hubble(d, z_index)
                               * Dsq(d, z_index) * Pk_linear(d, k)
                               * gsl_pow_2(twoh_z[start+ii]);
        }

        if (d->n->monotonize)
        {
            for (int ii=0; ii<d->ps->Nell; ii++)
            {
                oneh[ii] += oneh_z[ii] * d->n->gr->zweights[z_index];
                twoh[ii] += twoh_z[ii] * d->n->gr->zweights[z_index];
            }
        }
    }

    if (!(d->n->monotonize))
    {
        for (int ii=0; ii<d->ps->Nell; ii++)
        {
            oneh[ii] = integr_real(d->n->gr->Nz, d->n->gr->dz, d->ps->Nell, oneh_z+ii);
            twoh[ii] = integr_real(d->n->gr->Nz, d->n->gr->dz, d->ps->Nell, twoh_z+ii);
        }
    }
}//}}}

int find_Nell(all_data *d)
{//{{{
    int N = 0;
    if (d->f->pixelside > 0.0)
    // pixel sidelength is defined, which gives us a natural choice of minimum angle
    {
        N = (int)round(d->n->gr->phimax / d->f->pixelside );
        // now check if the Bessel function zero is less than pixelside,
        // otherwise increase N
        while (1)
        {
            double phimin = d->n->gr->phimax * gsl_sf_bessel_zero_J0(1)
                            /gsl_sf_bessel_zero_J0(N);
            if (phimin < d->f->pixelside)
            {
                break;
            }
            else
            // I don't actually think this path is ever taken given what the
            // zeros look like, but anyway.
            {
                N = (int)round(1.1 * (double)N);
            }
        }
    }
    N = GSL_MAX(N, PS_NELL);
    return N;
}//}}}

void create_Cell(all_data *d)
{//{{{
    if (d->ps->created_Cell) { return; }
    printf("\tcreate_Cell\n");
    d->ps->Nell = find_Nell(d);
    d->ps->dht_ws = gsl_dht_new(d->ps->Nell, 0, d->n->gr->phimax);
    d->ps->ell = (double *)malloc(d->ps->Nell * sizeof(double));
    d->ps->Cell_1h = (double *)malloc(d->ps->Nell * sizeof(double));
    d->ps->Cell_2h = (double *)malloc(d->ps->Nell * sizeof(double));
    d->ps->Cell_tot = (double *)malloc(d->ps->Nell * sizeof(double));
    for (int ii=0; ii<d->ps->Nell; ii++)
    {
        d->ps->ell[ii] = gsl_dht_k_sample(d->ps->dht_ws, ii);
    }
    // compute 1- and 2-halo terms
    ps_zint(d, d->ps->Cell_1h, d->ps->Cell_2h);
    for (int ii=0; ii<d->ps->Nell; ii++)
    {
        d->ps->Cell_tot[ii] = d->ps->Cell_1h[ii] + d->ps->Cell_2h[ii];
    }
    // filter
    apply_filters(d, d->ps->Nell, d->ps->ell, d->ps->Cell_1h, d->ps->Cell_1h, filter_ps);
    apply_filters(d, d->ps->Nell, d->ps->ell, d->ps->Cell_2h, d->ps->Cell_2h, filter_ps);
    apply_filters(d, d->ps->Nell, d->ps->ell, d->ps->Cell_tot, d->ps->Cell_tot, filter_ps);

    d->ps->created_Cell = 1;

    // TODO TEST
    savetxt("Cell_dblres.dat", d->ps->Nell, 4, d->ps->ell, d->ps->Cell_1h, d->ps->Cell_2h, d->ps->Cell_tot);
}//}}}

void create_Cphi(all_data *d)
{//{{{
    if (d->ps->created_Cphi) { return; }
    printf("\tcreate_Cphi\n");
    d->ps->phi = (double *)malloc(d->ps->Nell * sizeof(double));
    d->ps->Cphi = (double *)malloc(d->ps->Nell * sizeof(double));
    for (int ii=0; ii<d->ps->Nell; ii++)
    {
        d->ps->phi[ii] = gsl_dht_x_sample(d->ps->dht_ws, ii);
    }
    double hankel_norm = gsl_pow_2(gsl_dht_k_sample(d->ps->dht_ws, 0)
                                   / gsl_dht_x_sample(d->ps->dht_ws, 0));
    gsl_dht_apply(d->ps->dht_ws, d->ps->Cell_tot, d->ps->Cphi);
    for (int ii=0; ii<d->ps->Nell; ii++)
    {
        d->ps->Cphi[ii] *= 0.5 * M_1_PI * hankel_norm;
    }
    d->ps->Cphi_interp = gsl_interp_alloc(gsl_interp_cspline, d->ps->Nell);
    d->ps->Cphi_accel = gsl_interp_accel_alloc();
    gsl_interp_init(d->ps->Cphi_interp, d->ps->phi, d->ps->Cphi, d->ps->Nell);

    d->ps->created_Cphi = 1;
}//}}}

void prepare_Cell(all_data *d)
{//{{{
    printf("In powerspectrum.h -> prepare_Cell.\n");
    create_conj_profiles(d);
    create_Cell(d);
}//}}}

void prepare_Cphi(all_data *d)
{//{{{
    printf("In powerspectrum.h -> prepare_Cphi :\n");
    prepare_Cell(d);
    create_Cphi(d);
}//}}}

typedef enum
{
    onehalo,
    twohalo,
    total,
}
Cell_mode;

void get_Cell(all_data *d, int Nell, double *ell, double *Cell, Cell_mode mode)
{//{{{
    double *C;
    switch (mode)
    {
        case onehalo : C = d->ps->Cell_1h; break;
        case twohalo : C = d->ps->Cell_2h; break;
        case total   : C = d->ps->Cell_tot; break;
        default      : printf("Invalid Cell_mode in get_Cell.\n"); return;
    }
    interp1d *interp = new_interp1d(d->ps->Nell, d->ps->ell, C, C[0], 0.0, CELL_INTERP_TYPE, NULL);

    for (int ii=0; ii<Nell; ii++)
    {
        Cell[ii] = interp1d_eval(interp, ell[ii]);
    }

    delete_interp1d(interp);
}//}}}

void get_Cphi(all_data *d, int Nphi, double *phi, double *Cphi)
{//{{{
    for (int ii=0; ii<Nphi; ii++)
    {
        double p = phi[ii];
        if (p < d->ps->phi[0])
        {
            Cphi[ii] = d->ps->Cphi[0];
        }
        else if (p > d->ps->phi[d->ps->Nell-1])
        {
            Cphi[ii] = 0.0;
        }
        else
        {
            Cphi[ii] = gsl_interp_eval(d->ps->Cphi_interp, d->ps->phi, d->ps->Cphi,
                                       p, d->ps->Cphi_accel);
        }
    }
}//}}}
#endif

#endif
