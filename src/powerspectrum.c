#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_sf_bessel.h>

#include "utils.h"
#include "configs.h"
#include "data.h"
#include "numerics.h"
#include "cosmology.h"
#include "power.h"
#include "profiles.h"
#include "filter.h"
#include "powerspectrum.h"

#include "hmpdf.h"

void null_powerspectrum(all_data *d)
{//{{{
    d->ps->created_Cell = 0;
    d->ps->ell = NULL;
    d->ps->Cell_1h = NULL;
    d->ps->Cell_2h = NULL;
    d->ps->Cell_tot = NULL;
    d->ps->created_Cphi = 0;
    d->ps->phi = NULL;
    d->ps->Cphi_1h = NULL;
    d->ps->Cphi_2h = NULL;
    d->ps->Cphi_tot = NULL;
}//}}}

void reset_powerspectrum(all_data *d)
{//{{{
    if (d->ps->ell != NULL) { free(d->ps->ell); }
    if (d->ps->Cell_1h != NULL) { free(d->ps->Cell_1h); }
    if (d->ps->Cell_2h != NULL) { free(d->ps->Cell_2h); }
    if (d->ps->Cell_tot != NULL) { free(d->ps->Cell_tot); }
    if (d->ps->phi != NULL) { free(d->ps->phi); }
    if (d->ps->Cphi_1h != NULL) { free(d->ps->Cphi_1h); }
    if (d->ps->Cphi_2h != NULL) { free(d->ps->Cphi_2h); }
    if (d->ps->Cphi_tot != NULL) { free(d->ps->Cphi_tot); }
}//}}}

int find_Nell(all_data *d)
{//{{{
    int N = 0;
    if (d->f->pixelside > 0.0)
    // pixel sidelength is defined, which gives us a natural choice of minimum angle
    {
        N = (int)round(2.0 * d->n->phimax / d->f->pixelside );
        // now check if the Bessel function zero is less than pixelside,
        // otherwise increase N
        while (1)
        {
            double phimin = 2.0 * d->n->phimax * gsl_sf_bessel_zero_J0(1)
                            / gsl_sf_bessel_zero_J0(N);
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

void ps_Mint(all_data *d, int z_index, double *oneh, double *twoh)
// oneh, twoh are d->ps->Nell long
// this function nulls oneh, twoh first
// function includes squaring of the two-halo term,
//     but not application of z-dependent filters
{//{{{
    zero_real(d->ps->Nell, oneh);
    zero_real(d->ps->Nell, twoh);
    double *temp = (double *)malloc(d->ps->Nell * sizeof(double));

    for (int M_index=0; M_index<d->n->NM; M_index++)
    {
        // evaluate the conjugate profile interpolator
        s_of_ell(d, z_index, M_index, d->ps->Nell, d->ps->ell, temp);
        
        // read dndlogM and bias
        double n = d->h->hmf[z_index][M_index];
        double b = d->h->bias[z_index][M_index];

        for (int ii=0; ii<d->ps->Nell; ii++)
        {
            oneh[ii] += n * gsl_pow_2(temp[ii])
                        * d->n->Mweights[M_index];
            twoh[ii] += n * b * temp[ii]
                        * d->n->Mweights[M_index];
        }
    }

    // square the two-halo term
    for (int ii=0; ii<d->ps->Nell; ii++)
    {
        twoh[ii] = gsl_pow_2(twoh[ii]);
    }

    free(temp);
}//}}}

void ps_zint(all_data *d, double *oneh, double *twoh)
// oneh, twoh are d->ps->Nell long
// this function nulls oneh, twoh first
{//{{{
    zero_real(d->ps->Nell, oneh);
    zero_real(d->ps->Nell, twoh);
    double *oneh_z = (double *)malloc(d->ps->Nell * sizeof(double));
    double *twoh_z = (double *)malloc(d->ps->Nell * sizeof(double));

    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        ps_Mint(d, z_index, oneh_z, twoh_z);

        // apply z-dependent filters
        apply_filters(d, d->ps->Nell, d->ps->ell, twoh_z, twoh_z, filter_ps, &z_index);
        apply_filters(d, d->ps->Nell, d->ps->ell, oneh_z, oneh_z, filter_ps, &z_index);

        for (int ii=0; ii<d->ps->Nell; ii++)
        {
            double k = d->ps->ell[ii] / d->c->comoving[z_index];
            #ifdef LOGK
            k = log(k);
            #endif

            oneh[ii] += gsl_pow_2(d->c->comoving[z_index]) / d->c->hubble[z_index]
                        * oneh_z[ii] * d->n->zweights[z_index];
            twoh[ii] += gsl_pow_2(d->c->comoving[z_index]) / d->c->hubble[z_index]
                        * d->c->Dsq[z_index] * Pk_linear(d, k)
                        * twoh_z[ii] * d->n->zweights[z_index];
        }
    }

    free(oneh_z);
    free(twoh_z);
}//}}}

void create_Cell(all_data *d)
{//{{{
    if (d->ps->created_Cell) { return; }
    fprintf(stdout, "\tcreate_Cell\n");
    fflush(stdout);
    d->ps->Nell = PS_NELL;
    d->ps->ell = (double *)malloc(d->ps->Nell * sizeof(double));
    logspace(d->ps->Nell, PS_ELLMIN, PS_ELLMAX, d->ps->ell);
    d->ps->Cell_1h = (double *)malloc(d->ps->Nell * sizeof(double));
    d->ps->Cell_2h = (double *)malloc(d->ps->Nell * sizeof(double));
    d->ps->Cell_tot = (double *)malloc(d->ps->Nell * sizeof(double));
    // compute 1- and 2-halo terms
    ps_zint(d, d->ps->Cell_1h, d->ps->Cell_2h);
    // sum the 1- and 2-halo term
    for (int ii=0; ii<d->ps->Nell; ii++)
    {
        d->ps->Cell_tot[ii] = d->ps->Cell_1h[ii] + d->ps->Cell_2h[ii];
    }
    // filter (z-independent)
    apply_filters(d, d->ps->Nell, d->ps->ell, d->ps->Cell_1h,
                  d->ps->Cell_1h, filter_ps, NULL);
    apply_filters(d, d->ps->Nell, d->ps->ell, d->ps->Cell_2h,
                  d->ps->Cell_2h, filter_ps, NULL);
    apply_filters(d, d->ps->Nell, d->ps->ell, d->ps->Cell_tot,
                  d->ps->Cell_tot, filter_ps, NULL);

    d->ps->created_Cell = 1;
}//}}}

static
void create_Cphi(all_data *d)
{//{{{
    if (d->ps->created_Cphi) { return; }
    fprintf(stdout, "\tcreate_Cphi\n");
    fflush(stdout);
    d->ps->Nell_corr = find_Nell(d);
    gsl_dht *t = gsl_dht_new(d->ps->Nell_corr, 0, 2.0 * d->n->phimax);
    
    d->ps->phi = (double *)malloc(d->ps->Nell_corr * sizeof(double));
    double *temp_ell = (double *)malloc(d->ps->Nell_corr * sizeof(double));
    double *temp_Cell = (double *)malloc(d->ps->Nell_corr * sizeof(double));
    d->ps->Cphi_1h = (double *)malloc(d->ps->Nell_corr * sizeof(double));
    d->ps->Cphi_2h = (double *)malloc(d->ps->Nell_corr * sizeof(double));
    d->ps->Cphi_tot = (double *)malloc(d->ps->Nell_corr * sizeof(double));
    
    for (int ii=0; ii<d->ps->Nell_corr; ii++)
    {
        d->ps->phi[ii] = gsl_dht_x_sample(t, ii);
        temp_ell[ii] = gsl_dht_k_sample(t, ii);
    }

    double hankel_norm = gsl_pow_2(gsl_dht_k_sample(t, 0)
                                   / gsl_dht_x_sample(t, 0));

    // 1halo term
    get_Cell(d, d->ps->Nell_corr, temp_ell, temp_Cell, onehalo);
    gsl_dht_apply(t, temp_Cell, d->ps->Cphi_1h);
    // 2halo term
    get_Cell(d, d->ps->Nell_corr, temp_ell, temp_Cell, twohalo);
    gsl_dht_apply(t, temp_Cell, d->ps->Cphi_2h);
    // total
    get_Cell(d, d->ps->Nell_corr, temp_ell, temp_Cell, total);
    gsl_dht_apply(t, temp_Cell, d->ps->Cphi_tot);
    // normalizetion
    for (int ii=0; ii<d->ps->Nell_corr; ii++)
    {
        d->ps->Cphi_1h[ii] *= 0.5 * M_1_PI * hankel_norm;
        d->ps->Cphi_2h[ii] *= 0.5 * M_1_PI * hankel_norm;
        d->ps->Cphi_tot[ii] *= 0.5 * M_1_PI * hankel_norm;
    }

    free(temp_ell);
    free(temp_Cell);
    gsl_dht_free(t);

    d->ps->created_Cphi = 1;
}//}}}

static
void prepare_Cell(all_data *d)
{//{{{
    fprintf(stdout, "In powerspectrum.h -> prepare_Cell.\n");
    fflush(stdout);
    create_conj_profiles(d);
    create_Cell(d);
}//}}}

static
void prepare_Cphi(all_data *d)
{//{{{
    fprintf(stdout, "In powerspectrum.h -> prepare_Cphi :\n");
    fflush(stdout);
    create_conj_profiles(d);
    create_Cell(d);
    create_Cphi(d);
}//}}}

void get_Cell(all_data *d, int Nell, double *ell, double *Cell, Cell_mode mode)
{//{{{
    prepare_Cell(d);

    double *C;
    switch (mode)
    {
        case onehalo : C = d->ps->Cell_1h; break;
        case twohalo : C = d->ps->Cell_2h; break;
        case total   : C = d->ps->Cell_tot; break;
        default      : fprintf(stderr, "Invalid Cell_mode in get_Cell.\n"); return;
                       fflush(stderr);
    }
    interp1d *interp = new_interp1d(d->ps->Nell, d->ps->ell, C, C[0], 0.0, CELL_INTERP_TYPE, NULL);

    for (int ii=0; ii<Nell; ii++)
    {
        Cell[ii] = interp1d_eval(interp, ell[ii]);
    }

    delete_interp1d(interp);
}//}}}

void get_Cphi(all_data *d, int Nphi, double *phi, double *Cphi, Cell_mode mode)
{//{{{
    prepare_Cphi(d);

    double *C;
    switch (mode)
    {
        case onehalo : C = d->ps->Cphi_1h; break;
        case twohalo : C = d->ps->Cphi_2h; break;
        case total   : C = d->ps->Cphi_tot; break;
        default      : fprintf(stderr, "Invalid Cphi_mode in get_Cphi.\n"); return;
                       fflush(stderr);
    }
    interp1d *interp = new_interp1d(d->ps->Nell_corr, d->ps->phi,
                                    C, C[0], 0.0, CPHI_INTERP_TYPE, NULL);

    for (int ii=0; ii<Nphi; ii++)
    {
        Cphi[ii] = interp1d_eval(interp, phi[ii]);
    }
    
    delete_interp1d(interp);
}//}}}

