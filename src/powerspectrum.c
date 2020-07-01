#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_result.h>

#include "utils.h"
#include "configs.h"
#include "object.h"
#include "numerics.h"
#include "cosmology.h"
#include "power.h"
#include "profiles.h"
#include "filter.h"
#include "powerspectrum.h"

#include "hmpdf.h"

int
null_powerspectrum(hmpdf_obj *d)
{//{{{
    STARTFCT

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

    ENDFCT
}//}}}

int
reset_powerspectrum(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\treset_powerspectrum\n")

    if (d->ps->ell != NULL) { free(d->ps->ell); }
    if (d->ps->Cell_1h != NULL) { free(d->ps->Cell_1h); }
    if (d->ps->Cell_2h != NULL) { free(d->ps->Cell_2h); }
    if (d->ps->Cell_tot != NULL) { free(d->ps->Cell_tot); }
    if (d->ps->phi != NULL) { free(d->ps->phi); }
    if (d->ps->Cphi_1h != NULL) { free(d->ps->Cphi_1h); }
    if (d->ps->Cphi_2h != NULL) { free(d->ps->Cphi_2h); }
    if (d->ps->Cphi_tot != NULL) { free(d->ps->Cphi_tot); }

    ENDFCT
}//}}}

int
find_Nell(hmpdf_obj *d, int *out)
{//{{{
    STARTFCT

    int N = 0;
    if (d->f->pixelside > 0.0)
    // pixel sidelength is defined, which gives us a natural choice of minimum angle
    {
        N = (int)round(2.0 * d->n->phimax / d->f->pixelside );
        // now check if the Bessel function zero is less than pixelside,
        // otherwise increase N
        while (1)
        {
            gsl_sf_result result1, resultN;
            SAFEGSL(gsl_sf_bessel_zero_J0_e(1, &result1))
            SAFEGSL(gsl_sf_bessel_zero_J0_e(N, &resultN))
            double phimin = 2.0 * d->n->phimax * result1.val
                            / resultN.val;
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
    *out = N;

    ENDFCT
}//}}}

int
ps_Mint(hmpdf_obj *d, int z_index, double *oneh, double *twoh)
// oneh, twoh are d->ps->Nell long
// this function nulls oneh, twoh first
// function includes squaring of the two-halo term,
//     but not application of z-dependent filters
{//{{{
    STARTFCT

    zero_real(d->ps->Nell, oneh);
    zero_real(d->ps->Nell, twoh);
    SAFEALLOC(double *, temp, malloc(d->ps->Nell * sizeof(double)))

    for (int M_index=0; M_index<d->n->NM; M_index++)
    {
        // evaluate the conjugate profile interpolator
        SAFEHMPDF(s_of_ell(d, z_index, M_index,
                           d->ps->Nell, d->ps->ell, temp))
        
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

    ENDFCT
}//}}}

int
ps_zint(hmpdf_obj *d, double *oneh, double *twoh)
// oneh, twoh are d->ps->Nell long
// this function nulls oneh, twoh first
{//{{{
    STARTFCT

    zero_real(d->ps->Nell, oneh);
    zero_real(d->ps->Nell, twoh);
    SAFEALLOC(double *, oneh_z, malloc(d->ps->Nell * sizeof(double)))
    SAFEALLOC(double *, twoh_z, malloc(d->ps->Nell * sizeof(double)))

    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        SAFEHMPDF(ps_Mint(d, z_index, oneh_z, twoh_z))

        // apply z-dependent filters
        SAFEHMPDF(apply_filters(d, d->ps->Nell, d->ps->ell,
                                twoh_z, twoh_z, 1, filter_ps, &z_index))
        SAFEHMPDF(apply_filters(d, d->ps->Nell, d->ps->ell,
                                oneh_z, oneh_z, 1, filter_ps, &z_index))

        for (int ii=0; ii<d->ps->Nell; ii++)
        {
            double k = d->ps->ell[ii] / d->c->comoving[z_index];
            #ifdef LOGK
            k = log(k);
            #endif

            oneh[ii] += gsl_pow_2(d->c->comoving[z_index]) / d->c->hubble[z_index]
                        * oneh_z[ii] * d->n->zweights[z_index];
            double Pk;
            SAFEHMPDF(Pk_linear(d, k, &Pk))
            twoh[ii] += gsl_pow_2(d->c->comoving[z_index]) / d->c->hubble[z_index]
                        * d->c->Dsq[z_index] * Pk
                        * twoh_z[ii] * d->n->zweights[z_index];
        }
    }

    free(oneh_z);
    free(twoh_z);

    ENDFCT
}//}}}

int
create_Cell(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->ps->created_Cell) { return hmpdf_status; }

    HMPDFPRINT(2, "\tcreate_Cell\n")
    
    d->ps->Nell = PS_NELL;
    SAFEALLOC(, d->ps->ell, malloc(d->ps->Nell * sizeof(double)))
    logspace(d->ps->Nell, PS_ELLMIN, PS_ELLMAX, d->ps->ell);
    SAFEALLOC(, d->ps->Cell_1h,  malloc(d->ps->Nell * sizeof(double)))
    SAFEALLOC(, d->ps->Cell_2h,  malloc(d->ps->Nell * sizeof(double)))
    SAFEALLOC(, d->ps->Cell_tot, malloc(d->ps->Nell * sizeof(double)))
    // compute 1- and 2-halo terms
    SAFEHMPDF(ps_zint(d, d->ps->Cell_1h, d->ps->Cell_2h))
    // sum the 1- and 2-halo term
    for (int ii=0; ii<d->ps->Nell; ii++)
    {
        d->ps->Cell_tot[ii] = d->ps->Cell_1h[ii] + d->ps->Cell_2h[ii];
    }
    // filter (z-independent)
    SAFEHMPDF(apply_filters(d, d->ps->Nell, d->ps->ell, d->ps->Cell_1h,
                            d->ps->Cell_1h, 1, filter_ps, NULL))
    SAFEHMPDF(apply_filters(d, d->ps->Nell, d->ps->ell, d->ps->Cell_2h,
                            d->ps->Cell_2h, 1, filter_ps, NULL))
    SAFEHMPDF(apply_filters(d, d->ps->Nell, d->ps->ell, d->ps->Cell_tot,
                            d->ps->Cell_tot, 1, filter_ps, NULL))

    d->ps->created_Cell = 1;

    ENDFCT
}//}}}

static int
create_Cphi(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->ps->created_Cphi) { return hmpdf_status; }

    HMPDFPRINT(2, "\tcreate_Cphi\n")
    
    SAFEHMPDF(find_Nell(d, &(d->ps->Nell_corr)))
    gsl_dht *t = gsl_dht_new(d->ps->Nell_corr, 0, 2.0 * d->n->phimax);
    
    SAFEALLOC(, d->ps->phi, malloc(d->ps->Nell_corr * sizeof(double)))
    SAFEALLOC(double *, temp_ell,  malloc(d->ps->Nell_corr * sizeof(double)))
    SAFEALLOC(double *, temp_Cell, malloc(d->ps->Nell_corr * sizeof(double)))
    SAFEALLOC(, d->ps->Cphi_1h,  malloc(d->ps->Nell_corr * sizeof(double)))
    SAFEALLOC(, d->ps->Cphi_2h,  malloc(d->ps->Nell_corr * sizeof(double)))
    SAFEALLOC(, d->ps->Cphi_tot, malloc(d->ps->Nell_corr * sizeof(double)))
    
    for (int ii=0; ii<d->ps->Nell_corr; ii++)
    {
        d->ps->phi[ii] = gsl_dht_x_sample(t, ii);
        temp_ell[ii] = gsl_dht_k_sample(t, ii);
    }

    double hankel_norm = gsl_pow_2(gsl_dht_k_sample(t, 0)
                                   / gsl_dht_x_sample(t, 0));

    // 1halo term
    SAFEHMPDF(hmpdf_get_Cell(d, d->ps->Nell_corr, temp_ell,
                             temp_Cell, hmpdf_onehalo))
    SAFEGSL(gsl_dht_apply(t, temp_Cell, d->ps->Cphi_1h))
    // 2halo term
    SAFEHMPDF(hmpdf_get_Cell(d, d->ps->Nell_corr, temp_ell,
                             temp_Cell, hmpdf_twohalo))
    SAFEGSL(gsl_dht_apply(t, temp_Cell, d->ps->Cphi_2h))
    // total
    SAFEHMPDF(hmpdf_get_Cell(d, d->ps->Nell_corr, temp_ell,
                             temp_Cell, hmpdf_total))
    SAFEGSL(gsl_dht_apply(t, temp_Cell, d->ps->Cphi_tot))
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

    ENDFCT
}//}}}

static int
prepare_Cell(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(1, "prepare_Cell\n")
    
    SAFEHMPDF(create_conj_profiles(d))
    SAFEHMPDF(create_Cell(d))

    ENDFCT
}//}}}

static int
prepare_Cphi(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(1, "prepare_Cphi")
    
    SAFEHMPDF(create_conj_profiles(d))
    SAFEHMPDF(create_Cell(d))
    SAFEHMPDF(create_Cphi(d))

    ENDFCT
}//}}}

int
hmpdf_get_Cell(hmpdf_obj *d, int Nell, double ell[Nell], double Cell[Nell], hmpdf_Cell_mode_e mode)
{//{{{
    STARTFCT

    SAFEHMPDF(prepare_Cell(d))

    double *C;
    switch (mode)
    {
        case hmpdf_onehalo : C = d->ps->Cell_1h;
                             break;
        case hmpdf_twohalo : C = d->ps->Cell_2h;
                             break;
        case hmpdf_total   : C = d->ps->Cell_tot;
                             break;
        default            : C = NULL;
                             HMPDFERR("invalid Cell_mode.")
    }
    interp1d *interp;
    SAFEHMPDF(new_interp1d(d->ps->Nell, d->ps->ell, C,
                           C[0], 0.0, CELL_INTERP_TYPE, NULL, &interp))

    for (int ii=0; ii<Nell; ii++)
    {
        SAFEHMPDF(interp1d_eval(interp, ell[ii], Cell+ii))
    }

    delete_interp1d(interp);

    ENDFCT
}//}}}

int
hmpdf_get_Cphi(hmpdf_obj *d, int Nphi, double phi[Nphi], double Cphi[Nphi], hmpdf_Cell_mode_e mode)
{//{{{
    STARTFCT

    SAFEHMPDF(prepare_Cphi(d))

    double *C;
    switch (mode)
    {
        case hmpdf_onehalo : C = d->ps->Cphi_1h;
                             break;
        case hmpdf_twohalo : C = d->ps->Cphi_2h;
                             break;
        case hmpdf_total   : C = d->ps->Cphi_tot;
                             break;
        default            : C = NULL;
                             HMPDFERR("invalid Cphi_mode.")
    }
    interp1d *interp;
    SAFEHMPDF(new_interp1d(d->ps->Nell_corr, d->ps->phi,
                           C, C[0], 0.0, CPHI_INTERP_TYPE, NULL, &interp))

    for (int ii=0; ii<Nphi; ii++)
    {
        SAFEHMPDF(interp1d_eval(interp, phi[ii], Cphi+ii))
    }
    
    delete_interp1d(interp);

    ENDFCT
}//}}}

