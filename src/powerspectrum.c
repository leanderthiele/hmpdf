#include <stdlib.h>
#include <string.h>
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
    d->ps->Cell_noisy = NULL;
    d->ps->created_Cphi = 0;
    d->ps->phi = NULL;
    d->ps->Cphi_1h = NULL;
    d->ps->Cphi_2h = NULL;
    d->ps->Cphi_tot = NULL;
    d->ps->Cphi_noisy = NULL;
    d->ps->created_Covell = 0;
    d->ps->Covell = NULL;
}//}}}

void reset_powerspectrum(all_data *d)
{//{{{
    if (d->ps->ell != NULL) { free(d->ps->ell); }
    if (d->ps->Cell_1h != NULL) { free(d->ps->Cell_1h); }
    if (d->ps->Cell_2h != NULL) { free(d->ps->Cell_2h); }
    if (d->ps->Cell_tot != NULL) { free(d->ps->Cell_tot); }
    if (d->ps->Cell_noisy != NULL) { free(d->ps->Cell_noisy); }
    if (d->ps->phi != NULL) { free(d->ps->phi); }
    if (d->ps->Cphi_1h != NULL) { free(d->ps->Cphi_1h); }
    if (d->ps->Cphi_2h != NULL) { free(d->ps->Cphi_2h); }
    if (d->ps->Cphi_tot != NULL) { free(d->ps->Cphi_tot); }
    if (d->ps->Cphi_noisy != NULL) { free(d->ps->Cphi_noisy); }
    if (d->ps->Covell != NULL) { free(d->ps->Covell); }
}//}}}

static
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

static
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

static
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
        apply_filters(d, d->ps->Nell, d->ps->ell, twoh_z, twoh_z, 1, filter_ps, &z_index);
        apply_filters(d, d->ps->Nell, d->ps->ell, oneh_z, oneh_z, 1, filter_ps, &z_index);

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

static
void add_psnoise(all_data *d, int N, double *Cell)
{//{{{
    double Omegapix = gsl_pow_2(d->f->pixelside);
    for (int ii=0; ii<N; ii++)
    {
        Cell[ii] += Omegapix * gsl_pow_2(d->ns->noise);
    }
}//}}}

static
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
                  d->ps->Cell_1h, 1, filter_ps, NULL);
    apply_filters(d, d->ps->Nell, d->ps->ell, d->ps->Cell_2h,
                  d->ps->Cell_2h, 1, filter_ps, NULL);
    apply_filters(d, d->ps->Nell, d->ps->ell, d->ps->Cell_tot,
                  d->ps->Cell_tot, 1, filter_ps, NULL);

    if (d->ns->noise>0.0 && d->f->pixelside>0.0)
    {
        d->ps->Cell_noisy = (double *)malloc(d->ps->Nell * sizeof(double));
        memcpy(d->ps->Cell_noisy, d->ps->Cell_tot, d->ps->Nell * sizeof(double));
        add_psnoise(d, d->ps->Nell, d->ps->Cell_noisy);
    }

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
    get_Cell(d, d->ps->Nell_corr, temp_ell, temp_Cell, onehalo, 0);
    gsl_dht_apply(t, temp_Cell, d->ps->Cphi_1h);
    // 2halo term
    get_Cell(d, d->ps->Nell_corr, temp_ell, temp_Cell, twohalo, 0);
    gsl_dht_apply(t, temp_Cell, d->ps->Cphi_2h);
    // total
    get_Cell(d, d->ps->Nell_corr, temp_ell, temp_Cell, total, 0);
    gsl_dht_apply(t, temp_Cell, d->ps->Cphi_tot);
    // normalizetion
    for (int ii=0; ii<d->ps->Nell_corr; ii++)
    {
        d->ps->Cphi_1h[ii] *= 0.5 * M_1_PI * hankel_norm;
        d->ps->Cphi_2h[ii] *= 0.5 * M_1_PI * hankel_norm;
        d->ps->Cphi_tot[ii] *= 0.5 * M_1_PI * hankel_norm;
    }

    if (d->ns->noise>0.0 && d->f->pixelside>0.0)
    {
        d->ps->Cphi_noisy = (double *)malloc(d->ps->Nell_corr * sizeof(double));
        get_Cell(d, d->ps->Nell_corr, temp_ell, temp_Cell, total, 1);
        gsl_dht_apply(t, temp_Cell, d->ps->Cphi_noisy);
    }

    free(temp_ell);
    free(temp_Cell);
    gsl_dht_free(t);

    d->ps->created_Cphi = 1;
}//}}}

static
void pscov_Mint(all_data *d, int z_index, double *oneh)
{//{{{
    zero_real(d->ps->Nell * d->ps->Nell, oneh);
    double *temp = (double *)malloc(d->ps->Nell * sizeof(double));

    for (int M_index=0; M_index<d->n->NM; M_index++)
    {
        // evalute conjugate profile interpolator
        s_of_ell(d, z_index, M_index, d->ps->Nell, d->ps->ell, temp);

        // read dndlogM
        double n = d->h->hmf[z_index][M_index];

        for (int ii=0; ii<d->ps->Nell; ii++)
        {
            for (int jj=0; jj<d->ps->Nell; jj++)
            {
                oneh[ii*d->ps->Nell+jj] += n * gsl_pow_2(temp[ii])
                                             * gsl_pow_2(temp[jj])
                                             * d->n->Mweights[M_index];
            }
        }
    }

    free(temp);
}//}}}

static
void pscov_zint(all_data *d, double *oneh)
{//{{{
    zero_real(d->ps->Nell * d->ps->Nell, oneh);
    double *oneh_z = (double *)malloc(d->ps->Nell * d->ps->Nell * sizeof(double));

    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        // perform integration over masses
        pscov_Mint(d, z_index, oneh_z);

        // apply z-dependent filters
        for (int ii=0; ii<d->ps->Nell; ii++)
        {
            // row-wise
            apply_filters(d, d->ps->Nell, d->ps->ell,
                          oneh_z+ii, oneh_z+ii,
                          d->ps->Nell/*stride*/, filter_ps, &z_index);
            // col-wise
            apply_filters(d, d->ps->Nell, d->ps->ell,
                          oneh_z+ii*d->ps->Nell, oneh_z+ii*d->ps->Nell,
                          1/*stride*/, filter_ps, &z_index);
        }

        // add to covariance matrix
        for (int ii=0; ii<d->ps->Nell*d->ps->Nell; ii++)
        {
            oneh[ii] += oneh_z[ii] * d->n->zweights[z_index];
        }
    }

    free(oneh_z);
}//}}}

static
void rescale_to_fsky1(int N, double *cov)
{//{{{
    for (int ii=0; ii<N*N; ii++)
    {
        cov[ii] /= 4.0 * M_PI;
    }
}//}}}

static
void add_shotnoise(int N, double *ell, double *cov, double *ps)
{//{{{
    for (int ii=0; ii<N; ii++)
    {
        cov[ii*(N+1)] += gsl_pow_2(ps[ii]) / (ell[ii]+0.5);
    }
}//}}}

static
void create_Covell(all_data *d)
{//{{{
    if (d->ps->created_Covell) { return; }
    fprintf(stdout, "\tcreate_Covell\n");
    fflush(stdout);

    d->ps->Covell = (double *)malloc(d->ps->Nell * d->ps->Nell * sizeof(double));
    // compute the z-integral
    pscov_zint(d, d->ps->Covell);
    // apply z-independent filters
    for (int ii=0; ii<d->ps->Nell; ii++)
    {
        // row-wise
        apply_filters(d, d->ps->Nell, d->ps->ell,
                      d->ps->Covell+ii, d->ps->Covell+ii,
                      d->ps->Nell/*stride*/, filter_ps, NULL);
        // col-wise
        apply_filters(d, d->ps->Nell, d->ps->ell,
                      d->ps->Covell+ii*d->ps->Nell, d->ps->Covell+ii*d->ps->Nell,
                      1/*stride*/, filter_ps, NULL);
    }

    rescale_to_fsky1(d->ps->Nell, d->ps->Covell);

    d->ps->created_Covell = 1;
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
void prepare_Covell(all_data *d)
{//{{{
    fprintf(stdout, "In powerspectrum.h -> prepare_Covell.\n");
    fflush(stdout);
    create_conj_profiles(d);
    create_Cell(d);
    create_Covell(d);
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

void get_Cell(all_data *d, int Nbins, double *binedges, double *Cell, Cell_mode mode, int noisy)
{//{{{
    if (noisy && mode!=total)
    {
        fprintf(stderr, "Error : Noisy power spectrum only makes sense in total Cell_mode.\n");
        fflush(stderr);
        return;
    }
    if (noisy && d->ns->noise<0.0)
    {
        fprintf(stderr, "Error : You requested noisy power spectrum "
                        "but no/invalid noise level was passed.\n");
        fflush(stderr);
        return;
    }
    if (noisy && d->f->pixelside<0.0)
    {
        fprintf(stderr, "Error : You requested noisy power spectrum "
                        "but no/invalid pixel sidelength was passed.\n");
        fflush(stderr);
        return;
    }

    prepare_Cell(d);
    
    double *C;
    switch (mode)
    {
        case onehalo : C = d->ps->Cell_1h; break;
        case twohalo : C = d->ps->Cell_2h; break;
        case total   : C = (noisy) ? d->ps->Cell_noisy : d->ps->Cell_tot; break;
        default      : fprintf(stderr, "Invalid Cell_mode in get_Cell.\n"); return;
                       fflush(stderr);
    }

    bin_1d(d->ps->Nell, d->ps->ell, C,
           Nbins, binedges, Cell, CELL_INTERP_TYPE);
}//}}}

void get_Cphi(all_data *d, int Nphi, double *phi, double *Cphi, Cell_mode mode, int noisy)
{//{{{
    if (noisy && mode!=total)
    {
        fprintf(stderr, "Error : Noisy correlation function only makes sense in total Cell_mode.\n");
        fflush(stderr);
        return;
    }
    if (noisy && d->ns->noise<0.0)
    {
        fprintf(stderr, "Error : You requested noisy correlation function "
                        "but no/invalid noise level was passed.\n");
        fflush(stderr);
        return;
    }
    if (noisy && d->f->pixelside<0.0)
    {
        fprintf(stderr, "Error : You requested noisy correlation function "
                        "but no/invalid pixel sidelength was passed.\n");
        fflush(stderr);
        return;
    }

    prepare_Cphi(d);

    double *C;
    switch (mode)
    {
        case onehalo : C = d->ps->Cphi_1h; break;
        case twohalo : C = d->ps->Cphi_2h; break;
        case total   : C = (noisy) ? d->ps->Cphi_noisy : d->ps->Cphi_tot; break;
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

void get_Cell_cov(all_data *d, int Nbins, double *binedges, double *Covell, int noisy)
{//{{{
    if (noisy && d->ns->noise<0.0)
    {
        fprintf(stderr, "Error : You requested noisy Cell covariance matrix "
                        "but no/invalid noise level was passed.\n");
        fflush(stderr);
        return;
    }
    if (noisy && d->f->pixelside<0.0)
    {
        fprintf(stderr, "Error : You requested noisy Cell covariance matrix "
                        "but no/invalid pixel sidelength was passed.\n");
        fflush(stderr);
        return;
    }

    // perform the computation
    prepare_Covell(d);

    bin_2d(d->ps->Nell, d->ps->ell, d->ps->Covell, PS_COVINTEGR_N,
           Nbins, binedges, Covell, interp2d_bilinear);

    double *temp = (double *)malloc(Nbins * sizeof(double));
    get_Cell(d, Nbins, binedges, temp, total, noisy);
    double *bincentres = (double *)malloc(Nbins * sizeof(double));
    // normalize correctly, and compute bin centres
    for (int ii=0; ii<Nbins; ii++)
    {
        temp[ii] *= binedges[ii+1] - binedges[ii];
        bincentres[ii] = 0.5 * (binedges[ii+1] + binedges[ii]);
    }
    add_shotnoise(Nbins, bincentres, Covell, temp);

    free(temp);
    free(bincentres);
}//}}}
