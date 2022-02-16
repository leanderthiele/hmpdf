#include <stdio.h>
#include <math.h>
#ifdef _OPENMP
#   include <omp.h>
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_sf_expint.h>

#include "configs.h"
#include "utils.h"
#include "object.h"
#include "cosmology.h"
#include "halo_model.h"
#include "bcm.h"
#include "power.h"
#include "pk.h"

#include "hmpdf.h"

int
null_pk(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->pk->inited_pk = 0;
    d->pk->dht_kgrid = NULL;
    d->pk->dht_rgrid = NULL;
    d->pk->dht_ws = NULL;

    ENDFCT
}//}}}

int
reset_pk(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->pk->dht_kgrid != NULL) { free(d->pk->dht_kgrid); }
    if (d->pk->dht_rgrid != NULL) { free(d->pk->dht_rgrid); }
    if (d->pk->dht_ws != NULL) { gsl_dht_free(d->pk->dht_ws); }

    ENDFCT
}//}}}

int
init_pk(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->pk->inited_pk) { return 0; }

    // TODO make customizable or at least confirm that these are good choices
    d->pk->dht_Nk = 256;

    SAFEALLOC(d->pk->dht_kgrid, malloc(d->pk->dht_Nk * sizeof(double)));
    SAFEALLOC(d->pk->dht_rgrid, malloc(d->pk->dht_Nk * sizeof(double)));
    SAFEALLOC(d->pk->dht_ws, gsl_dht_new(d->pk->dht_Nk, 0.5, 1.0));

    for (int ii=0; ii<d->pk->dht_Nk; ii++)
    {
        d->pk->dht_rgrid[ii] = gsl_dht_x_sample(d->pk->dht_ws, ii);
        d->pk->dht_kgrid[ii] = gsl_dht_k_sample(d->pk->dht_ws, ii);
    }

    d->pk->inited_pk = 1;

    ENDFCT
}//}}}

static int
profile(hmpdf_obj *d, int z_index, int M_index, int N, double *r, double *out)
{//{{{
    STARTFCT
    
    HMPDFCHECK(d->p->stype == hmpdf_tsz, "tSZ not supported at the moment");

    if (d->bcm->Arico20_params == NULL)
    // do standard NFW
    {
        double rhos, rs;
        SAFEHMPDF(NFW_fundamental(d, z_index, M_index, 1.0, &rhos, &rs));

        for (int ii=0; ii<N; ii++)
        {
            double x = r[ii] / rs;
            out[ii] = rhos/(x * gsl_pow_2(1.0+x));
        }
    }
    else
    // do BCM
    {
        bcm_ws *ws = d->bcm->ws[THIS_THREAD];
        SAFEHMPDF(bcm_init_ws(d, z_index, M_index, 1.0, ws));

        for (int ii=0; ii<N; ii++)
            SAFEHMPDF(bcm_density_profile(d, ws, r[ii], out+ii));
    }

    ENDFCT
}//}}}

static double
uk_nfw(double rhos, double rs, double Rout, double k)
{//{{{
    double kR  = k * Rout;
    double krs = k * rs;
    double kp  = k * (Rout + rs);
    return 4.0 * M_PI * gsl_pow_3(rs) * rhos
           * (  cos(krs) * ( gsl_sf_Ci(kp) - gsl_sf_Ci(krs) )
              + sin(krs) * ( gsl_sf_Si(kp) - gsl_sf_Si(krs) )
              - sin(kR) / kp );
}//}}}

static int
pk_at_z(hmpdf_obj *d, int z_index, int N, double *k, double *pk_1h, double *pk_2h)
{//{{{
    STARTFCT

    // zero the outputs
    zero_real(N, pk_1h);
    zero_real(N, pk_2h);

    const int uk_analytic = d->bcm->Arico20_params == NULL;

    #ifdef _OPENMP
    #   pragma omp parallel for num_threads(d->Ncores) schedule(static)
    #endif
    for (int M_index=0; M_index<d->n->NM; M_index++)
    {
        CONTINUE_IF_ERR

        // need to use if Hankel transform is not analytic
        double *rgrid = NULL;
        double *profile_real = NULL;
        double *profile_reci = NULL;
        interp1d *i = NULL;

        // compute the outer radius which we use to normalize
        double M, Rout, c;
        SAFEHMPDF_NORETURN(Mconv(d, z_index, M_index, d->p->rout_def, 1.0, &M, &Rout, &c));
        CONTINUE_IF_ERR
        Rout *= d->p->rout_scale;

        // correction factor used if !uk_analytic to rescale the analytic part of the profile
        double correction_factor = 1.0;

        double rhos, rs;

        SAFEHMPDF_NORETURN(NFW_fundamental(d, z_index, M_index, 1.0, &rhos, &rs));
        CONTINUE_IF_ERR

        if (!uk_analytic)
        // we need to do it numerically
        {
            SAFEALLOC_NORETURN(rgrid, malloc(d->pk->dht_Nk * sizeof(double)));
            SAFEALLOC_NORETURN(profile_real, malloc(d->pk->dht_Nk * sizeof(double)));
            SAFEALLOC_NORETURN(profile_reci, malloc(d->pk->dht_Nk * sizeof(double)));
            CONTINUE_IF_ERR

            for (int ii=0; ii<d->pk->dht_Nk; ii++)
                rgrid[ii] = d->pk->dht_rgrid[ii] * Rout;

            SAFEHMPDF_NORETURN(profile(d, z_index, M_index, d->pk->dht_Nk, rgrid, profile_real));
            CONTINUE_IF_ERR

            // for the DHT, we need to normalize the argument function
            for (int ii=0; ii<d->pk->dht_Nk; ii++)
                profile_real[ii] *= sqrt(d->pk->dht_rgrid[ii]);

            // now perform the DHT to get to k-space (note that dht_ws is const here so thread safe)
            SAFEGSL_NORETURN(gsl_dht_apply(d->pk->dht_ws, profile_real, profile_reci));
            CONTINUE_IF_ERR

            // now normalize properly
            for (int ii=0; ii<d->pk->dht_Nk; ii++)
                profile_reci[ii] *= gsl_pow_3(M_SQRT2 * M_SQRTPI * Rout)
                                    / sqrt(d->pk->dht_kgrid[ii]);

            // interpolate to get to the global k points
            SAFEHMPDF_NORETURN(new_interp1d(d->pk->dht_Nk, d->pk->dht_kgrid, profile_reci,
                                            profile_reci[0], 0.0, interp_steffen, NULL, &i));
            CONTINUE_IF_ERR

            // compute the rescaling factor that makes the profile continuous (TODO play with this)
            correction_factor = profile_reci[0] / uk_nfw(rhos, rs, Rout, d->pk->dht_kgrid[0]/Rout);
        }

        for (int ii=0; ii<N; ii++)
        {
            CONTINUE_IF_ERR

            double uk;

            // TODO not entirely sure if this is correct but seems very likely that our internal
            //      halo quantities are in physical Mpc
            double kphys = k[ii] * (1.0 + d->n->zgrid[z_index]);

            if (!uk_analytic
                && kphys*Rout>d->pk->dht_kgrid[0]
                && kphys*Rout<d->pk->dht_kgrid[d->pk->dht_Nk-1])
            // use the interpolator. For low k outside the interpolated region we revert back to NFW
            {
                SAFEHMPDF_NORETURN(interp1d_eval(i, kphys*Rout, &uk));
                CONTINUE_IF_ERR
            }
            else
                uk = uk_nfw(rhos, rs, Rout, kphys) * correction_factor;

            double this_1h = d->n->Mweights[M_index] * d->h->hmf[z_index][M_index]
                             * gsl_pow_2(uk / d->c->rho_m_0);
            double this_2h = d->n->Mweights[M_index] * d->h->hmf[z_index][M_index]
                             * uk / d->c->rho_m_0
                             * d->h->bias[z_index][M_index];

            #ifdef _OPENMP
            #   pragma omp atomic
            #endif
            pk_1h[ii] += this_1h;

            #ifdef _OPENMP
            #   pragma omp atomic
            #endif
            pk_2h[ii] += this_2h;
        }

        if (rgrid != NULL) free(rgrid);
        if (profile_real != NULL) free(profile_real);
        if (profile_reci != NULL) free(profile_reci);
        if (i != NULL) delete_interp1d(i);
    }

    // complete the two-halo term calculation
    for (int ii=0; ii<N; ii++)
    {
        double Plin;

        #ifdef LOGK
        double this_k = log(k[ii]);
        #else
        double this_k = k[ii];
        #endif

        SAFEHMPDF(Pk_linear(d, this_k, &Plin));

        pk_2h[ii] = gsl_pow_2(pk_2h[ii]) * Plin;
    }

    ENDFCT
}//}}}

int
hmpdf_get_pk(hmpdf_obj *d, double z, int Npoints, double k[Npoints],
             double Pk_1h[Npoints], double Pk_2h[Npoints])
{//{{{
    STARTFCT

    CHECKINIT;

    HMPDFCHECK(z<d->n->zgrid[0] || z>d->n->zgrid[d->n->Nz-1],
               "z outside possible region");

    // zero the output
    zero_real(Npoints, Pk_1h);
    zero_real(Npoints, Pk_2h);

    // find which z indices bracket the requested redshift
    int z_index_lo;
    for (z_index_lo=0; d->n->zgrid[z_index_lo]>z; z_index_lo++);
    int z_indices[2] = { z_index_lo, z_index_lo+1 };
    double z_points[2] = { d->n->zgrid[z_indices[0]], d->n->zgrid[z_indices[1]] };
    double delta_z = z_points[1] - z_points[0];
    double z_weights[2] = { (z_points[1]-z)/delta_z, (z-z_points[0])/delta_z };

    double *buffer_1h;
    double *buffer_2h;
    SAFEALLOC(buffer_1h,  malloc(Npoints * sizeof(double)));
    SAFEALLOC(buffer_2h,  malloc(Npoints * sizeof(double)));

    // compute the power spectra and take weighted average
    for (int ii=0; ii<2; ii++)
    {
        SAFEHMPDF(pk_at_z(d, z_indices[ii], Npoints, k, buffer_1h, buffer_2h));
        for (int jj=0; jj<Npoints; jj++)
        {
            Pk_1h[jj] += buffer_1h[jj] * z_weights[ii];
            Pk_2h[jj] += buffer_2h[jj] * z_weights[ii];
        }
    }

    free(buffer_1h);
    free(buffer_2h);

    ENDFCT
}//}}}
