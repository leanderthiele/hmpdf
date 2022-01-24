#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "utils.h"
#include "configs.h"
#include "object.h"
#include "halo_model.h"

int
null_halo_model(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->h->inited_halo = 0;
    d->h->hmf = NULL;
    d->h->bias = NULL;
    d->h->c_interp = NULL;
    d->h->c_accel = NULL;

    ENDFCT
}//}}}

int
reset_halo_model(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\treset_halo_model\n");

    if (d->h->hmf != NULL)
    {
        for (int z_index=0; z_index<d->n->Nz; z_index++)
        {
            if (d->h->hmf[z_index] != NULL) { free(d->h->hmf[z_index]); }
        }
        free(d->h->hmf);
    }
    if (d->h->bias != NULL)
    {
        for (int z_index=0; z_index<d->n->Nz; z_index++)
        {
            if (d->h->bias[z_index] != NULL) { free(d->h->bias[z_index]); }
        }
        free(d->h->bias);
    }
    if (d->h->c_interp != NULL) { gsl_spline_free(d->h->c_interp); }
    if (d->h->c_accel != NULL)
    {
        for (int ii=0; ii<d->Ncores; ii++)
        {
            if (d->h->c_accel[ii] != NULL)
            {
                gsl_interp_accel_free(d->h->c_accel[ii]);
            }
        }
        free(d->h->c_accel);
    }

    ENDFCT
}//}}}

static int 
DeltaVir_BryanNorman98(hmpdf_obj *d, int z_index, double *out)
{//{{{
    STARTFCT

    double x = d->c->Om[z_index] - 1.0;
    *out = 18.0*M_PI*M_PI + 82.0*x - 39.0*x*x;

    ENDFCT
}//}}}

static int
density_threshold(hmpdf_obj *d, int z_index, hmpdf_mdef_e mdef, double *out)
{//{{{{
    STARTFCT

    double dvir = 0.0; // to avoid maybe-uninitialized
    switch (mdef)
    {
        case hmpdf_mdef_c : *out = 200.0 * d->c->rho_c[z_index];
                            break;
        case hmpdf_mdef_v : SAFEHMPDF(DeltaVir_BryanNorman98(d, z_index, &dvir));
                            *out = dvir * d->c->rho_c[z_index];
                            break;
        case hmpdf_mdef_m : *out = 200.0 * d->c->rho_m[z_index];
                            break;
        default           : *out = 0.0; // to avoid maybe-uninitialized
                            HMPDFERR("Unknown mass definition.");
    }

    ENDFCT
}//}}}

static int 
RofM(hmpdf_obj *d, int z_index, int M_index, double *out,
     double mass_resc)
{//{{{
    STARTFCT

    double dt;
    SAFEHMPDF(density_threshold(d, z_index, MDEF_GLOBAL, &dt));
    *out = cbrt(3.0*mass_resc*d->n->Mgrid[M_index] / 4.0 / M_PI / dt);

    ENDFCT
}//}}}

static int
MofR(hmpdf_obj *d, int z_index, double R, hmpdf_mdef_e mdef, double *out)
{//{{{
    STARTFCT

    double dt;
    SAFEHMPDF(density_threshold(d, z_index, mdef, &dt));
    *out = 4.0 * M_PI * dt * gsl_pow_3(R) / 3.0;

    ENDFCT
}//}}}

static inline double
c_Duffy08_1(hmpdf_obj *d, double z, double M, hmpdf_mdef_e mdef)
{//{{{
    // convert to Msun/h
    M *= d->c->h / 2e12;
    const double *params = d->h->Duffy08_params + (int)(mdef)*3;

    switch (mdef)
    {
        // this case is special, since it has more parameters
        //     (the other definitions are not used for important things)
        case (hmpdf_mdef_m) :
            return params[0]
                   * pow(M,     params[1])
                   * pow(1.0+z, params[2])
                   * exp( params[3] * gsl_pow_2(z-params[4])
                                    * gsl_pow_2(log(M) - params[5]) );
        // use the Duffy08 parameterization
        default :
            return params[0]
                   * pow(M,      params[1])
                   * pow(1.0+z,  params[2]);
    }
}//}}}
static inline double
c_Duffy08(hmpdf_obj *d, int z_index, int M_index,
          double mass_resc)
{//{{{
    double out = c_Duffy08_1(d, d->n->zgrid[z_index],
                             mass_resc * d->n->Mgrid[M_index],
                             MDEF_GLOBAL);

    if (d->h->conc_resc != NULL)
    {
        out *= d->h->conc_resc(d->n->zgrid[z_index],
                               d->n->Mgrid[M_index] * d->c->h,
                               d->h->conc_resc_params);
    }

    return out;
}//}}}

int
NFW_fundamental(hmpdf_obj *d, int z_index, int M_index,
                double mass_resc,
                double *rhos, double *rs)
// returns rhos via function call and rs via return value
// this function is tested against Colossus --> everything here works
{//{{{
    STARTFCT

    double c = c_Duffy08(d, z_index, M_index, mass_resc);
    SAFEHMPDF(RofM(d, z_index, M_index, rs, mass_resc));
    *rs /= c;
    *rhos = mass_resc * d->n->Mgrid[M_index]/4.0/M_PI/gsl_pow_3(*rs)
            / (log1p(c)-c/(1.0+c));

    ENDFCT
}//}}}

static int
create_c_of_y(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\tcreate_c_of_y\n");

    double *logc_grid;
    double *logy_grid;
    SAFEALLOC(logc_grid, malloc(CINTERP_NC * sizeof(double)));
    SAFEALLOC(logy_grid, malloc(CINTERP_NC * sizeof(double)));
    for (int ii=0; ii<CINTERP_NC; ii++)
    {
        logc_grid[ii] = log(CINTERP_CMAX)
                        - (double)(ii)*log(CINTERP_CMAX/CINTERP_CMIN)/(double)(CINTERP_NC-1);
        double _c = exp(logc_grid[ii]);
        logy_grid[ii] = log(3.0) - 3.0*logc_grid[ii] + log(log1p(_c) - _c/(1.0+_c));
    }
    SAFEALLOC(d->h->c_interp, gsl_spline_alloc(gsl_interp_cspline, CINTERP_NC));
    SAFEALLOC(d->h->c_accel,  malloc(d->Ncores * sizeof(gsl_interp_accel *)));
    SETARRNULL(d->h->c_accel, d->Ncores);
    for (int ii=0; ii<d->Ncores; ii++)
    {
        SAFEALLOC(d->h->c_accel[ii], gsl_interp_accel_alloc());
    }
    SAFEGSL(gsl_spline_init(d->h->c_interp, logy_grid, logc_grid, CINTERP_NC));
    free(logc_grid);
    free(logy_grid);

    ENDFCT
}//}}}

static int
c_of_y(hmpdf_obj *d, double y, double *out)
// inverts the function y(c) = 3/c^3 * (log(1+c)-c/(1+c))
// this function is much smoother on loglog scale, so do the interpolation this way
{//{{{
    STARTFCT

    SAFEGSL(gsl_spline_eval_e(d->h->c_interp, log(y),
                              d->h->c_accel[THIS_THREAD], out));
    *out = exp(*out);

    ENDFCT
}//}}}

int
Mconv(hmpdf_obj *d, int z_index, int M_index, hmpdf_mdef_e mdef_out,
      double mass_resc,
      double *M, double *R, double *c)
// returns the converted mass via function call and the new radius and concentration via return value
// this function is tested against Colossus --> everything here works
{//{{{
    STARTFCT

    double rhos, rs;
    SAFEHMPDF(NFW_fundamental(d, z_index, M_index, mass_resc, &rhos, &rs));
    double dt;
    SAFEHMPDF(density_threshold(d, z_index, mdef_out, &dt));
    SAFEHMPDF(c_of_y(d, dt/rhos, c));
    *R = rs * *c;
    SAFEHMPDF(MofR(d, z_index, *R, mdef_out, M));

    ENDFCT
}//}}}

static inline double
fnu_Tinker10_primitive(hmpdf_obj *d, int n, double z)
{//{{{
    return d->h->Tinker10_params[n*2]
           *pow(1.0+z, d->h->Tinker10_params[n*2 + 1]);
}//}}}

static double
fnu_Tinker10(hmpdf_obj *d, double nu, double z)
{//{{{
    z = (z<3.0) ? z : 3.0;
    double beta  = fnu_Tinker10_primitive(d, 0, z);
    double phi   = fnu_Tinker10_primitive(d, 1, z);
    double eta   = fnu_Tinker10_primitive(d, 2, z);
    double gamma = fnu_Tinker10_primitive(d, 3, z);
    double alpha = fnu_Tinker10_primitive(d, 4, z);
    return nu * alpha*(1.0+pow(beta*nu, -2.0*phi))
           * pow(nu, 2.0*eta) * exp(-0.5*gamma*gsl_pow_2(nu));
}//}}}

static double
bnu_Tinker10(double nu)
{//{{{
    double y = 2.0 + M_LN2/M_LN10; // y = log_10(200)
    double A = 1.0 + 0.24 * y * exp(-gsl_pow_4(4.0/y));
    double a = 0.44 * y - 0.88;
    double B = 0.183;
    double b = 1.5;
    double C = 0.019 + 0.107 * y + 0.19 * exp(-gsl_pow_4(4.0/y));
    double c = 2.4;
    return 1.0 - A*pow(nu, a)/(pow(nu, a) + pow(1.686, a)) + B*pow(nu, b) + C*pow(nu, c);
}//}}}

static int
dndlogM(hmpdf_obj *d, int z_index, int M_index, double *hmf, double *bias)
{//{{{
    STARTFCT

    double sigma_squared = d->pwr->ssq[M_index][0];
    double sigma_squared_prime = d->pwr->ssq[M_index][1];
    double nu = 1.686/sqrt(d->c->Dsq[z_index] * sigma_squared);

    double fnu = fnu_Tinker10(d, nu, d->n->zgrid[z_index]);

    *hmf = -fnu * d->c->rho_m_0 * sigma_squared_prime
           / (2.0 * sigma_squared * d->n->Mgrid[M_index]);

    *bias = bnu_Tinker10(nu);

    if (d->h->massfunc_corr != NULL)
    {
        *hmf *= d->h->massfunc_corr(d->n->zgrid[z_index],
                                    d->n->Mgrid[M_index] * d->c->h,
                                    d->h->massfunc_corr_params);
    }

    ENDFCT
}//}}}

static int
create_dndlogM(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\tcreate_dndlogM\n");

    SAFEALLOC(d->h->hmf,  malloc(d->n->Nz * sizeof(double *)));
    SETARRNULL(d->h->hmf, d->n->Nz);
    SAFEALLOC(d->h->bias, malloc(d->n->Nz * sizeof(double *)));
    SETARRNULL(d->h->hmf, d->n->Nz);
    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        SAFEALLOC(d->h->hmf[z_index],  malloc(d->n->NM * sizeof(double)));
        SAFEALLOC(d->h->bias[z_index], malloc(d->n->NM * sizeof(double)));
        for (int M_index=0; M_index<d->n->NM; M_index++)
        {
            // if we are above the mass cut, set HMF to zero
            if (d->n->mass_cuts != NULL
                && d->n->Mgrid[M_index]
                   > d->n->mass_cuts(d->n->zgrid[z_index], d->n->mass_cuts_params)
                     / d->c->h)
            {
                d->h->hmf[z_index][M_index] = 0.0;
                d->h->bias[z_index][M_index] = 0.0;
            }
            else
            {
                SAFEHMPDF(dndlogM(d, z_index, M_index,
                                  d->h->hmf[z_index]+M_index,
                                  d->h->bias[z_index]+M_index));

                if (d->h->bias_resc != NULL)
                {
                    d->h->bias[z_index][M_index] *= d->h->bias_resc(d->n->zgrid[z_index],
                                                                    d->n->Mgrid[M_index] * d->c->h,
                                                                    d->h->bias_resc_params);
                }
            }
        }
    }

    ENDFCT
}//}}}

int
init_halo_model(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->h->inited_halo) { return 0; }

    HMPDFPRINT(1, "init_halo_model\n");

    SAFEHMPDF(create_c_of_y(d));
    SAFEHMPDF(create_dndlogM(d));
    d->h->inited_halo = 1;

    ENDFCT
}//}}}

