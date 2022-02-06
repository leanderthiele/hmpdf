#include <stdio.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>

#include "configs.h"
#include "utils.h"
#include "object.h"
#include "cosmology.h"
#include "halo_model.h"

#include "hmpdf.h"

/* BCM components:
 * 
 * -- bg: bound gas
 * -- cg: central galaxy
 * -- rg: re-accreted gas
 * -- eg: ejected gas
 * -- dm: (relaxed) dark matter
 *
 */

int
null_bcm(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->bcm->inited_bcm = 0;
    d->bcm->radii = NULL;

    ENDFCT
}//}}}

int
reset_bcm(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\treset_bcm\n");

    if (d->bcm->radii != NULL) { free(d->bcm->radii); }

    ENDFCT
}//}}}

int
init_bcm(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->bcm->inited_bcm) { return 0; }

    HMPDFPRINT(1, "init_bcm\n");

    // TODO allocate and fill the array with radii,
    //      perhaps other stuff necessary

    d->bcm->inited_bcm = 1;

    ENDFCT
}//}}}

int
bcm_new_ws(hmpdf_obj *d, bcm_ws *ws)
{//{{{
    STARTFCT

    SAFEALLOC(ws->dm_xi, malloc(d->bcm->Nradii * sizeof(double)));

    SAFEALLOC(ws->dm_r_accel, gsl_interp_accel_alloc());
    // TODO think about interpolation type here, will need to see examples
    SAFEALLOC(ws->dm_xi_interp, gsl_interp_alloc(gsl_interp_linear, d->bcm->Nradii));

    ENDFCT
}//}}}

int
bcm_delete_ws(bcm_ws *ws)
{//{{{
    STARTFCT

    free(ws->dm_xi);
    gsl_interp_accel_free(ws->dm_r_accel);
    gsl_interp_free(ws->dm_xi_interp);

    ENDFCT
}//}}}

inline static double 
rho_bg(double r, bcm_ws *ws)
// bound gas density profile
{//{{{
    if (r > ws->R200c)
        return 0.0;

    return ws->bg_y0 * pow(1.0+r/ws->bg_r_inn, -ws->bg_beta_i)
           / gsl_pow_2( 1.0+gsl_pow_2(r/ws->bg_r_out) );
}//}}}

inline static double
M_bg(double r, bcm_ws *ws)
// bound gas mass profile
{//{{{
    // TODO there is a closed form expression for this but it has special
    //      functions of complex arguments...

    return 0.0;
}//}}}

inline static double
rho_cg(double r, bcm_ws *ws)
// central galaxy density profile
{//{{{
    return ws->cg_y0 / (ws->cg_Rh * gsl_pow_2(r))
           * exp( -gsl_pow_2(0.5 * r / ws->cg_Rh) );
}//}}}

inline static double
M_cg(double r, bcm_ws *ws)
// central galaxy mass profile
{//{{{
    // y0 * 4 pi^(3/2) * Erf(r/2Rh)
    #define POW3(x) ((x)*(x)*(x))
    static double pref = 4.0 * POW3(M_SQRTPI);
    #undef POW3

    return ws->cg_y0 * pref * gsl_sf_erf(0.5*r/ws->cg_Rh);
}//}}}

inline static double
rho_rg(double r, bcm_ws *ws)
// re-accreted gas density profile
{//{{{
    static const double pref = 0.5 * M_SQRT1_2 * M_2_SQRTPI;

    if (r > ws->R200c)
        return 0.0;
    
    return ws->rg_y0 * pref / ws->rg_sigma
           * exp( -gsl_pow_2(0.5 * (r-ws->rg_mu) / ws->rg_sigma) );
}//}}}

inline static double
M_rg(double r, bcm_ws *ws)
// re-accreted gas mass profile
{//{{{
    double rout = GSL_MIN(r, ws->R200c);

    return ws->rg_y0 * 2.0 * M_SQRT2 * M_SQRTPI
           * ( 2.0*ws->rg_mu*ws->rg_sigma*exp(-gsl_pow_2(0.5*ws->rg_mu/ws->rg_sigma))
              -2.0*(rout+ws->rg_mu)*ws->rg_sigma*exp(-gsl_pow_2(0.5*(rout-ws->rg_mu)/ws->rg_sigma))
              +M_SQRTPI*(gsl_pow_2(ws->rg_mu)+2.0*gsl_pow_2(ws->rg_sigma))
               *(gsl_sf_erf(0.5*(rout-ws->rg_mu)/ws->rg_sigma)+gsl_sf_erf(0.5*ws->rg_mu/ws->rg_sigma)) );
}//}}}

inline static double
rho_eg(double r, bcm_ws *ws)
// ejected gas density profile
{//{{{
    #define POW3(x) ((x)*(x)*(x))
    static const double pref = POW3(0.5 * M_SQRT1_2 * M_2_SQRTPI);
    #undef POW3

    return ws->eg_f * ws->M200c * pref / gsl_pow_3(ws->eg_rej)
           * exp( -0.5*gsl_pow_2(r/ws->eg_rej) );
}//}}}

inline static double
M_eg(double r, bcm_ws *ws)
// ejected gas mass profile
{//{{{
    double x = r / ws->eg_rej;
    return ws->eg_f * ws->M200c
           * (  gsl_sf_erf(M_SQRT1_2*x)
              - M_SQRT1_2*M_2_SQRTPI*x*exp(-0.5*gsl_pow_2(x)) );
}//}}}

inline static double
rho_nfw(double r, bcm_ws *ws)
// original dark matter NFW density profile
{//{{{
    if (r > ws->R200c)
        return 0.0;

    double x = r / ws->rs;
    return ws->rhos / (x * gsl_pow_2(1.0 + x));
}//}}}

inline static double
M_nfw(double r, bcm_ws *ws)
// original dark matter NFW mass profile
{//{{{
    double rout = GSL_MIN(r, ws->R200c);

    return 4.0 * M_PI * ws->rhos * gsl_pow_3(ws->rs)
           * ( log1p(rout/ws->rs) - rout/(rout+ws->rs) );
}//}}}

inline static int 
rho_dm(hmpdf_obj *d, double r, bcm_ws *ws, double *out)
{//{{{
    STARTFCT

    double xi, xiprime;
    SAFEGSL(gsl_interp_eval_e(ws->dm_xi_interp, d->bcm->radii, ws->dm_xi,
                              r/ws->R200c, ws->dm_r_accel, &xi));
    SAFEGSL(gsl_interp_eval_deriv_e(ws->dm_xi_interp, d->bcm->radii, ws->dm_xi,
                                    r/ws->R200c, ws->dm_r_accel, &xiprime));

    // TODO check this expression
    *out = ws->dm_f * rho_nfw(r, ws) * (xi - r/ws->R200c*xiprime) / gsl_pow_4(xi);

    ENDFCT
}//}}}

inline static int
M_dm(hmpdf_obj *d, double r, bcm_ws *ws, double *out)
{//{{{
    STARTFCT

    double xi;
    SAFEGSL(gsl_interp_eval_e(ws->dm_xi_interp, d->bcm->radii, ws->dm_xi,
                              r/ws->R200c, ws->dm_r_accel, &xi));

    *out = ws->dm_f * M_nfw(r / xi, ws);

    ENDFCT
}//}}}

inline static double
f_cg_fit(hmpdf_obj *d, int z_index, double M200c)
{//{{{
    double z = d->n->zgrid[z_index];
    double a = 1.0 / (1.0+z);
    double nu = exp( -gsl_pow_2(2.0*a) );
    double log10_M1 = d->bcm->Arico20_params[hmpdf_Arico20_M_1_z0_cen]
                      + nu * (-1.793*(a-1.0)-0.251*z);
    double log10_eps = -1.6382721639824072 + nu*-0.006*(a-1.0) + -0.119*(a-1.0);
    double alpha = -1.779 + nu*0.731*(a-1.0);
    double delta = 4.394 + nu*(2.608*(a-1.0) + -0.043*z);
    double gamma = 0.547 + nu*(1.319*(a-1.0) +  0.279*z);
    double log10_M200_M1 = log10(M200c) - log10_M1;
    double gx = - log1p(exp(M_LN10*alpha*log10_M200_M1))/M_LN10
                + delta*pow(log1p(exp(log10_M200_M1))/M_LN10, gamma)
                  / ( 1.0 + exp(exp(-M_LN10*log10_M200_M1)) );
    double g0 = -M_LN2/M_LN10 + delta * pow(M_LN2/M_LN10, gamma) / (1.0+M_E);
    double log10_f_cg = log10_eps - log10_M200_M1 + gx - g0;
    return exp(M_LN10 * log10_f_cg);
}//}}}

inline static double
find_xi_at_rf(double rf, double xi_init, bcm_ws *ws)
// xi_init can be a guess for what xi should be, will be used to speed
//         up the search if chosen well.
{//{{{
    double M_bary = M_bg(rf, ws) + M_cg(rf, ws) + M_rg(rf, ws) + M_eg(rf, ws);

    double xi, diff;

    do
    {
        double Mi = M_nfw(rf/xi_init, ws);
        double Mf = ws->dm_f * Mi + M_bary;
        xi = 1.0 + 0.3 * ( gsl_pow_2(Mi/Mf) - 1.0 );
        diff = fabs(xi - xi_init);
        xi_init = xi;
    } while (diff > XI_SEARCH_TOL);

    return xi;
}//}}}

static int
interpolate_xi(hmpdf_obj *d, bcm_ws *ws)
{//{{{
    STARTFCT

    for (int r_index=0; r_index<d->bcm->Nradii; r_index++)
        ws->dm_xi[r_index] = find_xi_at_rf(d->bcm->radii[r_index]*ws->R200c,
                                           (r_index) ? ws->dm_xi[r_index-1] : 1.0,
                                           ws);

    SAFEGSL(gsl_interp_init(ws->dm_xi_interp, d->bcm->radii, ws->dm_xi, d->bcm->Nradii));

    ENDFCT
}//}}}

int
bcm_init_ws(hmpdf_obj *d, int z_index, int M_index, bcm_ws *ws)
{//{{{
    STARTFCT

    // magic numbers that were set to fixed values
    static const double beta_r = 2.0;

    // first compute basic NFW properties
    double _c;
    SAFEHMPDF(Mconv(d, z_index, M_index, hmpdf_mdef_c, 1.0, &(ws->M200c), &(ws->R200c), &_c));
    SAFEHMPDF(NFW_fundamental(d, z_index, M_index, 1.0, &(ws->rhos), &(ws->rs)));

    // compute the mass fractions
    double f_dm, f_cg, f_sg, f_hg, f_rg, f_bg, f_eg;

    f_dm = 1.0 - d->c->Ob_0 / d->c->Om_0;

    f_cg = f_cg_fit(d, z_index, ws->M200c);
    f_sg = f_cg; // according to Acori, this is OK and shouldn't matter much -- TODO test variations here

    f_hg = ( d->c->Ob_0/d->c->Om_0 - f_cg - f_sg )
           / ( 1.0 + pow(d->bcm->Arico20_params[hmpdf_Arico20_M_c]/ws->M200c,
                         d->bcm->Arico20_params[hmpdf_Arico20_beta]) );

    f_rg = f_hg * pow(d->bcm->Arico20_params[hmpdf_Arico20_M_c]/ws->M200c,
                      d->bcm->Arico20_params[hmpdf_Arico20_beta])
                / ( 1.0 + pow(d->bcm->Arico20_params[hmpdf_Arico20_M_r]/ws->M200c, beta_r) );

    f_bg = f_hg - f_rg;

    f_eg = d->c->Ob_0/d->c->Om_0 - f_cg - f_sg - f_hg;

    ws->dm_f = f_dm;
    ws->eg_f = f_eg;

    // populate parameters and compute normalization factors

    ws->bg_y0 = 1.0;
    ws->bg_r_inn = d->bcm->Arico20_params[hmpdf_Arico20_theta_inn]*ws->R200c;
    ws->bg_r_out = d->bcm->Arico20_params[hmpdf_Arico20_theta_out]*ws->R200c;
    ws->bg_beta_i = 3.0 - pow(d->bcm->Arico20_params[hmpdf_Arico20_M_inn]/ws->M200c, 0.31);
    ws->bg_y0 = f_bg * ws->M200c / M_bg(ws->R200c, ws);

    ws->cg_y0 = 1.0;
    ws->cg_Rh = 0.015 * ws->R200c;
    ws->cg_y0 = f_cg * ws->M200c / M_cg(ws->R200c, ws);

    ws->rg_y0 = 1.0;
    ws->rg_mu = 0.3 * ws->R200c;
    ws->rg_sigma = 0.1 * ws->R200c;
    ws->rg_y0 = f_rg * ws->M200c / M_rg(ws->R200c, ws);

    double resc = 0.5 * 10.0 * M_SQRT2 * ws->R200c;
    ws->eg_rej = d->bcm->Arico20_params[hmpdf_Arico20_eta] * 0.75 * resc;

    // compute the dark matter relaxation
    
    SAFEHMPDF(interpolate_xi(d, ws));

    ENDFCT
}//}}}

int
bcm_density_profile(hmpdf_obj *d, bcm_ws *ws, double r, double *out)
{//{{{
    STARTFCT

    double rho_bary = rho_bg(r, ws) + rho_cg(r, ws) + rho_rg(r, ws) + rho_eg(r, ws);
    double rho_rdm;
    SAFEHMPDF(rho_dm(d, r, ws, &rho_rdm));
    *out = rho_bary + rho_rdm;

    ENDFCT
}//}}}
