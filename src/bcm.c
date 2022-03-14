#include <stdio.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

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

#ifdef ARICO20
#   pragma message( "Using the updated Arico BCM from 2009.14225" )
#else
#   pragma message( "Using the original Arico BCM from 1911.08471" )
#endif

#ifndef M_SQRT5
#   define M_SQRT5 2.2360679774997896964091736687312762
#endif

int
null_bcm(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->bcm->inited_bcm = 0;
    d->bcm->radii = NULL;
    d->bcm->ws = NULL;
    d->bcm->profiles_indices = NULL;

    ENDFCT
}//}}}

int
reset_bcm(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\treset_bcm\n");

    if (d->bcm->radii != NULL) { free(d->bcm->radii); }
    if (d->bcm->ws != NULL)
    {
        for (int ii=0; ii<d->Ncores; ii++)
        {
            bcm_delete_ws(d->bcm->ws[ii]);
            free(d->bcm->ws[ii]);
        }
        free(d->bcm->ws);
    }
    if (d->bcm->profiles_indices != NULL) { free(d->bcm->profiles_indices); }

    ENDFCT
}//}}}

int
init_bcm(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->bcm->inited_bcm) { return 0; }

    HMPDFPRINT(1, "init_bcm\n");

    // TODO check that these settings make sense
    d->bcm->Nradii = 10000;
    d->bcm->rmin = 1e-5;
    d->bcm->rmax = 2e1;
    SAFEALLOC(d->bcm->radii, malloc(d->bcm->Nradii * sizeof(double)));
    SAFEHMPDF(logspace(d->bcm->Nradii, d->bcm->rmin, d->bcm->rmax, d->bcm->radii));

    // find the radius closest(ish) to R200c -- this is not super important and only
    // a small optimization in the computation of xi
    for (d->bcm->R200c_idx=0; d->bcm->radii[d->bcm->R200c_idx]<1.0; d->bcm->R200c_idx++);
    
    // allocate the workspaces
    SAFEALLOC(d->bcm->ws, malloc(d->Ncores * sizeof(bcm_ws *)));
    for (int ii=0; ii<d->Ncores; ii++)
    {
        SAFEALLOC(d->bcm->ws[ii], malloc(sizeof(bcm_ws)));
        SAFEHMPDF(bcm_new_ws(d, d->bcm->ws[ii]));
    }

    // if user requested, compute the profiles indices
    if (d->bcm->profiles_N)
    {
        HMPDFCHECK(!d->bcm->profiles_fnames, "profiles_fnames not provided");
        HMPDFCHECK(!d->bcm->profiles_where, "profiles_where not provided");
        HMPDFCHECK(!d->bcm->profiles_r, "profiles_r not provided");
        HMPDFCHECK(!d->bcm->profiles_Nr, "profiles_Nr not provided");

        SAFEALLOC(d->bcm->profiles_indices, malloc(2 * d->bcm->profiles_N * sizeof(int)));

        for (int ii=0; ii<d->bcm->profiles_N; ii++)
        {
            double ztarg = d->bcm->profiles_where[2*ii];
            double Mtarg = d->bcm->profiles_where[2*ii+1] / d->c->h;

            d->bcm->profiles_indices[2*ii] = find_closest(d->n->Nz, d->n->zgrid, ztarg);
            d->bcm->profiles_indices[2*ii+1] = find_closest(d->n->NM, d->n->Mgrid, Mtarg);
        }
    }

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

    SAFEALLOC(ws->bg_integr_ws, gsl_integration_workspace_alloc(BCM_BGINTEGR_LIMIT));

    ENDFCT
}//}}}

int
bcm_delete_ws(bcm_ws *ws)
{//{{{
    STARTFCT

    free(ws->dm_xi);
    gsl_interp_accel_free(ws->dm_r_accel);
    gsl_interp_free(ws->dm_xi_interp);
    gsl_integration_workspace_free(ws->bg_integr_ws);

    ENDFCT
}//}}}

inline static double
interpolate_param(hmpdf_obj *d, int idx, double ztarg, int logarithmic)
// do the redshift interpolation of the Arico20 parameters
{//{{{
    int Nz = d->bcm->Arico20_Nz;
    double *z = d->bcm->Arico20_z;
    double *params = d->bcm->Arico20_params;

    // trivial cases
    if (Nz==1 || ztarg<z[0])
        return params[idx];

    if (ztarg > z[Nz-1])
        return params[(Nz-1)*hmpdf_Arico20_Nparams + idx];

    double zlo, zhi, plo, phi;
    int found = 0;
    for (int ii=0; ii<Nz-1; ii++)
        if (ztarg > z[ii] && ztarg < z[ii+1])
        {
            zlo = z[ii];
            zhi = z[ii+1]; 
            plo = params[ii*hmpdf_Arico20_Nparams + idx];
            phi = params[(ii+1)*hmpdf_Arico20_Nparams + idx];
            found = 1;
            break;
        }

    if (!found)
        return NAN;
    
    if (logarithmic)
    {
        plo = log(plo);
        phi = log(phi);
    }

    double out = plo + (phi - plo) * (ztarg - zlo) / (zhi - zlo);

    if (logarithmic)
        out = exp(out);

    return out;
}//}}}

#ifndef ARICO20
inline static double
rho_bg_g0(double x, bcm_ws *ws)
{//{{{
    return pow(log1p(x)/x, ws->bg_Gamma);
}//}}}

inline static double
rho_bg_g1(double x, bcm_ws *ws)
{//{{{
    return 1.0/(x * gsl_pow_2(1.0+x));
}//}}}

inline static double
rho_bg_x2g0(double x, void *p)
{//{{{
    if (x < 1e-8)
        return 0.0; // prevent spurious divergence

    double Gamma = *(double *)p;

    return pow(log1p(x), Gamma) * pow(x, 2.0-Gamma);
}//}}}

static int
M_bg0(double r, bcm_ws *ws, double other_mass, double *out)
// return the mass enclosed in the inner part
{//{{{
    STARTFCT

    double rout = GSL_MIN(r, ws->R200c/M_SQRT5);

    gsl_function integrand;

    integrand.function = &rho_bg_x2g0;
    integrand.params = &(ws->bg_Gamma);

    double scaling = ws->bg_y0 * gsl_pow_3(ws->rs);

    double err;
    SAFEGSL(gsl_integration_qag(&integrand, 0.0, rout/ws->rs,
                                (other_mass > 0.0) ? BCM_BGINTEGR_EPSABS * other_mass / scaling : 0.0,
                                BCM_BGINTEGR_EPSREL,
                                BCM_BGINTEGR_LIMIT, BCM_BGINTEGR_KEY,
                                ws->bg_integr_ws, out, &err));

    // normalize properly
    *out *= scaling;

    ENDFCT
}//}}}

inline static double
M_bg1(double r, bcm_ws *ws)
// return mass enclosed in outer part
{//{{{
    const double rmin = ws->R200c/M_SQRT5;

    if (r < rmin)
        return 0.0;
    
    return ws->bg_y1 * 4.0 * M_PI * gsl_pow_3(ws->rs)
           * ( rmin/(rmin+ws->rs) - r/(r+ws->rs)
              +log((r+ws->rs)/(rmin+ws->rs)) );
}//}}}
#endif // ifndef ARICO20

inline static double 
rho_bg(double r, bcm_ws *ws)
// bound gas density profile
{//{{{
    if (r > ws->R200c)
        return 0.0;

#ifdef ARICO20
    return ws->bg_y0 * pow(1.0+r/ws->bg_r_inn, -ws->bg_beta_i)
           / gsl_pow_2( 1.0+gsl_pow_2(r/ws->bg_r_out) );
#else
    double x = r / ws->rs;
    if (r < ws->R200c/M_SQRT5)
        return ws->bg_y0 * rho_bg_g0(x, ws);
    else
        return ws->bg_y1 * rho_bg_g1(x, ws);
#endif
}//}}}

#ifdef ARICO20
inline static double
M_bg_integrand(double r, void *p)
// basically the same function as above only with 4pi r^2 Jacobian
// and GSL-compatible
{//{{{
    bcm_ws *ws = (bcm_ws *)p;
    return 4.0 * M_PI * gsl_pow_2(r) * rho_bg(r, ws);
}//}}}
#endif // ifdef Arico20

static int
M_bg(double r, bcm_ws *ws, double other_mass, double *out)
// bound gas mass profile
// The other_mass is the remaining baryonic mass, which should be computed beforehand.
// This can be used to robustly terminate the integration such that the total baryonic
// mass is well approximated.
// Can also pass a negative number for this parameter in case it is not known.
{//{{{
    STARTFCT

#ifdef ARICO20
    gsl_function integrand;

    integrand.function = &M_bg_integrand;
    integrand.params = ws;

    double err;
    SAFEGSL(gsl_integration_qag(&integrand, 0.0, r,
                                (other_mass > 0.0) ? BCM_BGINTEGR_EPSABS * other_mass : 0.0,
                                BCM_BGINTEGR_EPSREL,
                                BCM_BGINTEGR_LIMIT, BCM_BGINTEGR_KEY,
                                ws->bg_integr_ws, out, &err));
#else 
    SAFEHMPDF(M_bg0(r, ws, other_mass, out));
    *out += M_bg1(r, ws);
#endif

    ENDFCT
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
    const double pref = 4.0 * POW3(M_SQRTPI);
    #undef POW3

    return ws->cg_y0 * pref * erf(0.5*r/ws->cg_Rh);
}//}}}

#ifdef ARICO20
inline static double
rho_rg(double r, bcm_ws *ws)
// re-accreted gas density profile
{//{{{
    const double pref = 0.5 * M_SQRT1_2 * M_2_SQRTPI;

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
               *(erf(0.5*(rout-ws->rg_mu)/ws->rg_sigma)+erf(0.5*ws->rg_mu/ws->rg_sigma)) );
}//}}}
#endif // ifdef ARICO20

inline static double
rho_eg(double r, bcm_ws *ws)
// ejected gas density profile
{//{{{
    #define POW3(x) ((x)*(x)*(x))
    const double pref = POW3(0.5 * M_SQRT1_2 * M_2_SQRTPI); // 1/(2pi)^3/2
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
           * (  erf(M_SQRT1_2*x)
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
    *out = ws->dm_f * rho_nfw(r / xi, ws) * (xi - r/ws->R200c*xiprime) / gsl_pow_4(xi);

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

    double log10M1z0cen = interpolate_param(d, hmpdf_Arico20_M_1_z0_cen, z, 0);

    double log10_M1 = log10M1z0cen + nu * (-1.793*(a-1.0)-0.251*z);
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

inline static int
find_xi_at_rf(double rf, double xi_init, bcm_ws *ws, double *out)
// xi_init can be a guess for what xi should be, will be used to speed
//         up the search if chosen well.
{//{{{
    STARTFCT

#ifdef ARICO20
    double M_bary = M_cg(rf, ws) + M_rg(rf, ws) + M_eg(rf, ws);
#else
    double M_bary = M_cg(rf, ws) + M_eg(rf, ws);
#endif
    double m_bg;
    SAFEHMPDF(M_bg(rf, ws, M_bary, &m_bg));
    M_bary += m_bg;

    double xi, diff;

    do
    {
        double Mi = M_nfw(rf/xi_init, ws);
        double Mf = ws->dm_f * Mi + M_bary;

        xi = 1.0 + 0.3 * ( gsl_pow_2(Mi/Mf) - 1.0 );
        diff = fabs(xi - xi_init);
        xi_init = xi;
    } while (diff > BCM_XI_SEARCH_TOL);

    *out = xi;

    ENDFCT
}//}}}

static int
interpolate_xi(hmpdf_obj *d, bcm_ws *ws)
{//{{{
    STARTFCT

    // we know that at R200c xi=1, so this is a good place to start
    int start_idx = d->bcm->R200c_idx;

    SAFEHMPDF(find_xi_at_rf(d->bcm->radii[start_idx]*ws->R200c,
                            1.0, ws, ws->dm_xi+start_idx));

    // now go upwards in r
    for (int r_index=start_idx+1; r_index<d->bcm->Nradii; r_index++)
        SAFEHMPDF(find_xi_at_rf(d->bcm->radii[r_index]*ws->R200c,
                                ws->dm_xi[r_index-1], ws, ws->dm_xi+r_index));

    // now go downwards in r
    for (int r_index=start_idx-1; r_index>=0; r_index--)
        SAFEHMPDF(find_xi_at_rf(d->bcm->radii[r_index]*ws->R200c,
                                ws->dm_xi[r_index+1], ws, ws->dm_xi+r_index));

    // TODO we may need a smoothing function here to get rid of small-scale noise
    //      --> not sure if this is actually true, results look quite ok!
    SAFEGSL(gsl_interp_init(ws->dm_xi_interp, d->bcm->radii, ws->dm_xi, d->bcm->Nradii));

    ENDFCT
}//}}}

int
bcm_init_ws(hmpdf_obj *d, int z_index, int M_index, double mass_resc, bcm_ws *ws)
{//{{{
    STARTFCT

    // first compute basic NFW properties
    double c200c;
    SAFEHMPDF(Mconv(d, z_index, M_index, hmpdf_mdef_c, mass_resc, &(ws->M200c), &(ws->R200c), &c200c));
    SAFEHMPDF(NFW_fundamental(d, z_index, M_index, mass_resc, &(ws->rhos), &(ws->rs)));

    // compute Arico parameters at our redshift
    double z = d->n->zgrid[z_index];
    double M_c_this_z = interpolate_param(d, hmpdf_Arico20_M_c, z, 1);
    double beta_this_z = interpolate_param(d, hmpdf_Arico20_beta, z, 0);
    double eta_this_z = interpolate_param(d, hmpdf_Arico20_eta, z, 0);

#ifdef ARICO20
    double M_r_this_z = interpolate_param(d, hmpdf_Arico20_M_r, z, 1);
    double theta_inn_this_z = interpolate_param(d, hmpdf_Arico20_theta_inn, z, 0);
    double theta_out_this_z = interpolate_param(d, hmpdf_Arico20_theta_out, z, 0);
    double M_inn_this_z = interpolate_param(d, hmpdf_Arico20_M_inn, z, 1);
#endif

    // compute the mass fractions
#ifdef ARICO20
    double f_dm, f_cg, f_sg, f_hg, f_rg, f_bg, f_eg;

    f_dm = 1.0 - d->c->Ob_0 / d->c->Om_0;

    f_cg = f_cg_fit(d, z_index, ws->M200c);
    f_sg = f_cg; // according to Arico, this is OK and shouldn't matter much -- TODO test variations here

    f_hg = ( d->c->Ob_0/d->c->Om_0 - f_cg - f_sg )
           / ( 1.0 + pow(M_c_this_z/ws->M200c, beta_this_z) );

    static const double beta_r = 2.0;
    f_rg = f_hg * pow(M_c_this_z/ws->M200c, beta_this_z)
                / ( 1.0 + pow(M_r_this_z/ws->M200c, beta_r) );

    f_bg = f_hg - f_rg;

    f_eg = d->c->Ob_0/d->c->Om_0 - f_cg - f_sg - f_hg;
#else
    double f_dm, f_cg, f_bg, f_eg;

    f_dm = 1.0 - d->c->Ob_0/d->c->Om_0;

    f_cg = f_cg_fit(d, z_index, ws->M200c);

    f_bg = (d->c->Ob_0/d->c->Om_0 - f_cg)
           / (1.0 + pow(M_c_this_z/ws->M200c, beta_this_z) );

    f_eg = d->c->Ob_0/d->c->Om_0 - f_cg - f_bg;
#endif

    ws->dm_f = f_dm;
    ws->eg_f = f_eg;

    // populate parameters and compute normalization factors

#ifdef ARICO20
    ws->bg_y0 = 1.0;
    ws->bg_r_inn = theta_inn_this_z * ws->R200c;
    ws->bg_r_out = theta_out_this_z * ws->R200c;
    ws->bg_beta_i = 3.0 - pow(M_inn_this_z/ws->M200c, 0.31);
    double m_bg;
    SAFEHMPDF(M_bg(ws->R200c, ws, -1, &m_bg));
    ws->bg_y0 = f_bg * ws->M200c / m_bg;
#else
    ws->bg_y0 = 1.0;
    ws->bg_y1 = 1.0;
    double csqrt5 = c200c / M_SQRT5;
    ws->bg_Gamma = (1.0+3.0*csqrt5) * log1p(csqrt5)
                   / ( (1.0+csqrt5)*log1p(csqrt5) - csqrt5 );
    double m_bg0, m_bg1;
    SAFEHMPDF(M_bg0(ws->R200c/M_SQRT5, ws, -1, &m_bg0));
    m_bg1 = M_bg1(ws->R200c, ws);
    double g0, g1;
    g0 = rho_bg_g0(ws->R200c/ws->rs/M_SQRT5, ws);
    g1 = rho_bg_g1(ws->R200c/ws->rs/M_SQRT5, ws);
    // fix y0, y1 by requiring continuity and correct integral
    ws->bg_y0 = f_bg * ws->M200c * g1 / (g1*m_bg0 + g0*m_bg1);
    ws->bg_y1 = f_bg * ws->M200c * g0 / (g1*m_bg0 + g0*m_bg1);
#endif

    ws->cg_y0 = 1.0;
    ws->cg_Rh = 0.015 * ws->R200c;
    ws->cg_y0 = f_cg * ws->M200c / M_cg(ws->R200c, ws);

#ifdef ARICO20
    ws->rg_y0 = 1.0;
    ws->rg_mu = 0.3 * ws->R200c;
    ws->rg_sigma = 0.1 * ws->R200c;
    ws->rg_y0 = f_rg * ws->M200c / M_rg(ws->R200c, ws);
#endif

    double resc = 0.5 * 10.0 * M_SQRT2 * ws->R200c;
    ws->eg_rej = eta_this_z * 0.75 * resc;

    // compute the dark matter relaxation
    SAFEHMPDF(interpolate_xi(d, ws));

    // if requested, save corresponding profiles to file
    if (d->bcm->profiles_indices)
    {
        int do_it = -1;
        for (int ii=0; ii<d->bcm->profiles_N; ii++)
            if (   z_index == d->bcm->profiles_indices[2*ii]
                && M_index == d->bcm->profiles_indices[2*ii+1])
            {
                do_it = ii;
                break;
            }

        if (do_it >= 0)
            SAFEHMPDF(bcm_profiles_to_file(d, z_index, M_index, ws, d->bcm->profiles_fnames[do_it]));
    }

    ENDFCT
}//}}}

int
bcm_density_profile(hmpdf_obj *d, bcm_ws *ws, double r, double *out)
{//{{{
    STARTFCT

#ifdef ARICO20
    double rho_bary = rho_bg(r, ws) + rho_cg(r, ws) + rho_rg(r, ws) + rho_eg(r, ws);
#else
    double rho_bary = rho_bg(r, ws) + rho_cg(r, ws) + rho_eg(r, ws);
#endif
    double rho_rdm;
    SAFEHMPDF(rho_dm(d, r, ws, &rho_rdm));
    *out = rho_bary + rho_rdm;

    // do the DM particles that have not been translated (outside R200c)
    if (r > ws->R200c)
    {
        double x = r / ws->rs;
        *out += ws->dm_f * ws->rhos / (x * gsl_pow_2(1.0 + x));
    }

    ENDFCT
}//}}}

int
bcm_profiles_to_file(hmpdf_obj *d, int z_index, int M_index, bcm_ws *ws, char *fname)
{//{{{
    STARTFCT

    FILE *fp = fopen(fname, "w");
    HMPDFCHECK(!fp, "failed to open file %s", fname);

    // print information about the object
    fprintf(fp, "# z=%.18e\n# M200m[Msun/h]=%.18e\n# R200c[Mpc]=%.18e\n\n",
                d->n->zgrid[z_index], d->n->Mgrid[M_index]*d->c->h, ws->R200c); 

    // print header
#ifdef ARICO20
    fprintf(fp, "# r[Mpc], rho_bg, rho_cg, rho_rg, rho_eg, rho_dm [Msun/Mpc^3]\n");
#else
    fprintf(fp, "# r[Mpc], rho_bg, rho_cg, rho_eg, rho_dm [Msun/Mpc^3]\n");
#endif

    // now print the profiles
    for (int ii=0; ii<d->bcm->profiles_Nr; ii++)
    {
        double r = ws->R200c * d->bcm->profiles_r[ii];
        double bg = rho_bg(r, ws);
        double cg = rho_cg(r, ws);
#ifdef ARICO20
        double rg = rho_rg(r, ws);
#endif
        double eg = rho_eg(r, ws);
        double dm;
        SAFEHMPDF(rho_dm(d, r, ws, &dm));

        // add the undisturbed NFW component
        if (r > ws->R200c)
        {
            double x = r / ws->rs;
            dm += ws->dm_f * ws->rhos / (x * gsl_pow_2(1.0 + x));
        }

#ifdef ARICO20
        fprintf(fp, "%.18e %.18e %.18e %.18e %.18e %.18e\n", r, bg, cg, rg, eg, dm);
#else
        fprintf(fp, "%.18e %.18e %.18e %.18e %.18e\n", r, bg, cg, eg, dm);
#endif
    }

    fclose(fp);

    ENDFCT
}//}}}
