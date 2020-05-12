#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "utils.h"
#include "configs.h"
#include "data.h"
#include "halo_model.h"

static
double DeltaVir_BryanNorman98(all_data *d, int z_index)
{//{{{
    double x = d->c->Om[z_index] - 1.0;
    return 18.0*M_PI*M_PI + 82.0*x - 39.0*x*x;
}//}}}

static
double density_threshold(all_data *d, int z_index, mdef mdef)
{//{{{{
    switch (mdef)
    {
        case mdef_c : return 200.0 * d->c->rho_c[z_index];
        case mdef_v : return DeltaVir_BryanNorman98(d, z_index) * d->c->rho_c[z_index];
        case mdef_m : return 200.0 * d->c->rho_m[z_index];
        default     : return 0.0;
    }
}//}}}

static
double RofM(all_data *d, int z_index, int M_index)
{//{{{
    return cbrt(3.0*d->n->gr->Mgrid[M_index]/4.0/M_PI
                /density_threshold(d, z_index, MDEF_GLOBAL));
}//}}}

static
double MofR(all_data *d, int z_index, double R, mdef mdef)
{//{{{
    return 4.0*M_PI*density_threshold(d, z_index, mdef)*gsl_pow_3(R)/3.0;
}//}}}

                                  //   A       B      C
double Duffy08_fit_params[][3] = {{ 5.71, -0.087, -0.47},  // M200c//{{{
                                  { 7.85, -0.081, -0.71},  // Mvir
                                  {10.14, -0.081, -1.01}}; // M200m//}}}
static
double _c_Duffy08(all_data *d, double z, double M, mdef mdef)
{//{{{
    // convert to Msun/h
    M *= d->c->h;
    return Duffy08_fit_params[mdef][0]
           * pow(M/2.0e12, Duffy08_fit_params[mdef][1])
           * pow(1.0 + z,  Duffy08_fit_params[mdef][2]);
}//}}}
static
double c_Duffy08(all_data *d, int z_index, int M_index)
{//{{{
    return _c_Duffy08(d, d->n->gr->zgrid[z_index],
                      d->n->gr->Mgrid[M_index],
                      MDEF_GLOBAL);
}//}}}

double NFW_fundamental(all_data *d, int z_index, int M_index, double *rs)
// returns rhos via function call and rs via return value
// this function is tested against Colossus --> everything here works
{//{{{
    double c = c_Duffy08(d, z_index, M_index);
    *rs = RofM(d, z_index, M_index) / c;
    return d->n->gr->Mgrid[M_index]/4.0/M_PI/gsl_pow_3(*rs)
           / (log1p(c)-c/(1.0+c));
}//}}}

static
void create_c_of_y(all_data *d)
{//{{{
    printf("\tcreate_c_of_y\n");
    double *logc_grid = (double *)malloc(CINTERP_NC * sizeof(double));
    double *logy_grid = (double *)malloc(CINTERP_NC * sizeof(double));
    for (int ii=0; ii<CINTERP_NC; ii++)
    {
        logc_grid[ii] = log(CINTERP_CMAX)
                        - (double)(ii)*log(CINTERP_CMAX/CINTERP_CMIN)/(double)(CINTERP_NC-1);
        double _c = exp(logc_grid[ii]);
        logy_grid[ii] = log(3.0) - 3.0*logc_grid[ii] + log(log1p(_c) - _c/(1.0+_c));
    }
    d->h->c_interp = gsl_spline_alloc(gsl_interp_cspline, CINTERP_NC);
    d->h->c_accel = gsl_interp_accel_alloc();
    gsl_spline_init(d->h->c_interp, logy_grid, logc_grid, CINTERP_NC);
    free(logc_grid);
    free(logy_grid);
}//}}}

static
double c_of_y(all_data *d, double y)
// inverts the function y(c) = 3/c^3 * (log(1+c)-c/(1+c))
// this function is much smoother on loglog scale, so do the interpolation this way
{//{{{
    return exp(gsl_spline_eval(d->h->c_interp, log(y), d->h->c_accel));
}//}}}

double Mconv(all_data *d, int z_index, int M_index, mdef mdef_out, double *R, double *c)
// returns the converted mass via function call and the new radius and concentration via return value
// this function is tested against Colossus --> everything here works
{//{{{
    double rhos, rs;
    rhos = NFW_fundamental(d, z_index, M_index, &rs);
    *c = c_of_y(d, density_threshold(d, z_index, mdef_out)/rhos);
    *R = rs * *c;
    return MofR(d, z_index, *R, mdef_out);
}//}}}

static
double fnu_Tinker10(double nu, double z)
{//{{{
    z = (z<3.0) ? z : 3.0;
    double beta  =  0.589 * pow(1.0+z,  0.20);
    double phi   = -0.729 * pow(1.0+z, -0.08);
    double eta   = -0.243 * pow(1.0+z,  0.27);
    double gamma =  0.864 * pow(1.0+z, -0.01);
    double alpha =  0.368;
    return nu * alpha*(1.0+pow(beta*nu, -2.0*phi))*pow(nu, 2.0*eta)*exp(-0.5*gamma*gsl_pow_2(nu));
}//}}}

static
double bnu_Tinker10(double nu)
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

static
double _dndlogM(all_data *d, int z_index, int M_index, double *bias)
{//{{{
    double sigma_squared = d->pwr->ssq[M_index][0];
    double sigma_squared_prime = d->pwr->ssq[M_index][1];
    double nu = 1.686/sqrt(d->c->Dsq[z_index] * sigma_squared);

    double fnu = fnu_Tinker10(nu, d->n->gr->zgrid[z_index]);
    double hmf = -fnu * d->c->rho_m_0 * sigma_squared_prime
                 / (2.0 * sigma_squared * d->n->gr->Mgrid[M_index]);
    *bias = bnu_Tinker10(nu);
    return hmf;
}//}}}

static
void create_dndlogM(all_data *d)
{//{{{
    printf("\tcreate_dndlogM\n");
    d->h->hmf = (double *)malloc(d->n->gr->Nz * d->n->gr->NM * sizeof(double));
    d->h->bias = (double *)malloc(d->n->gr->Nz * d->n->gr->NM * sizeof(double));
    for (int z_index=0; z_index<d->n->gr->Nz; z_index++)
    {
        for (int M_index=0; M_index<d->n->gr->NM; M_index++)
        {
            d->h->hmf[z_index*d->n->gr->NM + M_index] =
                _dndlogM(d, z_index, M_index,
                         d->h->bias+(z_index*d->n->gr->NM+M_index));
        }
    }
}//}}}

void init_halo_model(all_data *d)
{//{{{
    if (d->h->inited_halo) { return; }
    printf("In halo_model.h -> init_halo_model :\n");
    create_c_of_y(d);
    create_dndlogM(d);
    d->h->inited_halo = 1;
}//}}}

