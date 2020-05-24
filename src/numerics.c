#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <complex.h>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>

#include "configs.h"
#include "utils.h"
#include "data.h"
#include "numerics.h"

void null_numerics(all_data *d)
{//{{{
    d->n->inited_numerics = 0;
    d->n->zgrid = NULL;
    d->n->zweights = NULL;
    d->n->Mgrid = NULL;
    d->n->Mweights = NULL;
    d->n->signalgrid = NULL;
    d->n->signalgrid_noisy = NULL;
    d->n->lambdagrid = NULL;
    d->n->phigrid = NULL;
    d->n->phiweights = NULL;
}//}}}

void reset_numerics(all_data *d)
{//{{{
    if (d->n->zgrid != NULL) { free(d->n->zgrid); }
    if (d->n->zweights != NULL) { free(d->n->zweights); }
    if (d->n->Mgrid != NULL) { free(d->n->Mgrid); }
    if (d->n->Mweights != NULL) { free(d->n->Mweights); }
    if (d->n->signalgrid != NULL) { free(d->n->signalgrid); }
    if (d->n->signalgrid_noisy != NULL) { free(d->n->signalgrid_noisy); }
    if (d->n->lambdagrid != NULL) { free(d->n->lambdagrid); }
    if (d->n->phigrid != NULL) { free(d->n->phigrid); }
    if (d->n->phiweights != NULL) { free(d->n->phiweights); }
}//}}}

static
double weight_fct(integr_mode m,
                  double a, double b, double alpha, double beta,
                  double x)
{//{{{
    switch (m)
    {
        case legendre    : return 1.0;
        case chebyshev   : return 1.0 / sqrt((b-x)*(x-a));
        case gegenbauer  : return pow((b-x)*(x-a), alpha);
        case jacobi      : return pow(b-x, alpha) * pow(x-a, beta);
        case laguerre    : return pow(x-a, alpha) * exp(-b*(x-a));
        case hermite     : return pow(fabs(x-a), alpha) * exp(-b*gsl_pow_2(x-a));
        case exponential : return pow(fabs(x-0.5*(a+b)), alpha);
        case rational    : return pow(x-a, alpha) * pow(x+b, beta);
        case chebyshev2  : return sqrt((b-x)*(x-a));
        default          : printf("Unknown gsl_integration_fixed_type.\n");
                           return 0.0;
    }
}//}}}

static
void gauss_fixed_point(integr_mode m, int N,
                       double a, double b, double alpha, double beta,
                       double *nodes, double *weights,
                       int neutralize_weights)
{//{{{
    const gsl_integration_fixed_type *T;
    switch (m)
    {
        case legendre    : T = gsl_integration_fixed_legendre; break;
        case chebyshev   : T = gsl_integration_fixed_chebyshev; break;
        case gegenbauer  : T = gsl_integration_fixed_gegenbauer; break;
        case jacobi      : T = gsl_integration_fixed_jacobi; break;
        case laguerre    : T = gsl_integration_fixed_laguerre; break;
        case hermite     : T = gsl_integration_fixed_hermite; break;
        case exponential : T = gsl_integration_fixed_exponential; break;
        case rational    : T = gsl_integration_fixed_rational; break;
        case chebyshev2  : T = gsl_integration_fixed_chebyshev2; break;
        default          : printf("Unknown gsl_integration_fixed_type.\n");
                           return;
    }
    gsl_integration_fixed_workspace *ws
        = gsl_integration_fixed_alloc(T, N, a, b, alpha, beta);
    double *_nodes = gsl_integration_fixed_nodes(ws);
    double *_weights = gsl_integration_fixed_weights(ws);
    for (int ii=0; ii<N; ii++)
    {
        nodes[ii] = _nodes[ii];
        if (neutralize_weights)
        {
            weights[ii] = _weights[ii]
                          / weight_fct(m, a, b, alpha, beta, nodes[ii]);
        }
        else
        {
            weights[ii] = _weights[ii];
        }
    }
    gsl_integration_fixed_free(ws);
}//}}}

static
void create_grids(all_data *d)
{//{{{
    printf("\tcreate_grids\n");
    d->n->zgrid = (double *)malloc(d->n->Nz * sizeof(double));
    d->n->zweights = (double *)malloc(d->n->Nz * sizeof(double));
    gauss_fixed_point(d->n->zintegr_type, d->n->Nz,
                      d->n->zmin, d->n->zmax,
                      d->n->zintegr_alpha, d->n->zintegr_beta,
                      d->n->zgrid, d->n->zweights,
                      1/*neutralize weights*/);

    d->n->Mgrid = (double *)malloc(d->n->NM * sizeof(double));
    d->n->Mweights = (double *)malloc(d->n->NM * sizeof(double));
    gauss_fixed_point(d->n->Mintegr_type, d->n->NM,
                      log(d->n->Mmin), log(d->n->Mmax),
                      d->n->Mintegr_alpha, d->n->Mintegr_beta,
                      d->n->Mgrid, d->n->Mweights,
                      1/*neutralize weights*/);
    for (int ii=0; ii<d->n->NM; ii++)
    {
        d->n->Mgrid[ii] = exp(d->n->Mgrid[ii]);
    }

    d->n->signalgrid = (double *)malloc(d->n->Nsignal * sizeof(double));
    linspace(d->n->Nsignal, d->n->signalmin, d->n->signalmax,
             d->n->signalgrid);

    d->n->lambdagrid = (double *)malloc((d->n->Nsignal/2+1) * sizeof(double));
    linspace(d->n->Nsignal/2+1,
             0.0, M_PI/(d->n->signalgrid[1]-d->n->signalgrid[0]),
             d->n->lambdagrid);
}//}}}

static
double simps_real(int N, double dx, int stride, double *f)
// real Simpson integration, N must be odd (even number of intervals)
{//{{{
    assert(GSL_IS_ODD(N));
    --N;
    double out = f[0] + 4.0*f[(N-1)*stride] + f[N*stride];
    for (int ii=1; ii<N-2; ii+=2)
    {
        out += 4.0*f[ii*stride] + 2.0*f[(ii+1)*stride];
    }
    return out * dx / 3.0;
}//}}}

static
complex simps_comp(int N, double dx, int stride, complex *f)
// complex Simpson integration, N must be odd (even number of intervals)
{//{{{
    assert(GSL_IS_ODD(N));
    --N;
    complex out = f[0] + 4.0*f[(N-1)*stride] + f[N*stride];
    for (int ii=1; ii<N-2; ii+=2)
    {
        out += 4.0*f[ii*stride] + 2.0*f[(ii+1)*stride];
    }
    return out * dx / 3.0;
}//}}}

static
double romb_real(int k, double dx, int stride, double *f)
// real Romberg integration. N = 2^k+1 sample points, real integrand
{//{{{
    int N = 1<<k;
    double R1[k];
    double R2[k];
    double *Rc = &R1[0];
    double *Rp = &R2[0];
    double h = (double)(N) * dx;
    Rp[0] = 0.5 * h * (f[0] + f[N*stride]);

    for (int nn=1; nn <= k; nn++)
    {
        h *= 0.5;
        double temp = 0.0;
        for (int jj=1; jj<=(1<<(nn-1)); jj++)
        {
            temp += f[((2*jj-1)*(N/(1<<nn)))*stride];
        }
        Rc[0] = h * temp + 0.5 * Rp[0];
        for (int jj=1; jj<=nn; jj++)
        {
            Rc[jj] = ((1<<(2*jj))*Rc[jj-1] - Rp[jj-1]) / (double)((1<<(2*jj)) - 1);
        }
        double *Rt = Rp;
        Rp = Rc;
        Rc = Rt;
    }
    return Rp[k - 1];
}//}}}

static
complex romb_comp(int k, double dx, int stride, complex *f)
// complex Romberg integration. N = 2^k subdivisions, complex integrand
{//{{{
    int N = 1<<k;
    complex R1[k];
    complex R2[k];
    complex *Rc = &R1[0];
    complex *Rp = &R2[0];
    double h = (double)(N) * dx;
    Rp[0] = 0.5 * h * (f[0] + f[N*stride]);

    for (int nn=1; nn <= k; nn++)
    {
        h *= 0.5;
        complex temp = 0.0;
        for (int jj=1; jj<=(1<<(nn-1)); jj++)
        {
            temp += f[((2*jj-1)*(N/(1<<nn)))*stride];
        }
        Rc[0] = h * temp + 0.5 * Rp[0];
        for (int jj=1; jj<=nn; jj++)
        {
            Rc[jj] = ((1<<(2*jj))*Rc[jj-1] - Rp[jj-1]) / (double)((1<<(2*jj)) - 1);
        }
        complex *Rt = Rp;
        Rp = Rc;
        Rc = Rt;
    }
    return Rp[k - 1];
}//}}}

double integr_real(int N, double dx, int stride, double *f)
{//{{{
    assert(GSL_IS_ODD(N));
    int k;
    if (ispwr2(N-1, &k))
    {
        return romb_real(k, dx, stride, f);
    }
    else
    {
        return simps_real(N, dx, stride, f);
    }
}//}}}

complex integr_comp(int N, double dx, int stride, complex *f)
{//{{{
    assert(GSL_IS_ODD(N));
    int k;
    if (ispwr2(N-1, &k))
    {
        return romb_comp(k, dx, stride, f);
    }
    else
    {
        return simps_comp(N, dx, stride, f);
    }
}//}}}

void init_numerics(all_data *d)
{//{{{
    if (d->n->inited_numerics) { return; }
    printf("In numerics.h -> init_numerics :\n");
    create_grids(d);

    // TODO this is somewhat awkward here
    d->n->zsource = (d->n->zsource < 0.0) ?
                    d->n->zmax+0.001
                    : d->n->zsource;

    d->n->inited_numerics = 1;
}//}}}

