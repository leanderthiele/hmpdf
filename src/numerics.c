#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <complex.h>

#include <gsl/gsl_integration.h>

#include "configs.h"
#include "utils.h"
#include "data.h"
#include "numerics.h"

int
null_numerics(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

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

    CHECKERR
    return hmpdf_status;
}//}}}

int
reset_numerics(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

    if (d->n->zgrid != NULL) { free(d->n->zgrid); }
    if (d->n->zweights != NULL) { free(d->n->zweights); }
    if (d->n->Mgrid != NULL) { free(d->n->Mgrid); }
    if (d->n->Mweights != NULL) { free(d->n->Mweights); }
    if (d->n->signalgrid != NULL) { free(d->n->signalgrid); }
    if (d->n->signalgrid_noisy != NULL) { free(d->n->signalgrid_noisy); }
    if (d->n->lambdagrid != NULL) { free(d->n->lambdagrid); }
    if (d->n->phigrid != NULL) { free(d->n->phigrid); }
    if (d->n->phiweights != NULL) { free(d->n->phiweights); }

    CHECKERR
    return hmpdf_status;
}//}}}

static int
weight_fct(hmpdf_integr_mode_e m,
           double a, double b, double alpha, double beta,
           double x, double *out)
{//{{{
    int hmpdf_status = 0;

    switch (m)
    {
        case hmpdf_legendre    : *out = 1.0; break;
        case hmpdf_chebyshev   : *out = 1.0 / sqrt((b-x)*(x-a)); break;
        case hmpdf_gegenbauer  : *out = pow((b-x)*(x-a), alpha); break;
        case hmpdf_jacobi      : *out = pow(b-x, alpha) * pow(x-a, beta); break;
        case hmpdf_laguerre    : *out = pow(x-a, alpha) * exp(-b*(x-a)); break;
        case hmpdf_hermite     : *out = pow(fabs(x-a), alpha) * exp(-b*gsl_pow_2(x-a)); break;
        case hmpdf_exponential : *out = pow(fabs(x-0.5*(a+b)), alpha); break;
        case hmpdf_rational    : *out = pow(x-a, alpha) * pow(x+b, beta); break;
        case hmpdf_chebyshev2  : *out = sqrt((b-x)*(x-a)); break;
        default                : *out = 0.0; // to avoid maybe-unitialized
                                 fprintf(stderr, "Error: Unknown gsl_integration_fixed_type.\n");
                                 fflush(stderr);
                                 ERRLOC
                                 hmpdf_status |= 1;
                                 break;
    }

    CHECKERR
    return hmpdf_status;
}//}}}

static int
gauss_fixed_point(hmpdf_integr_mode_e m, int N,
                  double a, double b, double alpha, double beta,
                  double *nodes, double *weights,
                  int neutralize_weights)
{//{{{
    int hmpdf_status = 0;

    const gsl_integration_fixed_type *T;
    switch (m)
    {
        case hmpdf_legendre    : T = gsl_integration_fixed_legendre; break;
        case hmpdf_chebyshev   : T = gsl_integration_fixed_chebyshev; break;
        case hmpdf_gegenbauer  : T = gsl_integration_fixed_gegenbauer; break;
        case hmpdf_jacobi      : T = gsl_integration_fixed_jacobi; break;
        case hmpdf_laguerre    : T = gsl_integration_fixed_laguerre; break;
        case hmpdf_hermite     : T = gsl_integration_fixed_hermite; break;
        case hmpdf_exponential : T = gsl_integration_fixed_exponential; break;
        case hmpdf_rational    : T = gsl_integration_fixed_rational; break;
        case hmpdf_chebyshev2  : T = gsl_integration_fixed_chebyshev2; break;
        default                : T = NULL;
                                 fprintf(stderr, "Error: Unknown gsl_integration_fixed_type.\n");
                                 fflush(stderr);
                                 ERRLOC;
                                 hmpdf_status = 1;
                                 return hmpdf_status;
    }
    SAFEALLOC(gsl_integration_fixed_workspace *, ws,
              gsl_integration_fixed_alloc(T, N, a, b, alpha, beta))
    double *_nodes = gsl_integration_fixed_nodes(ws);
    double *_weights = gsl_integration_fixed_weights(ws);
    for (int ii=0; ii<N; ii++)
    {
        nodes[ii] = _nodes[ii];
        if (neutralize_weights)
        {
            double temp;
            SAFEHMPDF(weight_fct(m, a, b, alpha, beta, nodes[ii], &temp))
            weights[ii] = _weights[ii] / temp;
        }
        else
        {
            weights[ii] = _weights[ii];
        }
    }
    gsl_integration_fixed_free(ws);

    CHECKERR
    return hmpdf_status;
}//}}}

static int
create_grids(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

    fprintf(stdout, "\tcreate_grids\n");
    fflush(stdout);
    SAFEALLOC(, d->n->zgrid,    malloc(d->n->Nz * sizeof(double)))
    SAFEALLOC(, d->n->zweights, malloc(d->n->Nz * sizeof(double)))
    SAFEHMPDF(gauss_fixed_point(d->n->zintegr_type, d->n->Nz,
                                d->n->zmin, d->n->zmax,
                                d->n->zintegr_alpha,
                                d->n->zintegr_beta,
                                d->n->zgrid, d->n->zweights,
                                1/*neutralize weights*/))

    SAFEALLOC(, d->n->Mgrid,    malloc(d->n->NM * sizeof(double)))
    SAFEALLOC(, d->n->Mweights, malloc(d->n->NM * sizeof(double)))
    SAFEHMPDF(gauss_fixed_point(d->n->Mintegr_type, d->n->NM,
                                log(d->n->Mmin), log(d->n->Mmax),
                                d->n->Mintegr_alpha,
                                d->n->Mintegr_beta,
                                d->n->Mgrid, d->n->Mweights,
                                1/*neutralize weights*/))

    for (int ii=0; ii<d->n->NM; ii++)
    {
        d->n->Mgrid[ii] = exp(d->n->Mgrid[ii]);
    }

    SAFEALLOC(, d->n->signalgrid, malloc(d->n->Nsignal * sizeof(double)))
    linspace(d->n->Nsignal, d->n->signalmin, d->n->signalmax,
             d->n->signalgrid);

    SAFEALLOC(, d->n->lambdagrid, malloc((d->n->Nsignal/2+1) * sizeof(double)))
    linspace(d->n->Nsignal/2+1,
             0.0, M_PI/(d->n->signalgrid[1]-d->n->signalgrid[0]),
             d->n->lambdagrid);

    CHECKERR
    return hmpdf_status;
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

int
init_numerics(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

    if (d->n->inited_numerics) { return hmpdf_status; }
    fprintf(stdout, "In numerics.h -> init_numerics :\n");
    fflush(stdout);

    SAFEHMPDF(create_grids(d))

    d->n->inited_numerics = 1;

    CHECKERR
    return hmpdf_status;
}//}}}

