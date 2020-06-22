#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>

#include "utils.h"
#include "data.h"
#include "numerics.h"
#include "noise.h"

#include "hmpdf.h"

int
null_noise(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

    d->ns->toepl = NULL;

    CHECKERR
    return hmpdf_status;
}//}}}

int
reset_noise(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

    if (d->ns->toepl != NULL) { free(d->ns->toepl); }

    CHECKERR
    return hmpdf_status;
}//}}}

int
create_noisy_grids(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

    fprintf(stdout, "\tcreate_noisy_grids\n");
    fflush(stdout);

    // can in principle make this a user setting
    d->ns->len_kernel = d->n->Nsignal;
    // construct the new signal grid
    d->n->Nsignal_noisy = d->n->Nsignal+2*d->ns->len_kernel;
    SAFEALLOC(, d->n->signalgrid_noisy, malloc(d->n->Nsignal_noisy
                                               * sizeof(double)))
    double extra_signal = (double)(d->ns->len_kernel)/(double)(d->n->Nsignal-1)
                          *(d->n->signalmax - d->n->signalmin);
    double smin = d->n->signalmin - extra_signal;
    double smax = d->n->signalmax + extra_signal;
    linspace(d->n->Nsignal_noisy, smin, smax, d->n->signalgrid_noisy);

    CHECKERR
    return hmpdf_status;
}//}}}

int
create_toepl(hmpdf_obj *d)
// creates Toeplitz matrix nrows=Nsignal, ncols=Nsignal_noisy
{//{{{
    int hmpdf_status = 0;

    fprintf(stdout, "\tcreate_toepl\n");
    fflush(stdout);

    SAFEALLOC(, d->ns->toepl, malloc(d->n->Nsignal * d->n->Nsignal_noisy
                                     * sizeof(double)))
    zero_real(d->n->Nsignal*d->n->Nsignal_noisy, d->ns->toepl);
    double sigma = d->ns->noise
                   / (d->n->signalgrid[1] - d->n->signalgrid[0]);
    // fill the first row
    for (int ii=-d->ns->len_kernel; ii<=d->ns->len_kernel; ii++)
    {
        d->ns->toepl[ii+d->ns->len_kernel] = exp(-0.5*gsl_pow_2((double)(ii)/sigma))
                                      /sqrt(2.0*M_PI)/sigma;
    }
    // fill the remaining rows
    for (int ii=1; ii<d->n->Nsignal; ii++)
    {
        memcpy(d->ns->toepl + ii*(d->n->Nsignal_noisy+1), d->ns->toepl,
               (2*d->ns->len_kernel+1) * sizeof(double));
    }

    CHECKERR
    return hmpdf_status;
}//}}}

int
noise_vect(hmpdf_obj *d, double *in, double *out)
// assumes len(in) = Nsignal, len(out) = Nsignal_noisy
{//{{{
    int hmpdf_status = 0;

    if (d->ns->toepl == NULL)
    {
        ERRLOC
        fprintf(stderr, "Error : Toeplitz matrix not computed. "
                        "Invalid call to noise_vect.\n");
        fflush(stderr);
        hmpdf_status |= 1;
        return hmpdf_status;
    }

    cblas_dgemv(CblasRowMajor, CblasTrans/*toepl matrix needs to be transposed*/,
                d->n->Nsignal/*rows*/, d->n->Nsignal_noisy/*cols*/,
                1.0/*alpha*/, d->ns->toepl/*matrix A*/, d->n->Nsignal_noisy/*lda*/,
                in/*input vector X*/, 1/*stride of X*/, 0.0/*beta*/,
                out/*output vector Y*/, 1/*stride of Y*/);
    
    CHECKERR
    return hmpdf_status;
}//}}}

int
noise_matr(hmpdf_obj *d, double *in, double *out)
// assumes [in] = Nsignal*Nsignal, [out] = Nsignal_noisy*Nsignal_noisy
{//{{{
    int hmpdf_status = 0;

    if (d->ns->toepl == NULL)
    {
        ERRLOC
        fprintf(stderr, "Error : Toeplitz matrix not computed. "
                        "Invalid call to noise_matr.\n");
        fflush(stderr);
        hmpdf_status |= 1;
        return hmpdf_status;
    }

    // no aliasing allowed, so we need intermediate storage
    SAFEALLOC(double *, temp, malloc(d->n->Nsignal * d->n->Nsignal_noisy
                                     * sizeof(double)))

    // multiply from the left with Toeplitz matrix
    cblas_dgemm(CblasRowMajor, CblasTrans/*toepl matrix needs to be transposed*/,
                CblasNoTrans/*the right matrix is not transposed*/,
                d->n->Nsignal_noisy/*rows of left matrix & output*/,
                d->n->Nsignal/*cols of right matrix & output*/,
                d->n->Nsignal/*cols of left matrix, rows of right matrix*/,
                1.0/*alpha*/, d->ns->toepl/*left matrix*/,
                d->n->Nsignal_noisy/*lda*/, in/*right matrix*/,
                d->n->Nsignal/*ldb*/, 0.0/*beta*/, temp/*output matrix*/,
                d->n->Nsignal/*ldc*/);

    // multiply from the right with Toeplitz matrix (transposed)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                d->n->Nsignal_noisy, d->n->Nsignal_noisy, d->n->Nsignal,
                1.0, temp, d->n->Nsignal, d->ns->toepl, d->n->Nsignal_noisy,
                0.0, out, d->n->Nsignal_noisy);

    free(temp);

    CHECKERR
    return hmpdf_status;
}//}}}

int
init_noise(hmpdf_obj *d)
{//{{{
    int hmpdf_status = 0;

    if (d->ns->noise > 0.0)
    {
        fprintf(stdout, "In noise.h -> init_noise.\n");
        fflush(stdout);

        SAFEHMPDF(create_noisy_grids(d))
        SAFEHMPDF(create_toepl(d))
    }

    CHECKERR
    return hmpdf_status;
}//}}}
