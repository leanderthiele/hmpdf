#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>

#include "utils.h"
#include "object.h"
#include "numerics.h"
#include "noise.h"

#include "hmpdf.h"

int
null_noise(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->ns->toepl = NULL;

    ENDFCT
}//}}}

int
reset_noise(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\treset_noise\n");

    if (d->ns->toepl != NULL) { free(d->ns->toepl); }

    ENDFCT
}//}}}

int
create_noisy_grids(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\tcreate_noisy_grids\n");

    // can in principle make this a user setting
    d->ns->len_kernel = d->n->Nsignal;
    // construct the new signal grid
    d->n->Nsignal_noisy = d->n->Nsignal+2*d->ns->len_kernel;
    SAFEALLOC(d->n->signalgrid_noisy, malloc(d->n->Nsignal_noisy
                                             * sizeof(double)));
    double extra_signal = (double)(d->ns->len_kernel)/(double)(d->n->Nsignal-1)
                          *(d->n->signalmax - d->n->signalmin);
    double smin = d->n->signalmin - extra_signal;
    double smax = d->n->signalmax + extra_signal;
    SAFEHMPDF(linspace(d->n->Nsignal_noisy, smin, smax, d->n->signalgrid_noisy));

    ENDFCT
}//}}}

int
create_toepl(hmpdf_obj *d)
// creates Toeplitz matrix nrows=Nsignal, ncols=Nsignal_noisy
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\tcreate_toepl\n");

    SAFEALLOC(d->ns->toepl, malloc(d->n->Nsignal * d->n->Nsignal_noisy
                                     * sizeof(double)));
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

    ENDFCT
}//}}}

int
noise_vect(hmpdf_obj *d, double *in, double *out)
// assumes len(in) = Nsignal, len(out) = Nsignal_noisy
{//{{{
    STARTFCT

    HMPDFCHECK(d->ns->toepl == NULL, "Toeplitz matrix not computed.");

    cblas_dgemv(CblasRowMajor, CblasTrans/*toepl matrix needs to be transposed*/,
                d->n->Nsignal/*rows*/, d->n->Nsignal_noisy/*cols*/,
                1.0/*alpha*/, d->ns->toepl/*matrix A*/, d->n->Nsignal_noisy/*lda*/,
                in/*input vector X*/, 1/*stride of X*/, 0.0/*beta*/,
                out/*output vector Y*/, 1/*stride of Y*/);
    
    ENDFCT
}//}}}

int
noise_matr(hmpdf_obj *d, double *in, double *out)
// assumes [in] = Nsignal*Nsignal, [out] = Nsignal_noisy*Nsignal_noisy
{//{{{
    STARTFCT

    HMPDFCHECK(d->ns->toepl == NULL, "Toeplitz matrix not computed.");

    // no aliasing allowed, so we need intermediate storage
    double *temp;
    SAFEALLOC(temp, malloc(d->n->Nsignal * d->n->Nsignal_noisy
                                     * sizeof(double)));

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

    ENDFCT
}//}}}

int
init_noise(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->ns->noise > 0.0)
    {
        HMPDFPRINT(1, "init_noise\n");

        SAFEHMPDF(create_noisy_grids(d));
        SAFEHMPDF(create_toepl(d));
    }

    ENDFCT
}//}}}
