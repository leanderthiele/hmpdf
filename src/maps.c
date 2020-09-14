#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#ifdef _OPENMP
#   include <omp.h>
#endif

#include <fftw3.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "configs.h"
#include "utils.h"
#include "object.h"
#include "filter.h"
#include "profiles.h"
#include "maps.h"

#include "hmpdf.h"

int
null_maps(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->m->created_map = 0;
    d->m->ellgrid = NULL;
    d->m->map_real = NULL;
    d->m->p_r2c = NULL;
    d->m->p_c2r = NULL;
    d->m->created_map_ws = 0;
    d->m->ws = NULL;

    ENDFCT
}//}}}

int
reset_maps(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(2, "\treset_maps\n");

    if (d->m->ellgrid != NULL) { free(d->m->ellgrid); }
    if (d->m->map_real != NULL) { fftw_free(d->m->map_real); }
    if (d->m->p_r2c != NULL) { fftw_destroy_plan(*(d->m->p_r2c)); free(d->m->p_r2c); }
    if (d->m->p_c2r != NULL) { fftw_destroy_plan(*(d->m->p_c2r)); free(d->m->p_c2r); }
    if (d->m->ws != NULL)
    {
        for (int ii=0; ii<d->m->Nws; ii++)
        {
            if (d->m->ws[ii] != NULL)
            {
                if (d->m->ws[ii]->map != NULL) { free(d->m->ws[ii]->map); };
                if (d->m->ws[ii]->pos != NULL) { free(d->m->ws[ii]->pos); };
                if (d->m->ws[ii]->buf != NULL) { free(d->m->ws[ii]->buf); };
                if (d->m->ws[ii]->rng != NULL) { gsl_rng_free(d->m->ws[ii]->rng); }
                free(d->m->ws[ii]);
            }
        }
        free(d->m->ws);
    }

    ENDFCT
}//}}}

#define NEWMAPWS_SAFEALLOC(var, expr)  \
    do {                               \
        var = expr;                    \
        if (UNLIKELY(!(var)))          \
        {                              \
            if (ws->map != NULL)       \
            { free(ws->map); }         \
            if (ws->pos != NULL)       \
            { free(ws->pos); }         \
            if (ws->buf != NULL)       \
            { free(ws->buf); }         \
            free(*out);                \
            if (ws->rng != NULL)       \
            { gsl_rng_free(ws->rng); } \
            return 1;                  \
        }                              \
    } while (0)

static int
new_map_ws(hmpdf_obj *d, map_ws **out)
// allocates a new map workspace
{//{{{
    STARTFCT

    SAFEALLOC(*out, malloc(sizeof(map_ws)));

    map_ws *ws = *out; // for convenience

    // initialize to NULL so we can free realiably in case an alloc fails
    ws->map = NULL;
    ws->pos = NULL;
    ws->buf = NULL;
    ws->rng = NULL;

    NEWMAPWS_SAFEALLOC(ws->map, malloc(d->m->Nside * d->m->Nside
                                       * sizeof(double)));
    NEWMAPWS_SAFEALLOC(ws->pos, malloc(d->m->buflen
                                       * sizeof(double)));
    NEWMAPWS_SAFEALLOC(ws->buf, malloc(d->m->buflen
                                       * sizeof(double)));
    NEWMAPWS_SAFEALLOC(ws->rng, gsl_rng_alloc(gsl_rng_ranlux));

    ENDFCT
}//}}}

#undef NEWMAPWS_SAFEALLOC

static int
create_map_ws(hmpdf_obj *d)
// allocates as many workspaces as necessary/possible
{//{{{
    STARTFCT

    if (d->m->created_map_ws) { return 0; }

    HMPDFPRINT(2, "\tcreate_map_ws\n");
    HMPDFPRINT(3, "\t\ttrying to allocate workspaces for %d threads.\n", d->Ncores);

    SAFEALLOC(d->m->ws, malloc(d->Ncores * sizeof(map_ws *)));
    SETARRNULL(d->m->ws, d->Ncores);
    d->m->Nws = 0;
    for (int ii=0; ii<d->Ncores; ii++)
    {
        int alloc_failed = new_map_ws(d, d->m->ws+ii);
        if (alloc_failed)
        {
            d->m->ws[ii] = NULL;
            break;
        }
        else
        {
            ++d->m->Nws;
        }
    }

    if (d->m->Nws < d->Ncores)
    {
        HMPDFPRINT(1, "Allocated only %d workspaces, "
                      "because memory ran out.\n", d->m->Nws);
    }

    HMPDFCHECK(d->m->Nws<1, "Failed to allocate any workspaces.");

    d->m->created_map_ws = 1;

    ENDFCT
}//}}}

static int
reset_map_ws(hmpdf_obj *d, map_ws *ws)
{//{{{
    STARTFCT

    // seed the random number generator
    //     be careful that the seed is as close to random as possible
    gsl_rng_set(ws->rng, (unsigned long)time(NULL)
                         + (unsigned long)rand()
                         + (unsigned long)(ws->buf));

    zero_real(d->m->Nside * d->m->Nside, ws->map);

    ENDFCT
}//}}}

static int
fill_buf(hmpdf_obj *d, int z_index, int M_index, map_ws *ws)
// creates a map of the given object in the buffer
{//{{{
    STARTFCT

    // theta_out in units of the pixel spacing
    double tout = d->p->profiles[z_index][M_index][0]
                  / d->f->pixelside;

    // compute how large this specific map needs to be
    //     the map is of size (2*w+1)^2
    long w = (long)round(tout) + 1;
    ws->bufside = 2 * w + 1;

    // draw random displacement of the center of the halo
    double dx = gsl_rng_uniform(ws->rng) - 0.5;
    double dy = gsl_rng_uniform(ws->rng) - 0.5;

    // compute the dimensions of the buffer
    //     buf = [ w, w, ld, ld ]
    long ld = 2 * d->m->pxlgrid + 1; // innermost

    // safety check
    HMPDFCHECK(ws->bufside * ws->bufside * ld * ld >= d->m->buflen,
               "did not allocate enough buffer memory.");

    // fill the position grid
    long idx = 0;
    for (long xx= -w; xx<= w; xx++)
    {
        for (long yy= -w; yy<= w; yy++)
        {
            for (long xp= -d->m->pxlgrid;
                 xp<= d->m->pxlgrid; xp++)
            {
                for (long yp= -d->m->pxlgrid;
                     yp<= d->m->pxlgrid; yp++)
                {
                    double xpos = (double)xx + (double)(2*xp)/(double)ld + dx;
                    double ypos = (double)yy + (double)(2*yp)/(double)ld + dy;
                    // compute distance from center in units of the outer radius
                    ws->pos[idx] = hypot(xpos, ypos) / tout;
                    ++idx;
                }
            }
        }
    }

    // evaluate the profile interpolator
    SAFEHMPDF(s_of_t(d, z_index, M_index,
                     ws->bufside * ws->bufside * ld * ld,
                     ws->pos, ws->buf));

    // compute the average over pixels
    for (long ii=0; ii< ws->bufside * ws->bufside; ii++)
    {
        double avg = 0.0;
        for (long jj=0; jj< ld * ld; jj++)
        {
            avg += ws->buf[ii*ld*ld + jj];
        }
        avg /= (double)(ld * ld);
        ws->buf[ii] = avg;
    }

    ENDFCT
}//}}}

// convenience macro to reduce typing
#define INNERLOOP_OP                     \
    ws->map[ixx*d->m->Nside + iyy]       \
        += ws->buf[xx*ws->bufside + yy];

inline static void 
add_buf_inner_loop(hmpdf_obj *d, map_ws *ws, long y0, long xx, long ixx)
{//{{{
    for (long yy=0, iyy=y0;
         yy< GSL_MIN(ws->bufside, d->m->Nside - y0);
         yy++, iyy++)
    {
        INNERLOOP_OP
    }
    for (long yy= GSL_MIN(ws->bufside, d->m->Nside - y0), iyy=0;
         UNLIKELY(yy< ws->bufside);
         yy++, iyy++)
    {
        INNERLOOP_OP
    }
}//}}}

#define OUTERLOOP_OP                        \
    add_buf_inner_loop(d, ws, y0, xx, ixx);

static int 
add_buf(hmpdf_obj *d, map_ws *ws)
// picks random position in the map and adds buffer map once,
//     satisfies periodic boundary conditions
{//{{{
    STARTFCT

    // pick a random point in the map
    long x0 = gsl_rng_uniform_int(ws->rng, d->m->Nside);
    long y0 = gsl_rng_uniform_int(ws->rng, d->m->Nside);

    // add the pixel values from the buffer
    //     we 'unroll' the loops slightly for better efficiency
    //     with the period boundary conditions
    // A word on notation : xx, yy are coordinates in this map
    //                             (the one stored in ws->buf)
    //                      ixx, iyy are coordinates in the total map
    //                             (the one stored in ws->map)
    for (long xx=0, ixx=x0;
         xx< GSL_MIN(ws->bufside, d->m->Nside - x0);
         xx++, ixx++)
    {
        OUTERLOOP_OP
    }
    for (long xx= GSL_MIN(ws->bufside, d->m->Nside - x0), ixx=0;
         UNLIKELY(xx< ws->bufside); xx++, ixx++)
    {
        OUTERLOOP_OP
    }

    ENDFCT
}//}}}

#undef OUTERLOOP_OP
#undef INNERLOOP_OP

static int
draw_N_halos(hmpdf_obj *d, int z_index, int M_index, map_ws *ws, unsigned *N)
// draws the number of halos in the given bin
{//{{{
    STARTFCT
    
    // TODO think about how we compute this!
    double dlogM;
    if (M_index == 0)
    {
        dlogM = log(d->n->Mgrid[M_index+1])
                - log(d->n->Mgrid[M_index]);
    }
    else if (M_index == d->n->NM - 1)
    {
        dlogM = log(d->n->Mgrid[M_index])
                - log(d->n->Mgrid[M_index-1]);
    }
    else
    {
        dlogM = 0.5 * (log(d->n->Mgrid[M_index+1]
                       - log(d->n->Mgrid[M_index-1])));
    }

    // expected number of halos in this bin
    double n = d->h->hmf[z_index][M_index] // dn_3d / dlogM
               * gsl_pow_2(d->c->comoving[z_index])
               / d->c->hubble[z_index]
               * dlogM
               * d->m->area;

    *N = gsl_ran_poisson(ws->rng, n);

    ENDFCT
}//}}}

static int
do_this_bin(hmpdf_obj *d, int z_index, int M_index, map_ws *ws)
// draws random integer from correct distribution
// if ==0, return
// else, fill_buf and then integer x add_buf
{//{{{
    STARTFCT

    unsigned N;
    SAFEHMPDF(draw_N_halos(d, z_index, M_index, ws, &N));

    if (N == 0)
    {
        return 0;
    }
    else
    {
        SAFEHMPDF(fill_buf(d, z_index, M_index, ws));
        for (unsigned ii=0; ii<N; ii++)
        {
            SAFEHMPDF(add_buf(d, ws));
        }
    }

    ENDFCT
}//}}}

static int
grf(hmpdf_obj *d)
// adds random GRF realization to the Fourier space map
{//{{{
    STARTFCT

    ENDFCT
}//}}}

static int
create_map(hmpdf_obj *d)
// loops in parallel over flattened, shuffled array of (z_index, M_index)
// for each step, call do_this_bin
//
// performs r2c FFT
// multiplies with filters (CAREFUL: do not include the pixel window fct!)
// adds random GRF realization with the given noise power spectrum
// performs c2r FFT
{//{{{
    STARTFCT

    if (d->m->created_map) { return 0; }

    // reset the workspaces
    for (int ii=0; ii<d->m->Nws; ii++)
    {
        SAFEHMPDF(reset_map_ws(d, d->m->ws[ii]));
    }

    // create the array of bins
    int *bins;
    SAFEALLOC(bins, malloc(d->n->Nz * d->n->NM * sizeof(int)));
    for (int ii=0; ii<d->n->Nz * d->n->NM; ii++)
    {
        bins[ii] = ii;
    }
    // shuffle to equalize load
    gsl_ran_shuffle(d->m->ws[0]->rng, bins, d->n->Nz * d->n->NM, sizeof(int));

    // perform the loop
    #ifdef _OPENMP
    #   pragma omp parallel for num_threads(d->m->Nws) schedule(dynamic)
    #endif
    for (int ii=0; ii<d->n->Nz * d->n->NM; ii++)
    {
        CONTINUE_IF_ERR

        int z_index = bins[ii] / d->n->NM;
        int M_index = bins[ii] % d->n->NM;
        SAFEHMPDF_NORETURN(do_this_bin(d, z_index, M_index,
                                       d->m->ws[THIS_THREAD]));
    }

    free(bins);

    d->m->created_map = 1;

    ENDFCT
}//}}}

static int
find_Nside(hmpdf_obj *d, long *N)
{//{{{
    STARTFCT

    double map_side = sqrt(d->m->area);
    *N = (long)round(map_side/d->f->pixelside);

    ENDFCT
}//}}}

static int
find_buflen(hmpdf_obj *d, long *N)
// finds the maximum buffer length we could need
//     corresponds to the halo with the largest angular extent
{//{{{
    STARTFCT

    double max_t_out = 0.0;
    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        for (int M_index=0; M_index<d->n->NM; M_index++)
        {
            max_t_out = GSL_MAX(max_t_out,
                                d->p->profiles[z_index][M_index][0]);
        }
    }
    *N = (long)round(max_t_out/d->f->pixelside);
    *N += 4; // some safety buffer
    *N *= (long)(2*d->m->pxlgrid+1);
    *N *= *N;

    ENDFCT
}//}}}

static int
prepare_maps(hmpdf_obj *d)
// compute map dimensions
// allocate workspaces
{//{{{
    STARTFCT

    SAFEHMPDF(find_Nside(d, &(d->m->Nside)));
    SAFEHMPDF(find_buflen(d, &(d->m->buflen)));
    SAFEHMPDF(create_map_ws(d));

    ENDFCT
}//}}}

int
hmpdf_get_map_hist(hmpdf_obj *d, int Nbins, double binedges[Nbins+1], double *hist, int new_map)
// if (new_map), create one
// else, if not available, create one
//       else, use the existing one
{//{{{
    STARTFCT

    ENDFCT
}//}}}
