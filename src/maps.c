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
#include <gsl/gsl_histogram.h>

#include "configs.h"
#include "utils.h"
#include "object.h"
#include "filter.h"
#include "profiles.h"
#include "onepoint.h"
#include "maps.h"

#include "hmpdf.h"

int
null_maps(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->m->created_sidelengths = 0;

    d->m->created_mem = 0;

    d->m->created_ellgrid = 0;
    d->m->ellgrid = NULL;

    d->m->created_map = 0;
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
    if (d->m->map_real != NULL)
    {
        if (d->m->need_ft)
        {
            fftw_free(d->m->map_real);
        }
        else
        {
            free(d->m->map_real);
        }
    }
    if (d->m->p_r2c != NULL) { fftw_destroy_plan(*(d->m->p_r2c)); free(d->m->p_r2c); }
    if (d->m->p_c2r != NULL) { fftw_destroy_plan(*(d->m->p_c2r)); free(d->m->p_c2r); }
    if (d->m->ws != NULL)
    {
        for (int ii=0; ii<d->m->Nws; ii++)
        {
            if (d->m->ws[ii] != NULL)
            {
                if (d->m->ws[ii]->map != NULL)
                {
                    if (d->m->ws[ii]->for_fft)
                    {
                        fftw_free(d->m->ws[ii]->map);
                    }
                    else
                    {
                        free(d->m->ws[ii]->map);
                    }
                }
                if (d->m->ws[ii]->pos != NULL) { free(d->m->ws[ii]->pos); }
                if (d->m->ws[ii]->buf != NULL) { free(d->m->ws[ii]->buf); }
                if (d->m->ws[ii]->rng != NULL) { gsl_rng_free(d->m->ws[ii]->rng); }
                if (d->m->ws[ii]->p_r2c != NULL)
                {
                    fftw_destroy_plan(*(d->m->ws[ii]->p_r2c));
                    free(d->m->ws[ii]->p_r2c);
                }
                if (d->m->ws[ii]->p_c2r != NULL)
                {
                    fftw_destroy_plan(*(d->m->ws[ii]->p_c2r));
                    free(d->m->ws[ii]->p_c2r);
                }
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
new_map_ws(hmpdf_obj *d, int idx, map_ws **out)
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
    ws->p_r2c = NULL;
    ws->p_c2r = NULL;

    if (idx == 0 && d->f->has_z_dependent)
    {
        ws->for_fft = 1;
    }
    else
    {
        ws->for_fft = 0;
    }

    NEWMAPWS_SAFEALLOC(ws->pos, malloc(d->m->buflen
                                       * sizeof(double)));
    NEWMAPWS_SAFEALLOC(ws->buf, malloc(d->m->buflen
                                       * sizeof(double)));
    NEWMAPWS_SAFEALLOC(ws->rng, gsl_rng_alloc(gsl_rng_taus));

    if (ws->for_fft)
    {
        NEWMAPWS_SAFEALLOC(ws->map, fftw_malloc(d->m->Nside * (d->m->Nside+2)
                                                * sizeof(double)));
        ws->map_comp = (double complex *)ws->map;
        NEWMAPWS_SAFEALLOC(ws->p_r2c, malloc(sizeof(fftw_plan)));
        NEWMAPWS_SAFEALLOC(ws->p_c2r, malloc(sizeof(fftw_plan)));
        *(ws->p_r2c) = fftw_plan_dft_r2c_2d(d->m->Nside, d->m->Nside,
                                            ws->map, ws->map_comp, FFTW_MEASURE);
        *(ws->p_c2r) = fftw_plan_dft_c2r_2d(d->m->Nside, d->m->Nside,
                                            ws->map_comp, ws->map, FFTW_MEASURE);
    }
    else
    {
        NEWMAPWS_SAFEALLOC(ws->map, malloc(d->m->Nside * d->m->Nside
                                           * sizeof(double)));
    }

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
        int alloc_failed = new_map_ws(d, ii, d->m->ws+ii);
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
create_ellgrid(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->m->created_ellgrid) { return 0; }

    HMPDFPRINT(2, "\tcreate_ellgrid\n");

    SAFEALLOC(d->m->ellgrid, malloc((d->m->Nside/2+1) * sizeof(double)));
    SAFEHMPDF(linspace(d->m->Nside/2+1,
                       0.0, M_PI/d->f->pixelside,
                       d->m->ellgrid));

    d->m->created_ellgrid = 1;

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
    long w = (long)ceil(tout);
    ws->bufside = 2 * w + 1;
    long pixside = 2 * d->m->pxlgrid + 1;

    // draw random displacement of the center of the halo
    double dx = gsl_rng_uniform(ws->rng) - 0.5;
    double dy = gsl_rng_uniform(ws->rng) - 0.5;

    long Npix_filled = 0;
    while (Npix_filled < ws->bufside * ws->bufside)
    {
        long Npix_here
            = GSL_MIN(ws->bufside * ws->bufside - Npix_filled,           // physical constraint
                      (d->m->buflen - Npix_filled) / (pixside*pixside)); // memory constraint
        HMPDFCHECK(Npix_here <= 0, "no buffer left. this is a bug.");
        
        long posidx = 0;
        // fill the position buffer
        for (long ii=0; ii<Npix_here; ii++)
        {
            // figure out pixel coordinates in the map
            long xx = (ii + Npix_filled) / ws->bufside - w;
            long yy = (ii + Npix_filled) % ws->bufside - w;

            // loop over sample points within the pixel
            for (long xp= -d->m->pxlgrid; xp<= d->m->pxlgrid; xp++)
            {
                for (long yp= -d->m->pxlgrid; yp<= d->m->pxlgrid; yp++, posidx++)
                {
                    double xpos = (double)xx + (double)(2*xp)/(double)pixside + dx;
                    double ypos = (double)yy + (double)(2*yp)/(double)pixside + dy;
                    ws->pos[posidx] = hypot(xpos, ypos) / tout;
                }
            }
        }

        // evaluate the profile interpolator
        SAFEHMPDF(s_of_t(d, z_index, M_index, posidx, ws->pos, ws->buf+Npix_filled));

        posidx = 0;
        // perform the average
        for (long ii=0; ii<Npix_here; ii++)
        {
            double temp = 0.0;
            for (long jj=0; jj<pixside*pixside; jj++, posidx++)
            {
                temp += ws->buf[Npix_filled + posidx];
            }
            temp /= (double)(pixside*pixside);
            ws->buf[Npix_filled + ii] = temp;
        }

        Npix_filled +=  Npix_here;
    }

    ENDFCT
}//}}}

// convenience macro to reduce typing
#define INNERLOOP_OP                      \
    ws->map[ixx*ldmap + iyy]              \
        += ws->buf[xx * ws->bufside + yy];

static inline void 
add_buf_inner_loop(hmpdf_obj *d, map_ws *ws, long y0, long xx, long ixx)
{//{{{
    long ldmap = (ws->for_fft) ? (d->m->Nside+2) : d->m->Nside;

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
    //     with the periodic boundary conditions
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
         UNLIKELY(xx< ws->bufside);
         xx++, ixx++)
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

    // expected number of halos in this bin
    double n = d->h->hmf[z_index][M_index] // dn_3d / dlogM
               * gsl_pow_2(d->c->comoving[z_index])
               / d->c->hubble[z_index]
               * d->n->zweights[z_index]
               * d->n->Mweights[M_index]
               * d->m->area;

    if (d->m->mappoisson)
    {
        // use this funny construct because this function
        //    frequently sets errno (due to underflows I think)
        //    which is not critical but would be caught
        //    by ENDFCT
        SAFEGSL((*N = gsl_ran_poisson(ws->rng, n), 0));
    }
    else
    {
        double r = gsl_rng_uniform(ws->rng);
        double w_ceil = n - floor(n);
        *N = (unsigned)round((r < w_ceil) ? ceil(n) : floor(n));
    }

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

        HMPDFCHECK(ws->bufside >= d->m->Nside,
                   "attempting to add a halo that is larger than the map. "
                   "You should make the map larger.");

        for (unsigned ii=0; ii<N; ii++)
        {
            SAFEHMPDF(add_buf(d, ws));
        }
    }

    ENDFCT
}//}}}

static int
add_grf(hmpdf_obj *d, double (*pwr_spec)(double, void *), void *pwr_spec_params)
// adds random GRF realization to the Fourier space map
// if pwr_spec == NULL, no computation is performed
{//{{{
    STARTFCT

    if (pwr_spec == NULL)
    {
        return 0;
    }
    else
    {
        HMPDFPRINT(2, "\tadding Gaussian random field to the map\n");

        for (long ii=0; ii<d->m->Nside; ii++)
        // loop over the long direction (rows)
        {
            double ell1 = WAVENR(d->m->Nside, d->m->ellgrid, ii);

            for (long jj=0; jj<d->m->Nside/2+1; jj++)
            {
                double ell2 = WAVENR(d->m->Nside, d->m->ellgrid, jj);
                double ellmod = hypot(ell1, ell2);

                double Cl = pwr_spec(ellmod, pwr_spec_params);

                HMPDFCHECK(Cl < 0.0, "power spectrum must be positive everywhere.");

                double complex ampl
                    = gsl_ran_gaussian(d->m->ws[0]->rng, 1.0)
                      + _Complex_I * gsl_ran_gaussian(d->m->ws[0]->rng, 1.0);
                ampl *= sqrt(0.5 * Cl) / d->f->pixelside
                        * (double)(d->m->Nside);

                d->m->map_comp[ii*(d->m->Nside/2+1)+jj] += ampl;
            }
        }
    }

    ENDFCT
}//}}}

static int
filter_map(hmpdf_obj *d, double complex *map_comp, int *z_index)
{//{{{
    STARTFCT

    if (z_index == NULL)
    {
        HMPDFPRINT(3, "\t\tapplying filters to the map\n");
    }
    else
    {
        HMPDFPRINT(4, "\t\t\tapplying redshift dependent filters to the map\n");
    }

    double *ellmod;
    SAFEALLOC(ellmod, malloc((d->m->Nside/2+1) * sizeof(double)));

    for (long ii=0; ii<d->m->Nside; ii++)
    // loop over long direction (rows)
    {
        double ell1 = WAVENR(d->m->Nside, d->m->ellgrid, ii);

        for (long jj=0; jj<d->m->Nside/2+1; jj++)
        // loop over short direction (cols)
        {
            double ell2 = WAVENR(d->m->Nside, d->m->ellgrid, jj);
            ellmod[jj] = hypot(ell1, ell2);
        }

        SAFEHMPDF(apply_filters_map(d, d->m->Nside/2+1, ellmod,
                                    map_comp + ii * (d->m->Nside/2+1),
                                    map_comp + ii * (d->m->Nside/2+1),
                                    z_index));
    }

    free(ellmod);

    ENDFCT
}//}}}

static int
loop_no_z_dependence(hmpdf_obj *d)
// the loop if there are no z-dependent filters
// NOTE : the map that comes out of this is in conjugate space
//        if d->m->need_ft
{//{{{
    STARTFCT

    HMPDFPRINT(3, "\t\tloop_no_z_dependence\n");

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

        // TODO write status updates (similar to covariance.c)

        int z_index = bins[ii] / d->n->NM;
        int M_index = bins[ii] % d->n->NM;
        SAFEHMPDF_NORETURN(do_this_bin(d, z_index, M_index,
                                       d->m->ws[THIS_THREAD]));
    }

    free(bins);

    // add to the total map
    for (int ii=0; ii<d->m->Nws; ii++)
    {
        for (long jj=0; jj<d->m->Nside; jj++)
        {
            for (long kk=0; kk<d->m->Nside; kk++)
            {
                d->m->map_real[jj*d->m->ldmap + kk]
                    += d->m->ws[ii]->map[jj*d->m->Nside+kk];
            }
        }
    }

    if (d->m->need_ft)
    {
        // transform to conjugate space
        fftw_execute(*(d->m->p_r2c));
    }

    ENDFCT
}//}}}

static int
loop_w_z_dependence(hmpdf_obj *d)
// the less efficient loop we need to perform if there is a z-dependent filter
// NOTE : the map that comes out of this is in conjugate space!
{//{{{
    STARTFCT

    HMPDFPRINT(3, "\t\tloop_w_z_dependence\n");

    int *bins;
    SAFEALLOC(bins, malloc(d->n->NM * sizeof(int)));
    for (int ii=0; ii<d->n->NM; ii++)
    {
        bins[ii] = ii;
    }

    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        // reset the workspaces
        for (int ii=0; ii<d->m->Nws; ii++)
        {
            SAFEHMPDF(reset_map_ws(d, d->m->ws[ii]));
        }

        // shuffle to equalize load
        gsl_ran_shuffle(d->m->ws[0]->rng, bins, d->n->NM, sizeof(int));

        #ifdef _OPENMP
        #   pragma omp parallel for num_threads(d->m->Nws) schedule(dynamic)
        #endif
        for (int ii=0; ii<d->n->NM; ii++)
        {
            int M_index = bins[ii];
            SAFEHMPDF_NORETURN(do_this_bin(d, z_index, M_index,
                                           d->m->ws[THIS_THREAD]));
        }

        // sum all sub-maps in the 0th one (which always exists)
        for (int ii=1; ii<d->m->Nws; ii++)
        {
            for (long jj=0; jj<d->m->Nside; jj++)
            {
                for (long kk=0; kk<d->m->Nside; kk++)
                {
                    d->m->ws[0]->map[jj*(d->m->Nside+2) + kk]
                        += d->m->ws[ii]->map[jj*d->m->Nside + kk];
                }
            }
        }

        // transform to conjugate space
        fftw_execute(*(d->m->ws[0]->p_r2c));

        // apply the z-dependent filters
        SAFEHMPDF(filter_map(d, d->m->ws[0]->map_comp, &z_index));

        // add to the total map
        for (long ii=0; ii<d->m->Nside * (d->m->Nside/2+1); ii++)
        {
            d->m->map_comp[ii] += d->m->ws[0]->map_comp[ii];
        }
    }

    free(bins);

    ENDFCT
}//}}}

static int
create_mem(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->m->created_mem) { return 0; }

    HMPDFPRINT(2, "\tcreate_mem\n");

    if (d->f->Nfilters > 1 // the pixelization is done in real space,
                           //     which is more accurate
        || d->ns->have_noise)
    {
        d->m->ldmap = d->m->Nside + 2;
        d->m->need_ft = 1;
    }
    else
    {
        d->m->ldmap = d->m->Nside;
        d->m->need_ft = 0;
    }

    SAFEALLOC(d->m->map_real, ((d->m->need_ft) ?
                               fftw_malloc : malloc)(d->m->Nside * d->m->ldmap
                                                     * sizeof(double)));

    if (d->m->need_ft)
    {
        d->m->map_comp = (double complex *)d->m->map_real;

        SAFEALLOC(d->m->p_r2c, malloc(sizeof(fftw_plan)));
        *(d->m->p_r2c) = fftw_plan_dft_r2c_2d(d->m->Nside, d->m->Nside,
                                              d->m->map_real, d->m->map_comp,
                                              FFTW_MEASURE);

        SAFEALLOC(d->m->p_c2r, malloc(sizeof(fftw_plan)));
        *(d->m->p_c2r) = fftw_plan_dft_c2r_2d(d->m->Nside, d->m->Nside,
                                              d->m->map_comp, d->m->map_real,
                                              FFTW_MEASURE);
    }

    d->m->created_mem = 1;

    ENDFCT
}//}}}

static int
get_map_mean(hmpdf_obj *d)
{//{{{
    STARTFCT
    
    double sum = 0.0;
    for (long ii=0; ii<d->m->Nside; ii++)
    {
        for (long jj=0; jj<d->m->Nside; jj++)
        {
            sum += d->m->map_real[ii*d->m->ldmap+jj];
        }
    }
    
    d->m->mean = sum / (double)(d->m->Nside * d->m->Nside);

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

    HMPDFPRINT(2, "\tcreate_map\n");

    // zero the map
    zero_real(d->m->Nside * d->m->ldmap, d->m->map_real);

    // run the loop
    if (d->f->has_z_dependent)
    {
        SAFEHMPDF(loop_w_z_dependence(d));
    }
    else
    {
        SAFEHMPDF(loop_no_z_dependence(d));
    }

    if (d->m->need_ft)
    {
        // TODO in which order do the next two operations happen???
        // add the Gaussian random field
        SAFEHMPDF(add_grf(d, d->ns->noise_pwr, d->ns->noise_pwr_params));

        // apply the filters (not z-dependent)
        SAFEHMPDF(filter_map(d, d->m->map_comp, NULL));

        // transform back to real space
        fftw_execute(*(d->m->p_c2r));

        // normalize properly
        for (long ii=0; ii<d->m->Nside * (d->m->Nside+2); ii++)
        {
            d->m->map_real[ii] /= (double)(d->m->Nside * d->m->Nside);
        }
    }

    if (d->p->stype == hmpdf_kappa)
    {
        SAFEHMPDF(get_map_mean(d));
    }

    d->m->created_map = 1;

    ENDFCT
}//}}}

static int
create_sidelengths(hmpdf_obj *d)
{//{{{
    STARTFCT

    if (d->m->created_sidelengths) { return 0; }

    HMPDFPRINT(2, "\tcreate_sidelengths\n");

    double map_side = sqrt(d->m->area);
    d->m->Nside = (long)round(map_side/d->f->pixelside);

    double max_t_out = 0.0;
    for (int z_index=0; z_index<d->n->Nz; z_index++)
    {
        for (int M_index=0; M_index<d->n->NM; M_index++)
        {
            max_t_out = GSL_MAX(max_t_out,
                                d->p->profiles[z_index][M_index][0]);
        }
    }

    long temp = (long)round(max_t_out/d->f->pixelside);
    temp *= 2;
    temp += 4; // some safety buffer
    d->m->buflen = 2 * temp * temp; // not sufficient to do all halos
                                    // in one go, but long enough that it
                                    // is reasonably efficient

    d->m->created_sidelengths = 1;

    HMPDFPRINT(3, "\t\tmap = %ld x %ld <=> %g GB\n",
                  d->m->Nside, d->m->Nside,
                  1e-9 * (double)(d->m->Nside * d->m->Nside
                                  * sizeof(double)));
    HMPDFPRINT(3, "\t\tbuffer <=> %g GB\n",
                  1e-9 * (double)(d->m->buflen * sizeof(double)));

    ENDFCT
}//}}}

static int
prepare_maps(hmpdf_obj *d)
// compute map dimensions
// allocate workspaces
{//{{{
    STARTFCT

    HMPDFPRINT(1, "prepare_maps\n");

    SAFEHMPDF(create_sidelengths(d));
    SAFEHMPDF(create_mem(d));
    SAFEHMPDF(create_ellgrid(d));
    SAFEHMPDF(create_map_ws(d));
    SAFEHMPDF(create_map(d));

    ENDFCT
}//}}}

static int
common_input_processing(hmpdf_obj *d, int new_map)
{//{{{
    STARTFCT

    HMPDFCHECK(d->m->area < 0.0,
               "no/invalid sky fraction passed.");
    HMPDFCHECK(d->f->pixelside < 0.0,
               "no/invalid pixel sidelength passed.");

    if (new_map)
    {
        d->m->created_map = 0;
    }

    SAFEHMPDF(prepare_maps(d));

    ENDFCT
}//}}}

int
hmpdf_get_map_op(hmpdf_obj *d, int Nbins, double binedges[Nbins+1], double op[Nbins], int new_map)
// if (new_map), create one
// else, if not available, create one
//       else, use the existing one
{//{{{
    STARTFCT

    HMPDFCHECK(not_monotonic(Nbins+1, binedges, 1),
               "binedges not monotonically increasing.");

    SAFEHMPDF(common_input_processing(d, new_map));
    
    double _binedges[Nbins+1];
    SAFEHMPDF(pdf_adjust_binedges(d, Nbins, binedges,
                                  _binedges, d->m->mean));

    // prepare a histogram
    gsl_histogram *h;
    SAFEALLOC(h, gsl_histogram_alloc(Nbins));
    SAFEGSL(gsl_histogram_set_ranges(h, _binedges, Nbins+1));

    // accumulate the histogram
    for (long ii=0; ii<d->m->Nside; ii++)
    {
        for (long jj=0; jj<d->m->Nside; jj++)
        {
            double val = d->m->map_real[ii*d->m->ldmap+jj];
            if (val < _binedges[0] || val >= _binedges[Nbins])
            {
                continue;
            }
            else
            {
                SAFEGSL(gsl_histogram_increment(h, val));
            }
        }
    }

    // normalize
    SAFEGSL(gsl_histogram_scale(h, 1.0/(double)(d->m->Nside * d->m->Nside)));

    // write into output
    for (int ii=0; ii<Nbins; ii++)
    {
        op[ii] = gsl_histogram_get(h, ii);
    }

    gsl_histogram_free(h);

    ENDFCT
}//}}}

int
hmpdf_get_map(hmpdf_obj *d, double **map, long *Nside, int new_map)
{//{{{
    STARTFCT

    SAFEHMPDF(common_input_processing(d, new_map));

    SAFEALLOC(*map, malloc(d->m->Nside * d->m->Nside * sizeof(double)));
    for (long ii=0; ii<d->m->Nside; ii++)
    {
        memcpy(*map + ii*d->m->Nside,
               d->m->map_real + ii*d->m->ldmap,
               d->m->Nside * sizeof(double));
    }

    if (d->p->stype == hmpdf_kappa)
    {
        for (long ii=0; ii<d->m->Nside*d->m->Nside; ii++)
        {
            (*map)[ii] -= d->m->mean;
        }
    }

    *Nside = d->m->Nside;

    ENDFCT
}//}}}
