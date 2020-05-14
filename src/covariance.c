#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>

#include "utils.h"
#include "configs.h"
#include "power.h"
#include "onepoint.h"
#include "twopoint.h"
#include "covariance.h"

#include "hmpdf.h"

static int comp_int(const void *a, const void *b)
// to qsort an array of ints
{//{{{
    int int_a = *((int *)(a));
    int int_b = *((int *)(b));
    if (int_a == int_b)
    {
        return 0;
    }
    else if (int_a > int_b)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}//}}}

static
void create_phigrid(all_data *d)
{//{{{
    printf("\tcreate_phigrid\n");
    double *_phigrid = (double *)malloc(2 * d->n->gr->Nphi * sizeof(double));
    double *_phiweights = (double *)malloc(2 * d->n->gr->Nphi * sizeof(double));
    // allocate twice as much as we probably need to be safe,
    // the continuum integral adds a bit of noise

    // first treat the exact pixelization part
    int *_rsq = (int *)malloc(d->n->gr->pixelexactmax * d->n->gr->pixelexactmax
                              * sizeof(int));
    int N = 0; // counts possibly repeated sums of two squares
    // loop over one quadrant
    for (int ii=0; ii<=d->n->gr->pixelexactmax; ii++)
    {
        for (int jj=1; jj*jj<=d->n->gr->pixelexactmax * d->n->gr->pixelexactmax - ii*ii; jj++)
        {
            _rsq[N++] = ii*ii+jj*jj;
        }
    }
    // sort and identify unique values
    qsort(_rsq, N, sizeof(int), comp_int);
    int Nexact = 0; // counts number of phi-values for which exact treatment is needed
    for (int ii=0; ii<N;)
    {
        int ctr = 0;
        int jj;
        for (jj=ii; _rsq[jj] == _rsq[ii]; jj++)
        {
            ctr += 1;
        }
        if (Nexact >= d->n->gr->Nphi)
        {
            printf("Error : N_phi = %d too small.\n", d->n->gr->Nphi);
        }
        _phigrid[Nexact] = sqrt((double)(_rsq[ii]))
                                    * d->f->pixelside;
        _phiweights[Nexact++] = (double)(4 * ctr);
        ii = jj;
    }
    free(_rsq);

    // use a falling density of sample points
    // we transform phi = x^phipwr and sample uniformly in x
    // then the "density of states" is
    //      rho(phi) = rho0 * phi^(1/a - 1)
    // We want to fix the constant rho0 such that \int_phimin^phimax rho(phi) dphi = Nphi
    double integral = d->n->gr->phipwr
                      * (pow(d->n->gr->phimax, 1.0/d->n->gr->phipwr)
                         - pow(d->f->pixelside, 1.0/d->n->gr->phipwr));
    double rho0 = (double)(d->n->gr->Nphi) / integral;
    int end = Nexact; // index of the first free element in phigrid
    for (int ii=0; ii<Nexact; ii++)
    {
        double deltaphi;
        if (ii == 0)
        {
            deltaphi = _phigrid[1] - _phigrid[0];
        }
        else if (ii == Nexact-1)
        {
            deltaphi = _phigrid[Nexact-1] - _phigrid[Nexact-2];
        }
        else
        {
            deltaphi = 0.5 * (_phigrid[ii+1] - _phigrid[ii-1]);
        }
        double rho_here = pow(_phigrid[ii], 1.0/d->n->gr->phipwr-1.0);
        int Nphi_here = (int)(ceil(deltaphi * rho0 * rho_here)) - 1;
        if (Nphi_here%2)
        // want Nphi_here to be even for convenience
        // also this averages out the ceil
        {
            --Nphi_here;
        }
        _phiweights[ii] /= (double)(Nphi_here+1);
        for (int jj=1; jj<=Nphi_here/2; jj++)
        {
            _phigrid[end+2*jj-2] = _phigrid[ii]
                                   + deltaphi * d->n->gr->phijitter
                                     * (double)(jj) / (double)(Nphi_here/2);
            _phigrid[end+2*jj-1] = _phigrid[ii]
                                   - deltaphi * d->n->gr->phijitter
                                     * (double)(jj) / (double)(Nphi_here/2);
            _phiweights[end+2*jj-2] = _phiweights[end+2*jj-1] = _phiweights[ii];
        }
        end += Nphi_here;
    }
    Nexact = end; // Nphi - Nexact is the number of continuum phi values

    if (Nexact > d->n->gr->Nphi)
    {
        printf("Error : N_phi = %d too small "
               "(suggested increase to at least %d)\n", d->n->gr->Nphi, 2*Nexact);
    }

    // now fill the rest of the phi-grid
    // the continuum version is
    //      \int_lo^\phi_max d\phi 2\pi\phi/(pixel sidelength)^2 * integrand
    double lo = pow((double)(d->n->gr->pixelexactmax) * d->f->pixelside,
                    1.0/d->n->gr->phipwr);
    double hi = pow(d->n->gr->phimax, 1.0/d->n->gr->phipwr);
    gsl_integration_fixed_workspace *t =
        gsl_integration_fixed_alloc(gsl_integration_fixed_legendre, d->n->gr->Nphi-Nexact,
                                    lo, hi, 0.0, 0.0);
    double *nodes = gsl_integration_fixed_nodes(t);
    double *weights = gsl_integration_fixed_weights(t);
    int nn = Nexact;
    for (int ii=0; ii<gsl_integration_fixed_n(t); ii++)
    {
        if (nodes[ii] < lo)
        {
            continue;
        }
        else if (nodes[ii] > hi)
        // nodes are sorted
        {
            break;
        }
        else
        {
            _phigrid[nn] = pow(nodes[ii], d->n->gr->phipwr);
            // fix the normalization
            _phiweights[nn] = weights[ii]
                                       * 2.0 * M_PI * _phigrid[nn]
                                       / gsl_pow_2(d->f->pixelside)
                                       // Jacobian of the transformation
                                       * d->n->gr->phipwr
                                       * pow(nodes[ii], d->n->gr->phipwr-1.0);
            ++nn;
        }
    }
    d->n->gr->Nphi = nn;
    printf("\t\tNphi=%d, Nexact=%d\n", d->n->gr->Nphi, Nexact);
    gsl_integration_fixed_free(t);

    // now copy into the main grids
    // including shuffling to make parallel execution more efficient
    d->n->gr->phigrid = (double *)malloc(d->n->gr->Nphi * sizeof(double));
    d->n->gr->phiweights = (double *)malloc(d->n->gr->Nphi * sizeof(double));
    int *_indices = (int *)malloc(d->n->gr->Nphi * sizeof(int));
    for (int ii=0; ii<d->n->gr->Nphi; ii++)
    {
        _indices[ii] = ii;
    }
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r, 10); // TODO seed this randomly
    gsl_ran_shuffle(r, _indices, d->n->gr->Nphi, sizeof(int));
    gsl_rng_free(r);
    for (int ii=0; ii<d->n->gr->Nphi; ii++)
    {
        d->n->gr->phigrid[ii] = _phigrid[_indices[ii]];
        d->n->gr->phiweights[ii] = _phiweights[_indices[ii]];
    }
    free(_phigrid);
    free(_phiweights);
    free(_indices);

    // TODO TEST
    /*
    {
        FILE *f = fopen("phigrid.dat", "w");
        for (int ii=0; ii<d->n->gr->Nphi; ii++)
        {
            fprintf(f, "%.8e %.8e\n", d->n->gr->phigrid[ii]*180.0*60.0/M_PI,
                                      d->n->gr->phiweights[ii]);
        }
        fclose(f);
    }
    */
}//}}}

static
void alloc_tp_ws(all_data *d)
{//{{{
    printf("\talloc_tp_ws\n");
    printf("\t\tIn alloc_tp_ws : Trying to allocate workspaces for %d cores.\n", d->Ncores);
    d->cov->ws = (twopoint_workspace **)malloc(d->Ncores * sizeof(twopoint_workspace *));
    d->cov->Nws = 0;
    // allocate workspaces until we run out of memory
    for (int ii=0; ii<d->Ncores; ii++)
    {
        d->cov->ws[ii] = new_tp_ws(d->n->gr->Nsignal);
        if (d->cov->ws[ii] == NULL)
        {
            break;
        }
        else
        {
            ++d->cov->Nws;
        }
    }
    // NULL the remaining worspaces
    for (int ii=d->cov->Nws; ii<d->Ncores; ii++)
    {
        d->cov->ws[ii] = NULL;
    }
    printf("\t\tAllocated %d workspaces.\n", d->cov->Nws);
}//}}}

static
double corr_diagn(all_data *d, twopoint_workspace *ws)
{//{{{
    double out = 0.0;
    for (int ii=0; ii<d->n->gr->Nsignal; ii++)
    {
        for (int jj=0; jj<d->n->gr->Nsignal; jj++)
        {
            // assumes pdf_real to be properly normalized!
            out += ws->pdf_real[ii*(d->n->gr->Nsignal+2)+jj]
                   * (d->n->gr->signalgrid[ii] - d->op->signalmeanc)
                   * (d->n->gr->signalgrid[jj] - d->op->signalmeanc);
        }
    }
    return out;
}//}}}

static
void status_update(time_t t1, time_t t0, int done, int tot)
{//{{{
    double delta_time = difftime(t1, t0);
    double remains = delta_time/(double)done
                     *(double)(tot - done);
    int hrs = (int)floor(remains/60.0/60.0);
    int min = (int)round(remains/60.0 - 60.0*(double)(hrs));
    int done_perc = (int)round(100.0*(double)(done)/(double)(tot));
    printf("\t\t%3d %% done, %.2d hrs %.2d min remaining in create_cov.\n",
           done_perc, hrs, min);
}//}}}

static
void add_shotnoise(all_data *d)
// adds the zero-separation contribution to the covariance matrix
{//{{{
    for (int ii=0; ii<d->n->gr->Nsignal; ii++)
    {
        d->cov->Cov[ii*d->n->gr->Nsignal+ii] += d->op->PDFc[ii];
        for (int jj=0; jj<d->n->gr->Nsignal; jj++)
        {
            d->cov->Cov[ii*d->n->gr->Nsignal+jj]
                -= d->op->PDFc[ii] * d->op->PDFc[jj];
        }
    }
}//}}}

static
void rescale_to_fsky1(all_data *d)
// divides covariance matrix by number of pixels in the sky
{//{{{
    double Npixels = 4.0*M_PI/gsl_pow_2(d->f->pixelside);
    for (int ii=0; ii<d->n->gr->Nsignal*d->n->gr->Nsignal; ii++)
    {
        d->cov->Cov[ii] /= Npixels;
    }
}//}}}

static
void create_cov(all_data *d)
{//{{{
    printf("\tcreate_cov\n");
    // zero covariance
    zero_real(d->n->gr->Nsignal*d->n->gr->Nsignal, d->cov->Cov);

    // status
    int Nstatus = 0;
    time_t start_time = time(NULL);
    
    // loop over phi values
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(d->cov->Nws)
    #endif
    for (int pp=0; pp<d->n->gr->Nphi; pp++)
    {
        // create twopoint at this phi
        create_tp(d, d->n->gr->phigrid[pp], d->cov->ws[this_core()]);

        // compute the correlation function
        d->cov->corr_diagn[pp] = corr_diagn(d, d->cov->ws[this_core()]);
        
        // add to covariance
        for (int ii=0; ii<d->n->gr->Nsignal; ii++)
        {
            for (int jj=0; jj<d->n->gr->Nsignal; jj++)
            {
                #ifdef _OPENMP
                // make sure no two threads add to the same element simultaneously
                #pragma omp atomic
                #endif
                d->cov->Cov[ii*d->n->gr->Nsignal+jj]
                    += d->n->gr->phiweights[pp]
                       * (d->cov->ws[this_core()]->pdf_real[ii*(d->n->gr->Nsignal+2)+jj]
                          - d->op->PDFc[ii] * d->op->PDFc[jj]);
            }
        }

        // status update
        #ifdef _OPENMP
        #pragma omp critical(Status)
        #endif
        {
            ++Nstatus;
            if (Nstatus%COV_STATUS_PERIOD == 0)
            {
                time_t this_time = time(NULL);
                status_update(this_time, start_time, Nstatus, d->n->gr->Nphi);
            }
        }
    }

    add_shotnoise(d);
    rescale_to_fsky1(d);
}//}}}

static
void prepare_cov(all_data *d)
{//{{{
    printf("In covariance.h -> prepare_cov :\n");
    // run necessary code from other modules
    create_corr(d);
    create_phi_indep(d);
    create_op(d);
    
    // create phi grid
    create_phigrid(d);
    
    // allocate the workspaces
    alloc_tp_ws(d);
    
    // allocate the covariance matrix
    d->cov->Cov = (double *)malloc(d->n->gr->Nsignal*d->n->gr->Nsignal*sizeof(double));

    // allocate the diagnostic correlation function
    d->cov->corr_diagn = (double *)malloc(d->n->gr->Nphi * sizeof(double));

    // create covariance matrix
    create_cov(d);
}//}}}

static
void load_cov(all_data *d, char *fname)
{//{{{
    int Nlines;
    double **_x = fromfile(fname, &Nlines, 1);
    if (Nlines != d->n->gr->Nsignal*d->n->gr->Nsignal+1)
    {
        printf("In get_cov : cov matrix loaded from file %s "
               "not compatible with Nsignal. Aborting.\n", fname);
        return;
    }
    d->op->signalmeanc = _x[0][0];
    d->cov->Cov = _x[0]+1;
}//}}}

static
void save_cov(all_data *d, char *fname)
{//{{{
    FILE *f = fopen(fname, "wb");
    fwrite(&(d->op->signalmeanc), sizeof(double), 1, f);
    fwrite(d->cov->Cov, sizeof(double),
           d->n->gr->Nsignal*d->n->gr->Nsignal, f);
    fclose(f);
}//}}}

void get_cov(all_data *d, int Nbins, double *binedges, double *out, char *name)
{//{{{
    char covfile[512];
    char corrfile[512];
    if (Nbins == 0 && name == NULL)
    {
        printf("ERROR : both Nbins=0 and name=NULL,\nnothing to do in get_cov.\n");
        return;
    }

    if (name != NULL)
    {
        sprintf(covfile, "%s_cov.bin", name);
        sprintf(corrfile, "%s_corr.bin", name);
    }

    int to_compute = 1;
    if (name != NULL)
    {
        if (isfile(covfile))
        // covariance matrix has already been computed, we just need binning
        {
            if (Nbins == 0)
            {
                printf("ERROR : covariance matrix named %s already exists,\n"
                       "and you requested no binning.\n"
                       "Nothing to do in get_cov.\n", covfile);
                return;
            }
            load_cov(d, covfile);
            to_compute = 0;
        }
        else
        {
            printf("\t\tIn get_cov : file %s not found. Will compute.\n", covfile);
        }
    }

    // perform the computation if necessary
    if (to_compute)
    {
        prepare_cov(d);
    }

    if ((name != NULL) && (to_compute))
    {
        save_cov(d, covfile);
        tofile(corrfile, d->n->gr->Nphi, 3, d->n->gr->phigrid,
               d->n->gr->phiweights, d->cov->corr_diagn);
    }

    if (Nbins > 0)
    // binning requested
    {
        // if kappa, adjust the bins
        double _binedges[Nbins+1];
        memcpy(_binedges, binedges, (Nbins+1) * sizeof(double));
        if (d->p->stype == kappa)
        {
            for (int ii=0; ii<=Nbins; ii++)
            {
                _binedges[ii] += d->op->signalmeanc;
            }
        }

        // perform the binning
        printf("\t\tbinning the covariance matrix\n");
        bin_2d(d->n->gr->Nsignal, d->n->gr->signalgrid, d->cov->Cov, COVINTEGR_N,
               Nbins, _binedges, out, TPINTERP_TYPE);
    }

}//}}}
