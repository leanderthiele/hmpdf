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
#include "data.h"
#include "power.h"
#include "profiles.h"
#include "onepoint.h"
#include "twopoint.h"
#include "covariance.h"

#include "hmpdf.h"

void null_covariance(all_data *d)
{//{{{{
    d->cov->ws = NULL;
    d->cov->Cov = NULL;
    d->cov->corr_diagn = NULL;
    d->cov->created_tp_ws = 0;
    d->cov->created_phigrid = 0;
    d->cov->created_cov = 0;
}//}}}

void reset_covariance(all_data *d)
{//{{{
    if (d->cov->Cov != NULL) { free(d->cov->Cov); }
    if (d->cov->corr_diagn != NULL) { free(d->cov->corr_diagn); }
    if (d->cov->ws != NULL)
    {
        for (int ii=0; ii<d->cov->Nws; ii++)
        {
            if (d->cov->ws[ii] != NULL)
            {
                fftw_free(d->cov->ws[ii]->pdf_real);
                free(d->cov->ws[ii]->bc);
                fftw_free(d->cov->ws[ii]->tempc_real);
                fftw_destroy_plan(d->cov->ws[ii]->pu_r2c);
                fftw_destroy_plan(d->cov->ws[ii]->pc_r2c);
                fftw_destroy_plan(d->cov->ws[ii]->ppdf_c2r);
                free(d->cov->ws[ii]);
            }
        }
        free(d->cov->ws);
    }
}//}}}

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
    if (d->cov->created_phigrid) { return; }
    fprintf(stdout, "\tcreate_phigrid\n");
    fflush(stdout);
    double *_phigrid = (double *)malloc(2 * d->n->Nphi * sizeof(double));
    double *_phiweights = (double *)malloc(2 * d->n->Nphi * sizeof(double));
    // allocate twice as much as we probably need to be safe,
    // the continuum integral adds a bit of noise

    // first treat the exact pixelization part
    int *_rsq = (int *)malloc(d->n->pixelexactmax * d->n->pixelexactmax
                              * sizeof(int));
    int N = 0; // counts possibly repeated sums of two squares
    // loop over one quadrant
    for (int ii=0; ii<=d->n->pixelexactmax; ii++)
    {
        for (int jj=1; jj*jj<=d->n->pixelexactmax * d->n->pixelexactmax - ii*ii; jj++)
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
        if (Nexact >= d->n->Nphi)
        {
            fprintf(stdout, "Error : N_phi = %d too small.\n", d->n->Nphi);
            fflush(stdout);
            return;
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
    double integral = d->n->phipwr
                      * (pow(d->n->phimax, 1.0/d->n->phipwr)
                         - pow(d->f->pixelside, 1.0/d->n->phipwr));
    double rho0 = (double)(d->n->Nphi) / integral;
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
        double rho_here = pow(_phigrid[ii], 1.0/d->n->phipwr-1.0);
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
                                   + deltaphi * d->n->phijitter
                                     * (double)(jj) / (double)(Nphi_here/2);
            _phigrid[end+2*jj-1] = _phigrid[ii]
                                   - deltaphi * d->n->phijitter
                                     * (double)(jj) / (double)(Nphi_here/2);
            _phiweights[end+2*jj-2] = _phiweights[end+2*jj-1] = _phiweights[ii];
        }
        end += Nphi_here;
    }
    Nexact = end; // Nphi - Nexact is the number of continuum phi values

    if (Nexact > d->n->Nphi)
    {
        fprintf(stdout, "Error : N_phi = %d too small "
                        "(suggested increase to at least %d)\n",
                        d->n->Nphi, 2*Nexact);
        fflush(stdout);
    }

    // now fill the rest of the phi-grid
    // the continuum version is
    //      \int_lo^\phi_max d\phi 2\pi\phi/(pixel sidelength)^2 * integrand
    double lo = pow((double)(d->n->pixelexactmax) * d->f->pixelside,
                    1.0/d->n->phipwr);
    double hi = pow(d->n->phimax, 1.0/d->n->phipwr);
    gsl_integration_fixed_workspace *t =
        gsl_integration_fixed_alloc(gsl_integration_fixed_legendre, d->n->Nphi-Nexact,
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
            _phigrid[nn] = pow(nodes[ii], d->n->phipwr);
            // fix the normalization
            _phiweights[nn] = weights[ii]
                                       * 2.0 * M_PI * _phigrid[nn]
                                       / gsl_pow_2(d->f->pixelside)
                                       // Jacobian of the transformation
                                       * d->n->phipwr
                                       * pow(nodes[ii], d->n->phipwr-1.0);
            ++nn;
        }
    }
    d->n->Nphi = nn;
    fprintf(stdout, "\t\tNphi=%d, Nexact=%d\n", d->n->Nphi, Nexact);
    fflush(stdout);
    gsl_integration_fixed_free(t);

    // now copy into the main grids
    // including shuffling to make parallel execution more efficient
    d->n->phigrid = (double *)malloc(d->n->Nphi * sizeof(double));
    d->n->phiweights = (double *)malloc(d->n->Nphi * sizeof(double));
    int *_indices = (int *)malloc(d->n->Nphi * sizeof(int));
    for (int ii=0; ii<d->n->Nphi; ii++)
    {
        _indices[ii] = ii;
    }
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r, 10); // TODO seed this randomly
    gsl_ran_shuffle(r, _indices, d->n->Nphi, sizeof(int));
    gsl_rng_free(r);
    for (int ii=0; ii<d->n->Nphi; ii++)
    {
        d->n->phigrid[ii] = _phigrid[_indices[ii]];
        d->n->phiweights[ii] = _phiweights[_indices[ii]];
    }
    free(_phigrid);
    free(_phiweights);
    free(_indices);

    d->cov->created_phigrid = 1;
}//}}}

static
void create_tp_ws(all_data *d)
{//{{{
    if (d->cov->created_tp_ws) { return; }
    fprintf(stdout, "\tcreate_tp_ws\n");
    fflush(stdout);
    fprintf(stdout, "\t\tIn create_tp_ws : Trying to allocate workspaces for %d cores.\n", d->Ncores);
    fflush(stdout);
    d->cov->ws = (twopoint_workspace **)malloc(d->Ncores * sizeof(twopoint_workspace *));
    d->cov->Nws = 0;
    // allocate workspaces until we run out of memory
    for (int ii=0; ii<d->Ncores; ii++)
    {
        d->cov->ws[ii] = new_tp_ws(d->n->Nsignal);
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
    fprintf(stdout, "\t\tAllocated %d workspaces.\n", d->cov->Nws);
    fflush(stdout);

    d->cov->created_tp_ws = 1;
}//}}}

static
double corr_diagn(all_data *d, twopoint_workspace *ws)
{//{{{
    double out = 0.0;
    for (int ii=0; ii<d->n->Nsignal; ii++)
    {
        for (int jj=0; jj<d->n->Nsignal; jj++)
        {
            // assumes pdf_real to be properly normalized!
            out += ws->pdf_real[ii*(d->n->Nsignal+2)+jj]
                   * (d->n->signalgrid[ii] - d->op->signalmeanc)
                   * (d->n->signalgrid[jj] - d->op->signalmeanc);
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
    fprintf(stdout, "\t\t%3d %% done, %.2d hrs %.2d min remaining "
                    "in create_cov.\n", done_perc, hrs, min);
    fflush(stdout);
}//}}}

static
void add_shotnoise_offiag(int N, double *cov, double *p)
// adds the zero-separation contribution to the covariance matrix
{//{{{
    for (int ii=0; ii<N; ii++)
    {
        for (int jj=0; jj<N; jj++)
        {
            cov[ii*N + jj] -= p[ii] * p[jj];
        }
    }
}//}}}

static
void add_shotnoise_diag(int N, double *cov, double *p)
{//{{{
    for (int ii=0; ii<N; ii++)
    {
        cov[ii*(N+1)] += p[ii];
    }
}//}}}

static
void rescale_to_fsky1(all_data *d, int N, double *cov)
// divides covariance matrix by number of pixels in the sky
{//{{{
    double Npixels = 4.0*M_PI/gsl_pow_2(d->f->pixelside);
    for (int ii=0; ii<N*N; ii++)
    {
        cov[ii] /= Npixels;
    }
}//}}}

static
void create_cov(all_data *d)
{//{{{
    if (d->cov->created_cov) { return; }
    fprintf(stdout, "\tcreate_cov\n");
    fflush(stdout);
    // zero covariance
    zero_real(d->n->Nsignal*d->n->Nsignal, d->cov->Cov);

    // status
    int Nstatus = 0;
    time_t start_time = time(NULL);
    
    // loop over phi values
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(d->cov->Nws)
    #endif
    for (int pp=0; pp<d->n->Nphi; pp++)
    {
        // create twopoint at this phi
        create_tp(d, d->n->phigrid[pp], d->cov->ws[this_core()]);

        // compute the correlation function
        d->cov->corr_diagn[pp] = corr_diagn(d, d->cov->ws[this_core()]);
        
        // add to covariance
        for (int ii=0; ii<d->n->Nsignal; ii++)
        {
            for (int jj=0; jj<d->n->Nsignal; jj++)
            {
                #ifdef _OPENMP
                // make sure no two threads add to the same element simultaneously
                #pragma omp atomic
                #endif
                d->cov->Cov[ii*d->n->Nsignal+jj]
                    += d->n->phiweights[pp]
                       * (d->cov->ws[this_core()]->pdf_real[ii*(d->n->Nsignal+2)+jj]
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
                status_update(this_time, start_time, Nstatus, d->n->Nphi);
            }
        }
    }

    d->cov->created_cov = 1;
}//}}}

static
void prepare_cov(all_data *d)
{//{{{
    fprintf(stdout, "In covariance.h -> prepare_cov :\n");
    fflush(stdout);
    // run necessary code from other modules
    create_corr(d);
    if (d->f->Nfilters > 0)
    {
        create_conj_profiles(d);
        create_filtered_profiles(d);
    }
    create_breakpoints_or_monotonize(d);
    create_phi_indep(d);
    create_op(d);
    
    // create phi grid
    create_phigrid(d);
    
    // allocate the workspaces
    create_tp_ws(d);
    
    // create covariance matrix
    create_cov(d);
}//}}}

static
void load_cov(all_data *d, char *fname)
{//{{{
    int Nlines;
    double **_x = fromfile(fname, &Nlines, 1);
    if (Nlines != d->n->Nsignal*d->n->Nsignal+1)
    {
        fprintf(stdout, "In get_cov : cov matrix loaded from file %s "
                        "not compatible with Nsignal. Aborting.\n", fname);
        fflush(stdout);
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
           d->n->Nsignal*d->n->Nsignal, f);
    fclose(f);
}//}}}

void get_cov(all_data *d, int Nbins, double *binedges, double *out, int noisy, char *name)
{//{{{
    if (noisy && d->op->noise<0.0)
    {
        fprintf(stdout, "Error: noisy cov-matrix requested but no/invalid noise level passed.\n");
        fflush(stdout);
        return;
    }

    char covfile[512];
    char corrfile[512];
    if (Nbins == 0 && name == NULL)
    {
        fprintf(stdout, "ERROR : both Nbins=0 and name=NULL,\nnothing to do in get_cov.\n");
        fflush(stdout);
        return;
    }

    if (name != NULL)
    {
        sprintf(covfile, "%s_cov.bin", name);
        sprintf(corrfile, "%s_corr.bin", name);
    }

    // allocate the covariance matrix (if not allocated yet)
    if (d->cov->Cov == NULL)
    {
        d->cov->Cov = (double *)malloc(d->n->Nsignal*d->n->Nsignal*sizeof(double));
    }

    // allocate the diagnostic correlation function (if not allocated yet)
    if (d->cov->corr_diagn == NULL)
    {
        d->cov->corr_diagn = (double *)malloc(d->n->Nphi * sizeof(double));
    }

    int to_compute = 1;
    if (name != NULL)
    {
        if (isfile(covfile))
        // covariance matrix has already been computed, we just need binning
        {
            if (Nbins == 0)
            {
                fprintf(stdout, "ERROR : covariance matrix named %s already exists,\n"
                                "and you requested no binning.\n"
                                "Nothing to do in get_cov.\n", covfile);
                fflush(stdout);
                return;
            }
            load_cov(d, covfile);
            to_compute = 0;
        }
        else
        {
            fprintf(stdout, "\t\tIn get_cov : file %s not found. Will compute.\n", covfile);
            fflush(stdout);
        }
    }

    // perform the computation if necessary
    if (to_compute)
    {
        prepare_cov(d);
    }

    // save the raw data if necessary
    if ((name != NULL) && (to_compute))
    {
        save_cov(d, covfile);
        tofile(corrfile, d->n->Nphi, 3, d->n->phigrid,
               d->n->phiweights, d->cov->corr_diagn);
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
                _binedges[ii] += (noisy) ? d->op->signalmeanc_noisy
                                 : d->op->signalmeanc;
            }
        }

        // prepare final covariance matrix
        int N = (noisy) ? d->n->Nsignal_noisy : d->n->Nsignal;
        double *final_cov = (double *)malloc(N * N * sizeof(double));
        if (noisy)
        {
            zero_real(N * N, final_cov);
            int len_kernel = (d->n->Nsignal_noisy - d->n->Nsignal)/2;
            for (int ii=0; ii<d->n->Nsignal; ii++)
            {
                memcpy(final_cov + (len_kernel+ii)*d->n->Nsignal_noisy + len_kernel,
                       d->cov->Cov + ii*d->n->Nsignal, d->n->Nsignal * sizeof(double));
            }
        }
        else
        {
            memcpy(final_cov, d->cov->Cov, N * N * sizeof(double));
        }
        
        // add the off-diagonal shot-noise elements
        add_shotnoise_offiag((noisy) ? d->n->Nsignal_noisy : d->n->Nsignal,
                             final_cov,
                             (noisy) ? d->op->PDFc_noisy : d->op->PDFc);

        // perform the binning
        fprintf(stdout, "\t\tbinning the covariance matrix\n");
        fflush(stdout);
        bin_2d(N, (noisy) ? d->n->signalgrid_noisy : d->n->signalgrid,
               final_cov, COVINTEGR_N, Nbins, _binedges, out, TPINTERP_TYPE);

        // compute the shot noise term
        double *temp = (double *)malloc(Nbins * sizeof(double));
        get_op(d, Nbins, binedges, temp, cl, noisy);
        // normalize properly
        /*
        for (int ii=0; ii<Nbins; ii++)
        {
            temp[ii] *= binedges[ii+1] - binedges[ii];
        }
        */
        add_shotnoise_diag(Nbins, out, temp);
        free(temp);

        // normalize properly
        rescale_to_fsky1(d, Nbins, out);

        free(final_cov);
    }

}//}}}
