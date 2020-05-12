#include <stdbool.h>

#include <class.h>
#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_dht.h>

#include "utils.h"
#include "configs.h"
#include "class_interface.h"
#include "data.h"

#include "hmpdf.h"

void null_data(all_data *d)
{//{{{
    d->n->inited_numerics = 0;
    d->n->gr->zgrid = NULL;
    d->n->gr->zweights = NULL;
    d->n->gr->Mgrid = NULL;
    d->n->gr->Mweights = NULL;
    d->n->gr->signalgrid = NULL;
    d->n->gr->lambdagrid = NULL;
    d->n->gr->phigrid = NULL;
    d->n->gr->phiweights = NULL;

    /*
    d->cls->pr = NULL;
    d->cls->ba = NULL;
    d->cls->th = NULL;
    d->cls->pt = NULL;
    d->cls->tr = NULL;
    d->cls->pm = NULL;
    d->cls->sp = NULL;
    d->cls->nl = NULL;
    d->cls->le = NULL;
    d->cls->op = NULL;
    */
    d->cls = NULL;

    d->c->inited_cosmo = 0;
    d->c->hubble = NULL;
    d->c->comoving = NULL;
    d->c->angular_diameter = NULL;
    d->c->Scrit = NULL;
    d->c->Dsq = NULL;
    d->c->rho_m = NULL;
    d->c->rho_c = NULL;
    d->c->Om = NULL;

    d->pwr->inited_power = 0;
    d->pwr->k_arr = NULL;
    d->pwr->Pk_arr = NULL;
    d->pwr->Pk_interp = NULL;
    d->pwr->ssq = NULL;
    d->pwr->created_corr = 0;
    d->pwr->corr_interp = NULL;
    d->pwr->corr_accel = NULL;

    d->h->inited_halo = 0;
    d->h->hmf = NULL;
    d->h->bias = NULL;
    d->h->c_interp = NULL;
    d->h->c_accel = NULL;

    d->p->inited_profiles = 0;
    d->p->decr_tgrid = NULL;
    d->p->incr_tgrid = NULL;
    d->p->decr_tsqgrid = NULL;
    d->p->prtilde_thetagrid = NULL;
    d->p->reci_tgrid = NULL;
    d->p->breakpoints = NULL;
    d->p->dht_ws = NULL;
    d->p->profiles = NULL;
    d->p->created_conj_profiles = 0;
    d->p->conj_profiles = NULL;
    d->p->incr_tgrid_accel = NULL;
    d->p->reci_tgrid_accel = NULL;

    d->f->inited_filters = 0;
    d->f->ffilters = NULL;
    d->f->quadraticpixel_interp = NULL;
    d->f->quadraticpixel_accel = NULL;
    d->f->quadraticpixel_ellmin = NULL;
    d->f->quadraticpixel_ellmax = NULL;

    d->op->created_op = 0;
    d->op->PDFu = NULL;
    d->op->PDFc = NULL;

    d->tp->created_phi_indep = 0;
    d->tp->dtsq = NULL;
    d->tp->t = NULL;
    d->tp->ac = NULL;
    d->tp->au = NULL;
    d->tp->ws = NULL;

    /*
    d->ps->created_Cell = 0;
    d->ps->dht_ws = NULL;
    d->ps->ell = NULL;
    d->ps->Cell_1h = NULL;
    d->ps->Cell_2h = NULL;
    d->ps->Cell_tot = NULL;
    d->ps->created_Cphi = 0;
    d->ps->phi = NULL;
    d->ps->Cphi = NULL;
    d->ps->Cphi_interp = NULL;
    d->ps->Cphi_accel = NULL;
    */

    d->cov->ws = NULL;
    d->cov->Cov = NULL;
}//}}}

all_data *new_data(void)
{//{{{
    printf("In data.h -> new_data.\n");
    all_data *d = malloc(sizeof(all_data));
    d->n = malloc(sizeof(struct numerics));
//    d->cls = malloc(sizeof(struct class_interface));
    d->c = malloc(sizeof(struct cosmo));
    d->pwr = malloc(sizeof(struct power));
    d->h = malloc(sizeof(struct halo));
    d->f = malloc(sizeof(struct filters));
    d->p = malloc(sizeof(struct profiles));
    d->op = malloc(sizeof(struct onepoint));
    d->tp = malloc(sizeof(struct twopoint));
//    d->ps = malloc(sizeof(powerspectrum));
    d->cov = malloc(sizeof(struct covariance));

    d->n->gr = malloc(sizeof(struct grids));
    null_data(d);
    return d;
}//}}}

void reset_data(all_data *d)
{//{{{
    printf("In data.h -> reset_data.\n");
    if (d->n->gr->zgrid != NULL) { free(d->n->gr->zgrid); }
    if (d->n->gr->zweights != NULL) { free(d->n->gr->zweights); }
    if (d->n->gr->Mgrid != NULL) { free(d->n->gr->Mgrid); }
    if (d->n->gr->Mweights != NULL) { free(d->n->gr->Mweights); }
    if (d->n->gr->signalgrid != NULL) { free(d->n->gr->signalgrid); }
    if (d->n->gr->lambdagrid != NULL) { free(d->n->gr->lambdagrid); }
    if (d->n->gr->phigrid != NULL) { free(d->n->gr->phigrid); }
    if (d->n->gr->phiweights != NULL) { free(d->n->gr->phiweights); }

    /*
    if (d->cls->op != NULL) { free(d->cls->op); }
    if (d->cls->le != NULL) { lensing_free(d->cls->le); free(d->cls->le); }
    if (d->cls->sp != NULL) { spectra_free(d->cls->sp); free(d->cls->sp); }
    if (d->cls->tr != NULL) { transfer_free(d->cls->tr); free(d->cls->tr); }
    if (d->cls->nl != NULL) { nonlinear_free(d->cls->nl); free(d->cls->nl); }
    if (d->cls->pt != NULL) { perturb_free(d->cls->pt); free(d->cls->pt); }
    if (d->cls->pm != NULL) { primordial_free(d->cls->pm); free(d->cls->pm); }
    if (d->cls->th != NULL) { thermodynamics_free(d->cls->th); free(d->cls->th); }
    if (d->cls->ba != NULL) { background_free(d->cls->ba); free(d->cls->ba); }
    if (d->cls->pr != NULL) { free(d->cls->pr); }
    */

    if (d->c->hubble != NULL) { free(d->c->hubble); }
    if (d->c->comoving != NULL) { free(d->c->comoving); }
    if (d->c->angular_diameter != NULL) { free(d->c->angular_diameter); }
    if (d->c->Scrit != NULL) { free(d->c->Scrit); }
    if (d->c->Dsq != NULL) { free(d->c->Dsq); }
    if (d->c->rho_m != NULL) { free(d->c->rho_m); }
    if (d->c->rho_c != NULL) { free(d->c->rho_c); }
    if (d->c->Om != NULL) { free(d->c->Om); }

    if (d->pwr->k_arr != NULL) { free(d->pwr->k_arr); }
    if (d->pwr->Pk_arr != NULL) { free(d->pwr->Pk_arr); }
    if (d->pwr->Pk_interp != NULL) { delete_interp1d(d->pwr->Pk_interp); }
    if (d->pwr->ssq != NULL)
    {
        for (int M_index=0; M_index<d->n->gr->NM; M_index++)
        {
            if (d->pwr->ssq[M_index] != NULL)
            {
                free(d->pwr->ssq[M_index]);
            }
        }
        free(d->pwr->ssq);
    }
    if (d->pwr->corr_interp != NULL) { gsl_spline_free(d->pwr->corr_interp); }
    if (d->pwr->corr_accel != NULL)
    {
        for (int ii=0; ii<d->pwr->Ncorr_accel; ii++)
        {
            gsl_interp_accel_free(d->pwr->corr_accel[ii]);
        }
        free(d->pwr->corr_accel);
    }

    if (d->h->hmf != NULL) { free(d->h->hmf); }
    if (d->h->bias != NULL) { free(d->h->bias); }
    if (d->h->c_interp != NULL) { gsl_spline_free(d->h->c_interp); }
    if (d->h->c_accel != NULL) { gsl_interp_accel_free(d->h->c_accel); }

    if (d->p->decr_tgrid != NULL) { free(d->p->decr_tgrid); }
    if (d->p->incr_tgrid != NULL) { free(d->p->incr_tgrid); }
    if (d->p->decr_tsqgrid != NULL) { free(d->p->decr_tsqgrid); }
    if (d->p->reci_tgrid != NULL) { free(d->p->reci_tgrid); }
    if (d->p->prtilde_thetagrid != NULL) { free(d->p->prtilde_thetagrid); }
    if (d->p->incr_tgrid_accel != NULL) { gsl_interp_accel_free(d->p->incr_tgrid_accel); }
    if (d->p->reci_tgrid_accel != NULL) { gsl_interp_accel_free(d->p->reci_tgrid_accel); }
    if (d->p->profiles != NULL)
    {
        for (int z_index=0; z_index<d->n->gr->Nz; z_index++)
        {
            if (d->p->profiles[z_index] != NULL)
            {
                for (int M_index=0; M_index<d->n->gr->NM; M_index++)
                {
                    if (d->p->profiles[z_index][M_index] != NULL)
                    {
                        free(d->p->profiles[z_index][M_index]);
                    }
                }
                free(d->p->profiles[z_index]);
            }
        }
        free(d->p->profiles);
    }
    if (d->p->conj_profiles != NULL)
    {
        for (int z_index=0; z_index<d->n->gr->Nz; z_index++)
        {
            if (d->p->conj_profiles[z_index] != NULL)
            {
                for (int M_index=0; M_index<d->n->gr->NM; M_index++)
                {
                    if (d->p->conj_profiles[z_index][M_index] != NULL)
                    {
                        free(d->p->conj_profiles[z_index][M_index]);
                    }
                }
                free(d->p->conj_profiles[z_index]);
            }
        }
        free(d->p->conj_profiles);
    }
    if (d->p->dht_ws != NULL) { gsl_dht_free(d->p->dht_ws); }

    if (d->f->ffilters != NULL)
    {
        for (int ii=0; ii<d->f->Nfilters; ii++)
        {
            free(d->f->ffilters[ii].params);
        }
        free(d->f->ffilters);
    }
    if (d->f->quadraticpixel_interp != NULL)
    {
        if (d->f->quadraticpixel_interp[0] != NULL) { gsl_spline_free(d->f->quadraticpixel_interp[0]); }
        if (d->f->quadraticpixel_interp[1] != NULL) { gsl_spline_free(d->f->quadraticpixel_interp[1]); }
        free(d->f->quadraticpixel_interp);
    }
    if (d->f->quadraticpixel_accel != NULL)
    {
        if (d->f->quadraticpixel_accel[0] != NULL) { gsl_interp_accel_free(d->f->quadraticpixel_accel[0]); }
        if (d->f->quadraticpixel_accel[1] != NULL) { gsl_interp_accel_free(d->f->quadraticpixel_accel[1]); }
    }
    if (d->f->quadraticpixel_ellmin != NULL) { free(d->f->quadraticpixel_ellmin); }
    if (d->f->quadraticpixel_ellmax != NULL) { free(d->f->quadraticpixel_ellmax); }
    
    if (d->op->PDFu != NULL) { free(d->op->PDFu); }
    if (d->op->PDFc != NULL) { free(d->op->PDFc); }

    if (d->tp->dtsq != NULL)
    {
        for (int z_index=0; z_index<d->n->gr->Nz; z_index++)
        {
            if (d->tp->dtsq[z_index] != NULL)
            {
                for (int M_index=0; M_index<d->n->gr->NM; M_index++)
                {
                    if (d->tp->dtsq[z_index][M_index] != NULL)
                    {
                        free(d->tp->dtsq[z_index][M_index]);
                    }
                }
                free(d->tp->dtsq[z_index]);
            }
        }
        free(d->tp->dtsq);
    }
    if (d->tp->t != NULL)
    {
        for (int z_index=0; z_index<d->n->gr->Nz; z_index++)
        {
            if (d->tp->t[z_index] != NULL)
            {
                for (int M_index=0; M_index<d->n->gr->NM; M_index++)
                {
                    if (d->tp->t[z_index][M_index] != NULL)
                    {
                        free(d->tp->t[z_index][M_index]);
                    }
                }
                free(d->tp->t[z_index]);
            }
        }
        free(d->tp->t);
    }
    if (d->tp->ac != NULL)
    {
        for (int z_index=0; z_index<d->n->gr->Nz; z_index++)
        {
            if (d->tp->ac[z_index] != NULL)
            {
                free(d->tp->ac[z_index]);
            }
        }
        free(d->tp->ac);
    }
    if (d->tp->au != NULL) { fftw_free(d->tp->au); }
    if (d->tp->ws != NULL)
    {
        fftw_free(d->tp->ws->pdf_real);
        free(d->tp->ws->bc);
        fftw_free(d->tp->ws->tempc_real);
        fftw_destroy_plan(d->tp->ws->pu_r2c);
        fftw_destroy_plan(d->tp->ws->pc_r2c);
        fftw_destroy_plan(d->tp->ws->ppdf_c2r);
        free(d->tp->ws);
    }

    /*
    if (d->ps->dht_ws != NULL) { gsl_dht_free(d->ps->dht_ws); }
    if (d->ps->ell != NULL) { free(d->ps->ell); }
    if (d->ps->Cell_1h != NULL) { free(d->ps->Cell_1h); }
    if (d->ps->Cell_2h != NULL) { free(d->ps->Cell_2h); }
    if (d->ps->Cell_tot != NULL) { free(d->ps->Cell_tot); }
    if (d->ps->phi != NULL) { free(d->ps->phi); }
    if (d->ps->Cphi != NULL) { free(d->ps->Cphi); }
    if (d->ps->Cphi_interp != NULL) { gsl_interp_free(d->ps->Cphi_interp); }
    if (d->ps->Cphi_accel != NULL) { gsl_interp_accel_free(d->ps->Cphi_accel); }
    */

    if (d->cov->Cov != NULL) { free(d->cov->Cov); }
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

    null_data(d);
}//}}}

void delete_data(all_data *d)
{//{{{
    printf("In data.h -> delete_data.\n");
    reset_data(d);

    free_class(d);
    free(d->c);
    free(d->pwr);
    free(d->f);
    free(d->p);
    free(d->h);
    free(d->n->gr);
    free(d->n);
    free(d->op);

    free(d);
}//}}}

