#include <stdio.h>
#include <stdarg.h>
#include <math.h>

#include "utils.h"
#include "configs.h"
#include "data.h"
#include "class_interface.h"
#include "cosmology.h"
#include "numerics.h"
#include "power.h"
#include "halo_model.h"
#include "filter.h"
#include "profiles.h"
#include "init.h"

#include "hmpdf.h"

typedef enum
{//{{{
    str_type,
    int_type,
    dbl_type,
    dptr_type,
    mdef_type,
    integr_type,
    vptr_type, // void *
    lf_type, // ell filter
    kf_type, // k filter
}//}}}
dtype;

typedef struct
{//{{{
    configs indx;
    void *target;
    dtype dt;
    int depends_on_stype;
    void *def[2];
    int set;
}//}}}
param;

void init(all_data *d, char *class_ini, signaltype stype, ...)
{//{{{
    printf("In init.h -> init.\n");
    // this frees all the computed quantities,
    // since we assume that each call of init changes some
    // parameter (cosmological or numerical)
    reset_data(d);

    d->cls->class_ini = class_ini;
    d->p->stype = stype;

    param p[] = {//{{{
        //indx       target         dt       depst default { kappa, tsz }
        {hmpdf_N_cores,
            &(d->Ncores),           int_type,     0, {&(def.Ncores), }},
        {hmpdf_class_pre,
            &(d->cls->class_pre),        str_type,     0, {&(def.class_pre), }},
        {hmpdf_N_z,
            &(d->n->Nz),        int_type,     0, {&(def.Npoints_z), }},          
        {hmpdf_z_min,
            &(d->n->zmin),      dbl_type,     0, {&(def.z_min), }},
        {hmpdf_z_max,
            &(d->n->zmax),      dbl_type,     1, {NULL, &(def.z_max)}},
        {hmpdf_z_source,
            &(d->n->zsource),       dbl_type,     0, {&(def.z_source), }},
        {hmpdf_N_M,
            &(d->n->NM),        int_type,     0, {&(def.Npoints_M), }},
        {hmpdf_M_min,
            &(d->n->Mmin),      dbl_type,     0, {&(def.M_min), }},
        {hmpdf_M_max,
            &(d->n->Mmax),      dbl_type,     0, {&(def.M_max), }},
        {hmpdf_N_signal,
            &(d->n->Nsignal),   int_type,     0, {&(def.Npoints_signal), }},
        {hmpdf_signal_min,
            &(d->n->signalmin), dbl_type,     0, {&(def.signal_min), }},
        {hmpdf_signal_max,
            &(d->n->signalmax), dbl_type,     1, {&(def.max_kappa), &(def.max_tsz)}},
        {hmpdf_N_theta,
            &(d->p->Ntheta),        int_type,     0, {&(def.Npoints_theta), }},
        {hmpdf_rout_scale,
            &(d->p->rout_scale),    dbl_type,     0, {&(def.rout_scale), }},
        {hmpdf_rout_rdef,
            &(d->p->rout_def),      mdef_type,    0, {&(def.rout_rdef), }},
        {hmpdf_pixel_side,
            &(d->f->pixelside),     dbl_type,     0, {&(def.pixel_sidelength), }},
        {hmpdf_tophat_radius,
            &(d->f->tophat_radius), dbl_type,     0, {&(def.tophat_radius), }},
        {hmpdf_gaussian_fwhm,
            &(d->f->gaussian_sigma),dbl_type,     0, {&(def.gaussian_fwhm), }},
        {hmpdf_custom_ell_filter,
            &(d->f->custom_ell),    lf_type,      0, {&(def.custom_ell_filter), }},
        {hmpdf_custom_ell_filter_params,
            &(d->f->custom_ell_p),  vptr_type,    0, {&(def.custom_ell_filter_params), }},
        {hmpdf_custom_k_filter,
            &(d->f->custom_k),      kf_type,      0, {&(def.custom_k_filter), }},
        {hmpdf_custom_k_filter_params,
            &(d->f->custom_k_p),    vptr_type,    0, {&(def.custom_k_filter_params), }},
        {hmpdf_N_phi,
            &(d->n->Nphi),      int_type,     0, {&(def.Nphi), }},
        {hmpdf_phi_max,
            &(d->n->phimax),    dbl_type,     0, {&(def.phimax), }},
        {hmpdf_pixelexact_max,
            &(d->n->pixelexactmax), int_type, 0, {&(def.pixelexactmax), }},
        {hmpdf_phi_jitter,
            &(d->n->phijitter), dbl_type,     0, {&(def.phijitter), }},
        {hmpdf_phi_pwr,
            &(d->n->phipwr),    dbl_type,     0, {&(def.phipwr), }},
        {hmpdf_regularize_tp,
            &(d->tp->regularize),   int_type,     0, {&(def.regularize_tp), }},
        {hmpdf_monotonize,
            &(d->n->monotonize),    int_type,     0, {&(def.monotonize), }},
        {hmpdf_zintegr_type,
            &(d->n->zintegr_type),  integr_type,  0, {&(def.zintegr_type), }},
        {hmpdf_zintegr_alpha,
            &(d->n->zintegr_alpha), dbl_type,     0, {&(def.zintegr_alpha), }},
        {hmpdf_zintegr_beta,
            &(d->n->zintegr_beta),  dbl_type,     0, {&(def.zintegr_beta), }},
        {hmpdf_Mintegr_type,
            &(d->n->Mintegr_type),  integr_type,  0, {&(def.Mintegr_type), }},
        {hmpdf_Mintegr_alpha,
            &(d->n->Mintegr_alpha), dbl_type,     0, {&(def.Mintegr_alpha), }},
        {hmpdf_Mintegr_beta,
            &(d->n->Mintegr_beta),  dbl_type,     0, {&(def.Mintegr_beta), }},
        {hmpdf_Duffy08_conc_params,
            &(d->h->Duffy08_params), dptr_type,   0, {&(def.Duffy08_p), }},
        {hmpdf_Tinker10_hmf_params,
            &(d->h->Tinker10_params),dptr_type,   0, {&(def.Tinker10_p), }},
        {hmpdf_Battaglia12_tsz_params,
            &(d->p->Battaglia12_params), dbl_type,0, {&(def.Battaglia12_p), }},
    };//}}}

    for (int ii=0; ii<hmpdf_end_configs; ii++)
    {
        if (ii != p[ii].indx) { printf("ERROR : Indx not matching at %d.\n", ii); return; }
        p[ii].set = 0;
    }

    va_list valist;
    va_start(valist, stype);
    for (;;)
    {
        configs c = va_arg(valist, configs);
        if (c == hmpdf_end_configs)
        {
            break;
        }
        p[c].set = 1;
        switch (p[c].dt)
        {//{{{
            case (str_type)   : *((char **)(p[c].target)) = va_arg(valist, char *);
                                 break;
            case (int_type)    : *((int *)(p[c].target)) = va_arg(valist, int);
                                 break;
            case (dbl_type)    : *((double *)(p[c].target)) = va_arg(valist, double);
                                 break;
            case (dptr_type)   : *((double **)(p[c].target)) = va_arg(valist, double *);
                                 break;
            case (mdef_type)   : *((mdef *)(p[c].target)) = va_arg(valist, mdef);
                                 break;
            case (integr_type) : *((integr_mode *)(p[c].target)) = va_arg(valist, integr_mode);
                                 break;
            case (lf_type)     : *((ell_filter *)(p[c].target)) = va_arg(valist, ell_filter);
                                 break;
            case (kf_type)     : *((k_filter *)(p[c].target)) = va_arg(valist, k_filter);
                                 break;
            case (vptr_type)   : *((void **)(p[c].target)) = va_arg(valist, void *);
                                 break;
            default            : fprintf(stderr, "Not implemented : dtype.\n");
                                 break;
        }//}}}
    }
    va_end(valist);

    for (int ii=0; ii<hmpdf_end_configs; ii++)
    {
        if (p[ii].set)
        {
            continue;
        }
        if ((p[ii].depends_on_stype && p[ii].def[stype]==NULL)
            || (!p[ii].depends_on_stype && p[ii].def[0]==NULL))
        {
            fprintf(stderr, "Required setting %d for given signal_type not given.\n", ii);
        }
        p[ii].set = 1;
        switch (p[ii].dt)
        {//{{{
            case (str_type)      : *((char **)(p[ii].target)) =
                                        (p[ii].depends_on_stype) ? 
                                        *((char **)(p[ii].def[stype]))
                                        : *((char **)(p[ii].def[0]));
                                   break;
            case (int_type)      : *((int *)(p[ii].target)) =
                                        (p[ii].depends_on_stype) ? 
                                        *((int *)(p[ii].def[stype]))
                                        : *((int *)(p[ii].def[0]));
                                   break;
            case (dbl_type)      : *((double *)(p[ii].target)) =
                                        (p[ii].depends_on_stype) ? 
                                        *((double *)(p[ii].def[stype]))
                                        : *((double *)(p[ii].def[0]));
                                   break;
            case (dptr_type)     : *((double **)(p[ii].target)) =
                                        (p[ii].depends_on_stype) ? 
                                        *((double **)(p[ii].def[stype]))
                                        : *((double **)(p[ii].def[0]));
                                   break;
            case (mdef_type)     : *((mdef *)(p[ii].target)) =
                                        (p[ii].depends_on_stype) ? 
                                        *((mdef *)(p[ii].def[stype]))
                                        : *((mdef *)(p[ii].def[0]));
                                   break;
            case (integr_type)   : *((integr_mode *)(p[ii].target)) =
                                        (p[ii].depends_on_stype) ? 
                                        *((integr_mode *)(p[ii].def[stype]))
                                        : *((integr_mode *)(p[ii].def[0]));
                                   break;
            case (lf_type)       : *((ell_filter *)(p[ii].target)) =
                                        (p[ii].depends_on_stype) ? 
                                        *((ell_filter *)(p[ii].def[stype]))
                                        : *((ell_filter *)(p[ii].def[0]));
                                   break;
            case (kf_type)       : *((k_filter *)(p[ii].target)) =
                                        (p[ii].depends_on_stype) ? 
                                        *((k_filter *)(p[ii].def[stype]))
                                        : *((k_filter *)(p[ii].def[0]));
                                   break;
            case (vptr_type)     : *((void **)(p[ii].target)) =
                                        (p[ii].depends_on_stype) ? 
                                        *((void **)(p[ii].def[stype]))
                                        : *((void **)(p[ii].def[0]));
                                   break;
            default              : fprintf(stderr, "Not implemented : dtype.\n");
                                   break;
        }//}}}
    }

    for (int ii=0; ii<hmpdf_end_configs; ii++)
    {
        if (!(p[ii].set))
        {
            fprintf(stderr, "Required parameter %d not set.\n", ii);
        }
    }

    // sanity check
    #ifndef _OPENMP
    if (d->Ncores > 1)
    {
        printf("Warning : You requested N_cores = %d, "
               "but code is compiled without OpenMP.\n", d->Ncores);
    }
    #endif

    // do necessary conversions
    d->n->phimax *= M_PI/180.0/60.0;
    d->f->pixelside *= M_PI/180.0/60.0;
    d->f->tophat_radius *= M_PI/180.0/60.0;
    d->f->gaussian_sigma *= M_PI/180.0/60.0/sqrt(8.0*M_LN2); // convert FWHM (input) to sigma

    // compute things that we need for all output products
    init_numerics(d);
    init_class_interface(d);
    init_cosmology(d);
    init_power(d);
    init_halo_model(d);
    init_filters(d);
    init_profiles(d);
}//}}}
