#include <stdio.h>
#include <stdarg.h>
#include <math.h>

#include "utils.h"
#include "configs.h"
#include "object.h"
#include "class_interface.h"
#include "cosmology.h"
#include "numerics.h"
#include "power.h"
#include "halo_model.h"
#include "filter.h"
#include "profiles.h"
#include "noise.h"
#include "init.h"

#include "hmpdf.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

typedef enum
{//{{{
    int_type, // int
    long_type, // long
    dbl_type, // double
    mdef_type, // hmpdf_mdef_e
    integr_type, // hmpdf_integr_mode_e
    end_comparable_dtypes,
    str_type, // char *
    dptr_type, // double *
    vptr_type, // void *
    lf_type, // hmpdf_ell_filter_f
    kf_type, // hmpdf_k_filter_f
}//}}}
dtype;

// perform action depending on data type
//DT_DEP_ACTION{{{
#define DT_DEP_ACTION(dt, expr)                                \
    do {                                                       \
    switch (dt)                                                \
    {                                                          \
        case (str_type) : expr(char *); break;                 \
        case (int_type) : expr(int); break;                    \
        case (long_type) : expr(long); break;                \
        case (dbl_type) : expr(double); break;                 \
        case (dptr_type) : expr(double *); break;              \
        case (mdef_type) : expr(hmpdf_mdef_e); break;          \
        case (integr_type) : expr(hmpdf_integr_mode_e); break; \
        case (vptr_type) : expr(void *); break;                \
        case (lf_type) : expr(hmpdf_ell_filter_f); break;      \
        case (kf_type) : expr(hmpdf_k_filter_f); break;        \
        default : HMPDFERR("Unknown dtype.");                  \
                  break;                                       \
    }                                                          \
    } while (0)
//}}}

typedef struct
{//{{{
    char *name;
    void *target;
    dtype dt;
    void *def;
    void *lo;
    void *hi;
    int set;
}//}}}
param;

// convenience function to save typing
static void
init_p(param *p, char *name, void *target, dtype dt, void *def, void *lo, void *hi)
{//{{{
    p->name = name;
    p->target = target;
    p->dt = dt;
    p->def = def;
    p->lo = lo;
    p->hi = hi;
    p->set = 0;
}//}}}

// convenience macros to save even more typing
//INIT_P*{{{
#define INIT_P(indx, targ, dt, df)                         \
    do {                                                   \
    init_p(p+indx, #indx, &(targ), dt, &(df), NULL, NULL); \
    ++ctr;                                                 \
    } while (0)

#define INIT_P_B(indx, targ, dt, df)                                  \
    do {                                                              \
    init_p(p+indx, #indx, &(targ), dt, &(df[0]), &(df[1]), &(df[2])); \
    ++ctr;                                                            \
    } while (0)

#define INIT_P2(indx, targ, dt, dfk, dft)                \
    do {                                                 \
    init_p(p+indx, #indx, &(targ), dt,                   \
           (d->p->stype==hmpdf_kappa) ? &(dfk) : &(dft), \
           NULL, NULL);                                  \
    ++ctr;                                               \
    } while (0)

#define INIT_P2_BK(indx, targ, dt, dfk, dft)                \
    do {                                                    \
    init_p(p+indx, #indx, &(targ), dt,                      \
           (d->p->stype==hmpdf_kappa) ? &(dfk[0]) : &(dft), \
           (d->p->stype==hmpdf_kappa) ? &(dfk[1]) : NULL,   \
           (d->p->stype==hmpdf_kappa) ? &(dfk[2]) : NULL);  \
    ++ctr;                                                  \
    } while (0)

#define INIT_P2_BT(indx, targ, dt, dfk, dft)                 \
    do {                                                     \
    init_p(p+indx, #indx, &(targ), dt,                       \
           (d->p->stype==hmpdf_kappa) ? &(dfk) : &(dft[0]),  \
           (d->p->stype==hmpdf_kappa) ? NULL   : &(dft[1]),  \
           (d->p->stype==hmpdf_kappa) ? NULL   : &(dft[2])); \
    ++ctr;                                                   \
    } while (0)

#define INIT_P2_BKT(indx, targ, dt, dfk, dft)                   \
    do {                                                        \
    init_p(p+indx, #indx, &(targ), dt,                          \
           (d->p->stype==hmpdf_kappa) ? &(dfk[0]) : &(dft[0]),  \
           (d->p->stype==hmpdf_kappa) ? &(dfk[1]) : &(dft[1]),  \
           (d->p->stype==hmpdf_kappa) ? &(dfk[2]) : &(dft[2])); \
    ++ctr;                                                      \
    } while (0)
//}}}

static int 
init_params(hmpdf_obj *d, param *p)
{//{{{
    STARTFCT

    int ctr = 0;
    INIT_P_B(hmpdf_N_threads,
             d->Ncores, int_type, def.Ncores);
    INIT_P(hmpdf_verbosity,
           d->verbosity, int_type, def.verbosity);
    INIT_P(hmpdf_warn_is_err,
           d->warn_is_err, int_type, def.warn_is_err);
    INIT_P(hmpdf_class_pre,
           d->cls->class_pre, str_type, def.class_pre);
    INIT_P_B(hmpdf_N_z,
             d->n->Nz, int_type, def.Npoints_z);
    INIT_P_B(hmpdf_z_min,
             d->n->zmin, dbl_type, def.z_min);
    INIT_P2_BT(hmpdf_z_max,
               d->n->zmax, dbl_type, d->n->zsource, def.z_max);
    INIT_P_B(hmpdf_N_M,
             d->n->NM, int_type, def.Npoints_M);
    INIT_P_B(hmpdf_M_min,
             d->n->Mmin, dbl_type, def.M_min);
    INIT_P_B(hmpdf_M_max,
             d->n->Mmax, dbl_type, def.M_max);
    INIT_P_B(hmpdf_N_signal,
             d->n->Nsignal, long_type, def.Npoints_signal);
    INIT_P2_BKT(hmpdf_signal_min,
               d->n->signalmin, dbl_type, def.min_kappa, def.min_tsz);
    INIT_P2_BKT(hmpdf_signal_max,
                d->n->signalmax, dbl_type, def.max_kappa, def.max_tsz);
    INIT_P_B(hmpdf_N_theta,
             d->p->Ntheta, int_type, def.Npoints_theta);
    INIT_P_B(hmpdf_rout_scale,
             d->p->rout_scale, dbl_type, def.rout_scale);
    INIT_P_B(hmpdf_rout_rdef,
             d->p->rout_def, mdef_type, def.rout_rdef);
    INIT_P_B(hmpdf_pixel_side,
             d->f->pixelside, dbl_type, def.pixel_sidelength);
    INIT_P_B(hmpdf_tophat_radius,
             d->f->tophat_radius, dbl_type, def.tophat_radius);
    INIT_P_B(hmpdf_gaussian_fwhm,
             d->f->gaussian_sigma, dbl_type, def.gaussian_fwhm);
    INIT_P(hmpdf_custom_ell_filter,
           d->f->custom_ell, lf_type, def.custom_ell_filter);
    INIT_P(hmpdf_custom_ell_filter_params,
           d->f->custom_ell_p, vptr_type, def.custom_ell_filter_params);
    INIT_P(hmpdf_custom_k_filter,
           d->f->custom_k, kf_type, def.custom_k_filter);
    INIT_P(hmpdf_custom_k_filter_params,
           d->f->custom_k_p, vptr_type, def.custom_k_filter_params);
    INIT_P_B(hmpdf_N_phi,
             d->n->Nphi, int_type, def.Nphi);
    INIT_P_B(hmpdf_phi_max,
             d->n->phimax, dbl_type, def.phimax);
    INIT_P_B(hmpdf_pixelexact_max,
             d->n->pixelexactmax, int_type, def.pixelexactmax);
    INIT_P_B(hmpdf_phi_jitter,
             d->n->phijitter, dbl_type, def.phijitter);
    INIT_P(hmpdf_phi_pwr,
           d->n->phipwr, dbl_type, def.phipwr);
    INIT_P_B(hmpdf_zintegr_type,
             d->n->zintegr_type, integr_type, def.zintegr_type);
    INIT_P(hmpdf_zintegr_alpha,
           d->n->zintegr_alpha, dbl_type, def.zintegr_alpha);
    INIT_P(hmpdf_zintegr_beta,
           d->n->zintegr_beta, dbl_type, def.zintegr_beta);
    INIT_P_B(hmpdf_Mintegr_type,
             d->n->Mintegr_type, integr_type, def.Mintegr_type);
    INIT_P(hmpdf_Mintegr_alpha,
           d->n->Mintegr_alpha, dbl_type, def.Mintegr_alpha);
    INIT_P(hmpdf_Mintegr_beta,
           d->n->Mintegr_beta, dbl_type, def.Mintegr_beta);
    INIT_P(hmpdf_Duffy08_conc_params,
           d->h->Duffy08_params, dptr_type, def.Duffy08_p);
    INIT_P(hmpdf_Tinker10_hmf_params,
           d->h->Tinker10_params, dptr_type, def.Tinker10_p);
    INIT_P(hmpdf_Battaglia12_tsz_params,
           d->p->Battaglia12_params, dptr_type, def.Battaglia12_p);
    INIT_P_B(hmpdf_noise,
             d->ns->noise, dbl_type, def.noise);

    HMPDFCHECK(ctr != hmpdf_end_configs, "Not all params filled, ctr = %d.", ctr);

    ENDFCT
}//}}}

#undef INIT_P
#undef INIT_P_B
#undef INIT_P2
#undef INIT_P2_BK
#undef INIT_P2_BT
#undef INIT_P2_BKT

// assign parameter to value passed by user,
//     check for validity if bounds are present
//ASSIGN_SET{{{
#define ASSIGN_SET(dt)                         \
    do {                                       \
    *((dt*)(p->target)) = va_arg(*valist, dt); \
    } while (0)
//}}}

static int 
assign_set(param *p, va_list *valist)
{//{{{
    STARTFCT

    DT_DEP_ACTION(p->dt, ASSIGN_SET);

    ENDFCT
}//}}}

#undef ASSIGN_SET

// check if user setting is reasonable
//CHECK_VALIDITY{{{
#define CHECK_VALIDITY(dt)                             \
    do {                                               \
    if (comparable)                                    \
    {                                                  \
        if (p->lo != NULL)                             \
        {                                              \
            if (*((dt*)(p->target)) < *((dt*)(p->lo))) \
            {                                          \
                *invalid_param = 1;                    \
            }                                          \
        }                                              \
        if (p->hi != NULL)                             \
        {                                              \
            if (*((dt*)(p->target)) > *((dt*)(p->hi))) \
            {                                          \
                *invalid_param = 1;                    \
            }                                          \
        }                                              \
    }                                                  \
    } while (0)
//}}}

static int
check_validity(param *p, int *invalid_param)
{//{{{
    STARTFCT

    int comparable = (p->dt < end_comparable_dtypes) ? 1 : 0;
    DT_DEP_ACTION(p->dt, CHECK_VALIDITY);

    ENDFCT
}//}}}

#undef CHECK_VALIDITY

// assign parameter to default value
//ASSIGN_DEF{{{
#define ASSIGN_DEF(dt)                      \
    do {                                    \
    *((dt*)(p->target)) = *((dt*)(p->def)); \
    } while (0)
//}}}

static int
assign_def(param *p)
{//{{{
    STARTFCT

    DT_DEP_ACTION(p->dt, ASSIGN_DEF);

    ENDFCT
}//}}}

#undef ASSIGN_DEF

// print depending on data type
//FMT{{{
#define FMT(dt)                \
    (dt==str_type) ? "%s"      \
    : (dt==int_type) ? "%d"    \
    : (dt==long_type) ? "%ld"  \
    : (dt==dbl_type) ? "%g"    \
    : (dt==dptr_type) ? "%p"   \
    : (dt==mdef_type) ? "%d"   \
    : (dt==integr_type) ? "%d" \
    : (dt==vptr_type) ? "%p"   \
    : (dt==lf_type) ? "%p"     \
    : (dt==kf_type) ? "%p"     \
    : "%d"

#define PRINTVAL(dt1)                                \
    do {                                             \
    *printed += sprintf(f, FMT(dt), *((dt1*)(val))); \
    } while (0)
//}}}

static int
printval(char *f, void *val, dtype dt, int *printed)
{//{{{
    STARTFCT

    DT_DEP_ACTION(dt, PRINTVAL);

    ENDFCT
}//}}}

#undef FMT
#undef PRINTVAL
#undef DT_DEP_ACTION

static int
invalid_param_warn(hmpdf_obj *d, param *p)
{//{{{
    STARTFCT

    char msg[512];
    int printed = 0;
    printed += sprintf(msg+printed,
                       "option passed for %s "
                       "is out of recommended bounds:\n"
                       "\tyou passed ", p->name);
    SAFEHMPDF(printval(msg+printed, p->target, p->dt, &printed));
    printed += sprintf(msg+printed,
                       " which is outside the bounds ( ");
    SAFEHMPDF(printval(msg+printed, p->lo, p->dt, &printed));
    printed += sprintf(msg+printed, " , ");
    SAFEHMPDF(printval(msg+printed, p->hi, p->dt, &printed));
    printed += sprintf(msg+printed, " )\n");
    HMPDFWARN("%s", msg);

    ENDFCT
}//}}}

static int
unit_conversions(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->n->phimax         *= RADPERARCMIN;
    d->f->pixelside      *= RADPERARCMIN;
    d->f->tophat_radius  *= RADPERARCMIN;
    d->f->gaussian_sigma *= RADPERARCMIN/sqrt(8.0*M_LN2); // convert FWHM (input) to sigma

    ENDFCT
}//}}}

static int
sanity_checks(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFCHECK((d->p->stype != hmpdf_tsz) && (d->p->stype != hmpdf_kappa),
               "Invalid signal type %d.", d->p->stype);
    HMPDFCHECK((d->p->stype==hmpdf_kappa)
                && ((d->n->zsource < (d->n->zmax-0.001)) || (d->n->zsource > 1200.0)),
                "Invalid source redshift %g.", d->n->zsource);
    #ifndef _OPENMP
    HMPDFCHECK(d->Ncores>1, "You specified hmpdf_N_threads = %d, "
                            "but code is compiled without OpenMP.", d->Ncores);
    #endif
    HMPDFCHECK(d->n->signalmin >= d->n->signalmax,
               "hmpdf_signal_min must be less than hmpdf_signal_max.");
    HMPDFCHECK(d->n->zmin >= d->n->zmax,
               "hmpdf_z_min must be less than hmpdf_z_max.");
    HMPDFCHECK(d->n->Mmin >= d->n->Mmax,
               "hmpdf_M_min must be less than hmpdf_M_max.");

    ENDFCT
}//}}}

static int
compute_necessary_for_all(hmpdf_obj *d)
{//{{{
    STARTFCT

    SAFEHMPDF(init_numerics(d));
    SAFEHMPDF(init_class_interface(d));
    SAFEHMPDF(init_cosmology(d));
    SAFEHMPDF(init_power(d));
    SAFEHMPDF(init_halo_model(d));
    SAFEHMPDF(init_filters(d));
    SAFEHMPDF(init_profiles(d));
    SAFEHMPDF(init_noise(d));

    ENDFCT
}//}}}

int
hmpdf_init_fct(hmpdf_obj *d, char *class_ini, hmpdf_signaltype_e stype, ...)
{//{{{
    STARTFCT

    gsl_set_error_handler(&new_gsl_error_handler);

    d->inited = 0;

    d->cls->class_ini = class_ini;
    d->p->stype = stype;

    // this one is special because we need to have it set
    //     definitely when we want to check the validity
    //     of user inputs
    d->warn_is_err = def.warn_is_err;

    va_list valist;
    va_start(valist, stype);

    if (d->p->stype == hmpdf_kappa)
    {
        // zmax will be set to this value,
        // unless chosen differently.
        // Thus, do not move this past the init_params() call!
        d->n->zsource = va_arg(valist, double);
    }

    param *p;
    SAFEALLOC(p, malloc((int)(hmpdf_end_configs) * sizeof(param)));
    SAFEHMPDF(init_params(d, p));

    int read = 0;
    for (;;)
    {
        hmpdf_configs_e c = va_arg(valist, hmpdf_configs_e);
        if (c == hmpdf_end_configs)
        {
            break;
        }
        SAFEHMPDF(assign_set(p+(int)c, &valist));
        p[c].set = 1;

        ++read;
        HMPDFCHECK(read > hmpdf_end_configs,
                   "You likely forgot to end the variable argument list "
                   "with the mandatory argument hmpdf_end_configs, "
                   "or you passed options twice.");
    }

    va_end(valist);

    for (int ii=0; ii<hmpdf_end_configs; ii++)
    {
        if (p[ii].set)
        {
            int invalid_param = 0;
            SAFEHMPDF(check_validity(p+ii, &invalid_param));
            if (invalid_param)
            {
                SAFEHMPDF(invalid_param_warn(d, p+ii));
            }
        }
        else
        {
            SAFEHMPDF(assign_def(p+ii));
        }
    }

    free(p);

    // do necessary conversions
    SAFEHMPDF(unit_conversions(d));

    // perform basic sanity checks
    SAFEHMPDF(sanity_checks(d));

    // this frees all the computed quantities,
    // since we assume that each call of init changes some
    // parameter (cosmological or numerical)
    SAFEHMPDF(reset_obj(d));

    // compute things that we need for all output products
    SAFEHMPDF(compute_necessary_for_all(d));

    d->inited = 1;

    ENDFCT
}//}}}

