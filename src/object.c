#include "utils.h"
#include "object.h"

#include "hmpdf.h"

int
null_data(hmpdf_obj *d)
{//{{{
    STARTFCT

    d->inited = 0;

    SAFEHMPDF(null_numerics(d));
    SAFEHMPDF(null_class_interface(d));
    SAFEHMPDF(null_cosmology(d));
    SAFEHMPDF(null_power(d));
    SAFEHMPDF(null_halo_model(d));
    SAFEHMPDF(null_filters(d));
    SAFEHMPDF(null_profiles(d));
    SAFEHMPDF(null_noise(d));
    SAFEHMPDF(null_onepoint(d));
    SAFEHMPDF(null_twopoint(d));
    SAFEHMPDF(null_powerspectrum(d));
    SAFEHMPDF(null_covariance(d));
    SAFEHMPDF(null_maps(d));

    ENDFCT
}//}}}

#define HMPDFNEW_ALLOC(var,expr)     \
    do {                             \
        var = expr;                  \
        if (UNLIKELY(!(var)))        \
        {                            \
            return NULL;             \
        }                            \
    } while(0)

hmpdf_obj *
hmpdf_new(void)
{//{{{
    hmpdf_obj *d;

    HMPDFNEW_ALLOC(d, malloc(sizeof(hmpdf_obj)));

    d->inited = 0;

    // if user calls hmpdf_delete directly after hmpdf_new,
    //     we would like this to have a defined value
    d->verbosity = 0;

    HMPDFNEW_ALLOC(d->n,   malloc(sizeof(numerics_t)));
    HMPDFNEW_ALLOC(d->cls, malloc(sizeof(class_interface_t)));
    HMPDFNEW_ALLOC(d->c,   malloc(sizeof(cosmology_t)));
    HMPDFNEW_ALLOC(d->pwr, malloc(sizeof(power_t)));
    HMPDFNEW_ALLOC(d->h,   malloc(sizeof(halo_model_t)));
    HMPDFNEW_ALLOC(d->f,   malloc(sizeof(filters_t)));
    HMPDFNEW_ALLOC(d->p,   malloc(sizeof(profiles_t)));
    HMPDFNEW_ALLOC(d->ns,  malloc(sizeof(noise_t)));
    HMPDFNEW_ALLOC(d->op,  malloc(sizeof(onepoint_t)));
    HMPDFNEW_ALLOC(d->tp,  malloc(sizeof(twopoint_t)));
    HMPDFNEW_ALLOC(d->ps,  malloc(sizeof(powerspectrum_t)));
    HMPDFNEW_ALLOC(d->cov, malloc(sizeof(covariance_t)));
    HMPDFNEW_ALLOC(d->m,   malloc(sizeof(maps_t)));

    int status = null_data(d);
    
    if (UNLIKELY(status || errno))
    {
        return NULL;
    }
    else
    {
        return d;
    }
}//}}}

#undef HMPDFNEW_ALLOC

int
reset_obj(hmpdf_obj *d)
{//{{{
    STARTFCT

    SAFEHMPDF(reset_numerics(d));
    SAFEHMPDF(reset_cosmology(d));
    SAFEHMPDF(reset_class_interface(d));
    SAFEHMPDF(reset_power(d));
    SAFEHMPDF(reset_halo_model(d));
    SAFEHMPDF(reset_filters(d));
    SAFEHMPDF(reset_noise(d));
    SAFEHMPDF(reset_onepoint(d));
    SAFEHMPDF(reset_twopoint(d));
    SAFEHMPDF(reset_powerspectrum(d));
    SAFEHMPDF(reset_covariance(d));
    SAFEHMPDF(reset_profiles(d));
    SAFEHMPDF(reset_maps(d));

    SAFEHMPDF(null_data(d));

    ENDFCT
}//}}}

int
hmpdf_delete(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(1, "hmpdf_delete\n");

    SAFEHMPDF(reset_obj(d));

    free(d->cls);
    free(d->c);
    free(d->pwr);
    free(d->f);
    free(d->p);
    free(d->h);
    free(d->n);
    free(d->ns);
    free(d->op);
    free(d->tp);
    free(d->ps);
    free(d->cov);
    free(d->m);

    free(d);
     
    ENDFCT
}//}}}

