#include "utils.h"
#include "object.h"

#include "hmpdf.h"

int
null_data(hmpdf_obj *d)
{//{{{
    STARTFCT

    SAFEHMPDF(null_numerics(d))
    SAFEHMPDF(null_class_interface(d))
    SAFEHMPDF(null_cosmology(d))
    SAFEHMPDF(null_power(d))
    SAFEHMPDF(null_halo_model(d))
    SAFEHMPDF(null_filters(d))
    SAFEHMPDF(null_profiles(d))
    SAFEHMPDF(null_noise(d))
    SAFEHMPDF(null_onepoint(d))
    SAFEHMPDF(null_twopoint(d))
    SAFEHMPDF(null_powerspectrum(d))
    SAFEHMPDF(null_covariance(d))

    ENDFCT
}//}}}

hmpdf_obj *hmpdf_new(void)
{//{{{
    STARTFCT

    SAFEALLOC_NORETURN(hmpdf_obj *, d, malloc(sizeof(hmpdf_obj)))

    SAFEALLOC_NORETURN(, d->n, malloc(sizeof(numerics_t)))
    SAFEALLOC_NORETURN(, d->cls, malloc(sizeof(class_interface_t)))
    SAFEALLOC_NORETURN(, d->c, malloc(sizeof(cosmology_t)))
    SAFEALLOC_NORETURN(, d->pwr, malloc(sizeof(power_t)))
    SAFEALLOC_NORETURN(, d->h, malloc(sizeof(halo_model_t)))
    SAFEALLOC_NORETURN(, d->f, malloc(sizeof(filters_t)))
    SAFEALLOC_NORETURN(, d->p, malloc(sizeof(profiles_t)))
    SAFEALLOC_NORETURN(, d->ns, malloc(sizeof(noise_t)))
    SAFEALLOC_NORETURN(, d->op, malloc(sizeof(onepoint_t)))
    SAFEALLOC_NORETURN(, d->tp, malloc(sizeof(twopoint_t)))
    SAFEALLOC_NORETURN(, d->ps, malloc(sizeof(powerspectrum_t)))
    SAFEALLOC_NORETURN(, d->cov, malloc(sizeof(covariance_t)))
    SAFEHMPDF_NORETURN(null_data(d))
    
    if ((hmpdf_status) || (errno))
    {
        return NULL;
    }
    else
    {
        return d;
    }
}//}}}

int
reset_obj(hmpdf_obj *d)
{//{{{
    STARTFCT

    SAFEHMPDF(reset_numerics(d))
    SAFEHMPDF(reset_cosmology(d))
    SAFEHMPDF(reset_class_interface(d))
    SAFEHMPDF(reset_power(d))
    SAFEHMPDF(reset_halo_model(d))
    SAFEHMPDF(reset_filters(d))
    SAFEHMPDF(reset_noise(d))
    SAFEHMPDF(reset_onepoint(d))
    SAFEHMPDF(reset_twopoint(d))
    SAFEHMPDF(reset_powerspectrum(d))
    SAFEHMPDF(reset_covariance(d))
    SAFEHMPDF(reset_profiles(d))

    SAFEHMPDF(null_data(d))

    ENDFCT
}//}}}

int
hmpdf_delete(hmpdf_obj *d)
{//{{{
    STARTFCT

    HMPDFPRINT(1, "hmpdf_delete\n")

    SAFEHMPDF(reset_obj(d))

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

    free(d);
     
    ENDFCT
}//}}}

