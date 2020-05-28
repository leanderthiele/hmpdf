#include "data.h"

#include "hmpdf.h"

void null_data(all_data *d)
{//{{{
    fprintf(stdout, "\tnull_data\n");
    fflush(stdout);

    null_numerics(d);
    null_class_interface(d);
    null_cosmology(d);
    null_power(d);
    null_halo_model(d);
    null_filters(d);
    null_profiles(d);
    null_noise(d);
    null_onepoint(d);
    null_twopoint(d);
    null_powerspectrum(d);
    null_covariance(d);
}//}}}

all_data *new_data(void)
{//{{{
    fprintf(stdout, "In data.h -> new_data.\n");
    fflush(stdout);
    all_data *d = malloc(sizeof(all_data));

    d->n = malloc(sizeof(numerics_t));
    d->cls = malloc(sizeof(class_interface_t));
    d->c = malloc(sizeof(cosmology_t));
    d->pwr = malloc(sizeof(power_t));
    d->h = malloc(sizeof(halo_model_t));
    d->f = malloc(sizeof(filters_t));
    d->p = malloc(sizeof(profiles_t));
    d->ns = malloc(sizeof(noise_t));
    d->op = malloc(sizeof(onepoint_t));
    d->tp = malloc(sizeof(twopoint_t));
    d->ps = malloc(sizeof(powerspectrum_t));
    d->cov = malloc(sizeof(covariance_t));
    null_data(d);
    return d;
}//}}}

void reset_data(all_data *d)
{//{{{
    fprintf(stdout, "In data.h -> reset_data.\n");
    fflush(stdout);

    reset_numerics(d);
    reset_cosmology(d);
    reset_class_interface(d);
    reset_power(d);
    reset_halo_model(d);
    reset_filters(d);
    reset_profiles(d);
    reset_noise(d);
    reset_onepoint(d);
    reset_twopoint(d);
    reset_powerspectrum(d);
    reset_covariance(d);

    null_data(d);
}//}}}

void delete_data(all_data *d)
{//{{{
    fprintf(stdout, "In data.h -> delete_data.\n");
    fflush(stdout);
    reset_data(d);

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
}//}}}

