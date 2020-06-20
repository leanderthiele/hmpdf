#ifndef DATA_H
#define DATA_H

#include "utils.h"
#include "numerics.h"
#include "class_interface.h"
#include "cosmology.h"
#include "power.h"
#include "halo_model.h"
#include "profiles.h"
#include "filter.h"
#include "noise.h"
#include "onepoint.h"
#include "twopoint.h"
#include "powerspectrum.h"
#include "covariance.h"

#include "hmpdf.h"

struct hmpdf_obj_s
{//{{{
    int Ncores;

    numerics_t *n;
    class_interface_t *cls;
    cosmology_t *c;
    power_t *pwr;
    halo_model_t *h;
    filters_t *f;
    profiles_t *p;
    noise_t *ns;
    onepoint_t *op;
    twopoint_t *tp;
    powerspectrum_t *ps;
    covariance_t *cov;
};//}}}

hmpdf_obj *hmpdf_new(void);
void reset_data(hmpdf_obj *d);
void hmpdf_delete(hmpdf_obj *d);

#endif
