#ifndef OBJECT_H
#define OBJECT_H

#include "utils.h"
#include "numerics.h"
#include "class_interface.h"
#include "cosmology.h"
#include "power.h"
#include "halo_model.h"
#include "profiles.h"
#include "bcm.h"
#include "filter.h"
#include "noise.h"
#include "onepoint.h"
#include "twopoint.h"
#include "powerspectrum.h"
#include "covariance.h"
#include "maps.h"

#include "hmpdf.h"

struct hmpdf_obj_s
{//{{{
    int Ncores;
    int verbosity;
    int warn_is_err;
    int inited;

    numerics_t *n;
    class_interface_t *cls;
    cosmology_t *c;
    power_t *pwr;
    halo_model_t *h;
    filters_t *f;
    profiles_t *p;
    bcm_t *bcm;
    noise_t *ns;
    onepoint_t *op;
    twopoint_t *tp;
    powerspectrum_t *ps;
    covariance_t *cov;
    maps_t *m;
};//}}}

hmpdf_obj *hmpdf_new(void);
int reset_obj(hmpdf_obj *d);
int hmpdf_delete(hmpdf_obj *d);

#endif
