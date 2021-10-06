/*! @file hmpdf_mean.h */
#ifndef HMPDF_MEAN_H
#define HMPDF_MEAN_H

#include "hmpdf_object.h"

/*! Returns the sky-averaged mean signal.
 *
 * \param[in,out] d     hmpdf_init() must have been called on d
 * \param[out]    mean  return value
 *
 * \return error code
 *
 * \attention Does not take filters or noise into account.
 */
int hmpdf_get_mean(hmpdf_obj *d,
                   double *mean);


/*! Returns the sky-averaged Compton-y weighted temperature in keV.
 *
 * \param[in,out] d     hmpdf_init() must have been called on d
 * \param[out]    mean  return value
 *
 * \return error code
 *
 * \attention Does not take filters or noise into account.
 * \attention Only makes sense for signal type #hmpdf_tsz.
 */
int hmpdf_get_mean_T(hmpdf_obj *d,
                     double *mean);

#endif
