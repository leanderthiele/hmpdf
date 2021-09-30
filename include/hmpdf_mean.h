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

#endif
