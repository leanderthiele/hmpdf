/*! @file hmpdf_onepoint.h */
#ifndef HMPDF_ONEPOINT_H
#define HMPDF_ONEPOINT_H

#include "hmpdf_object.h"

/*! Returns the one-point PDF.
 *
 *  \param[in,out] d    hmpdf_init() must have been called on d
 *  \param[in] Nbins    number of bins the one-point PDF will be binned into
 *  \param[in] binedges monotonically increasing array of length Nbins+1
 *  \param[out] op      the binned one-point PDF will be written into the first Nbins elements of 
 *                      this output array
 *  \param[in] incl_2h  if set to non-zero, the two-halo term (i.e. halo clustering)
 *                      will be included in the output PDF
 *  \param[in] noisy    if set to non-zero, the one-point PDF will be convolved with a Gaussian
 *                      of standard deviation #hmpdf_noise
 *  \return error code
 *
 *  \remark if the one-point PDF has already been computed and since then no hmpdf_init()
 *          has been called on d, the pre-computed result is used and only the binning is performed.
 */
int hmpdf_get_op(hmpdf_obj *d,
                 int Nbins,
                 double binedges[Nbins+1],
                 double op[Nbins],
                 int incl_2h,
                 int noisy);

#endif
