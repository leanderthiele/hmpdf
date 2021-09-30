/*! @file hmpdf_twopoint.h */
#ifndef HMPDF_TWOPOINT_H
#define HMPDF_TWOPOINT_H

#include "hmpdf_object.h"

/*! Returns the two-point PDF.
 *
 *  \param[in,out] d    hmpdf_init() must have been called on d
 *  \param[in] phi      angular separation for which the two-point PDF will be computed (in arcmin).
 *                      Must be strictly positive, and less than #hmpdf_phi_max.
 *  \param[in] Nbins    number of bins the two-point PDF will be binned into
 *  \param[in] binedges monotonically increasing array of length Nbins+1
 *  \param[out] tp      the binned two-point PDF will be written into the first Nbins*Nbins elements of 
 *                      this output array
 *  \param[in] noisy    if set to non-zero, the two-point PDF will be convolved with a Gaussian
 *                      of covariance matrix determined from #hmpdf_noise_pwr
 *  \return error code
 *  
 *  \remark If the two-point PDF has already been computed with the same value of phi
 *          and since then no hmpdf_init() has been called on d
 *          or hmpdf_get_tp() has been called with a different value of phi,
 *          the pre-computed result is used and only the binning is performed.
 */
int hmpdf_get_tp(hmpdf_obj *d,
                 double phi,
                 int Nbins,
                 double binedges[Nbins+1],
                 double tp[Nbins*Nbins],
                 int noisy);


#endif
