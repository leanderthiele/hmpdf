/*! @file hmpdf_twopoint.h */
#ifndef HMPDF_TWOPOINT_H
#define HMPDF_TWOPOINT_H

/*! Returns the two-point PDF.
 *
 *  \param d        hmpdf_init() must have been called on d
 *  \param phi      angular separation for which the two-point PDF will be computed (in radians).
 *                  Must be strictly positive, and less than #hmpdf_phi_max
 *                  (note different units!)
 *  \param Nbins    number of bins the two-point PDF will be binned into
 *  \param binedges monotonically increasing array of length Nbins+1
 *  \param out      the binned two-point PDF will be written into the first Nbins*Nbins elements of 
 *                  this output array
 *  \param noisy    if set to non-zero, the two-point PDF will be convolved with a Gaussian
 *                  of standard deviation #hmpdf_noise
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
                 double out[Nbins*Nbins],
                 int noisy);


#endif
