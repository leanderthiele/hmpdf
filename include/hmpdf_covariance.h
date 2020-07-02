/*! @file hmpdf_covariance.h */
#ifndef HMPDF_COVARIANCE_H
#define HMPDF_COVARIANCE_H

/*! Returns the covariance matrix of the one-point PDF.
 *
 *  \param[in,out] d    hmpdf_init() must have been called on d
 *  \param[in] Nbins    number of bins the covariance matrix will be binned into
 *  \param[in] binedges monotonically increasing array of length Nbins+1
 *  \param[out] cov     the binned covariance matrix will be written into the first Nbins*Nbins elements of 
 *                      this output array
 *  \param[in] noisy    if set to non-zero, the covariance matrix will include Gaussian noise
 *                      of standard deviation #hmpdf_noise
 *  \return error code
 *
 *  \remark If the covariance matrix has already been computed and since then no hmpdf_init()
 *          has been called on d, the pre-computed result is used and only the binning is performed.
 *  \remark If #hmpdf_verbosity is set to a positive value, status updates with estimated remaining
 *          time will be given during execution.
 */
int hmpdf_get_cov(hmpdf_obj *d,
                  int Nbins,
                  double binedges[Nbins+1],
                  double cov[Nbins*Nbins],
                  int noisy);

/*! Returns diagnostic outputs for the covariance matrix computation.
 *  The main use of this function is to identify numerical instability at small
 *  pixel separations.
 *
 *  \param[in,out] d        hmpdf_init() must have been called on d
 *  \param[out] Nphi        the number of pixel separations will be written into the return value
 *  \param[out] phi         pointer will be set to an array of length Nphi,
 *                          containing the pixel separations used internally (in radians)
 *  \param[out] phiweights  pointer will be set to an array of length Nphi,
 *                          containing the weights assigned to each pixel separation
 *                          in the summation.
 *                          This can occasionally be used to tune the #hmpdf_pixelexact_max option.
 *  \param[out] corr_diagn  pointer will be set to an array of length Nphi,
 *                          containing the correlation function at the pixel separation sample points.
 *                          Noisy behaviour at small phi is a sign of numerical instability.
 *  \return error code
 *
 *  \remark the values in the phi-array are not ordered
 *  \remark you can pass NULL for any of the arguments (apart from the first of course),
 *          if you are not interested in that specific output
 *  \remark if the covariance matrix has not already been computed
 *          [hmpdf_get_cov() has not been called before],
 *          this function will do it
 *  \remark while the code does perform the memory allocation for phi, phiweights, and corr_diagn
 *          (so that the user does not have to figure out Nphi beforehand),
 *          the user is responsible for freeing these arrays, i.e. to call
 *          \code
 *          free(*phi); free(*phiweights); free(*corr_diagn);
 *          \endcode
 *          after use.
 */
int hmpdf_get_cov_diagnostics(hmpdf_obj *d,
                              int *Nphi,
                              double **phi,
                              double **phiweights,
                              double **corr_diagn);

#endif
