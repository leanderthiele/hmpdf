/*! @file hmpdf_maps.h */
#ifndef HMPDF_MAPS_H
#define HMPDF_MAPS_H

#include "hmpdf_object.h"

/*! Returns the histogram of a simplified simulation (map).
 *
 *  \param[in,out] d    hmpdf_init() must have been called on d
 *  \param[in] Nbins    number of bins the histogram will be binned into
 *  \param[in] binedges monotonically increasing array of length Nbins+1
 *  \param[out] op      the binned histogram, normalized
 *  \param[in] new_map  if set to non-zero, the simplified simulation will
 *                      be rerun even if a map has already been generated
 *  \return error code
 */
int hmpdf_get_map_op(hmpdf_obj *d,
                     int Nbins,
                     double binedges[Nbins+1],
                     double op[Nbins],
                     int new_map);

/*! Returns the power spectrum of a simplified simulation (map).
 *
 *  \param[in,out] d    hmpdf_init() must have been called on d
 *  \param[in] Nbins    number of bins the power spectrum will be binned into
 *  \param[in] binedges monotonically increasing array of length Nbins+1
 *  \param[out] ps      the binned, direction averaged power spectrum
 *  \param[in] new_map  if set to non-zero, the simplified simulation will
 *                      be rerun even if a map has already been generated
 *  \return error code
 */
int hmpdf_get_map_ps(hmpdf_obj *d,
                     int Nbins,
                     double binedges[Nbins+1],
                     double ps[Nbins],
                     int new_map);

/*! Returns a simplified simulation (map).
 *
 *  \param[in,out] d    hmpdf_init() must have been called on d
 *  \param[out] map     the map (flattened array of dimensions Nside x Nside)
 *  \param[out] Nside   sidelength of the map
 *  \param[in] new_map  if set to non-zero, the simplified simulation will
 *                      be rerun even if a map has already been generated
 *  \return error code
 *
 *  \remark while the code does perform the memory allocation for the map
 *          (so the user does not have to figure out Nside beforehand),
 *          the user is responsible for freeing the allocated array, i.e. to call
 *          \code
 *          free(*map)
 *          \endcode
 *          after use.
 */
int hmpdf_get_map(hmpdf_obj *d,
                  double **map,
                  long *Nside,
                  int new_map);

#endif
