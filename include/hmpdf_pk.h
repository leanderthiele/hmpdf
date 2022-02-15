/*! @file hmpdf_pk.h */
#ifndef HMPDF_PK_H
#define HMPDF_PK_H

#include "hmpdf_object.h"

/*! Returns the 3D power spectrum.
 *
 *  \param[in, out] d   hmpdf_init() must have been called on d
 *  \param[in] Npoints  number of k stencils
 *  \param[in] k        k stencils in 1/Mpc
 *  \param[out] Pk_1h   output array for 1-halo term
 *  \param[out] Pk_2h   output array for 2-halo term
 *  \return error code
 */
int hmpdf_get_pk(hmpdf_obj *d,
                 double z,
                 int Npoints, 
                 double k[Npoints],
                 double Pk_1h[Npoints],
                 double Pk_2h[Npoints]);

#endif
