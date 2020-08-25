/*! @file hmpdf_powerspectrum.h */
#ifndef HMPDF_POWERSPECTRUM_H
#define HMPDF_POWERSPECTRUM_H

#include "hmpdf_object.h"

/*! Power spectrum/correlation function output modes. */
typedef enum
{
    hmpdf_onehalo, /*!< output only the 1-halo term */
    hmpdf_twohalo, /*!< output only the 2-halo term */
    hmpdf_total,   /*!< output sum of both */
} hmpdf_Cell_mode_e;

/*! Returns the angular power spectrum.
 *
 *  \param[in,out] d    hmpdf_init() must have been called on d
 *  \param[in] Nell     number of ell-values the power spectrum is to be output at
 *  \param[in] ell      array of length Nell
 *  \param[out] Cell    output array, at least Nell long
 *  \param[in] mode     one of #hmpdf_Cell_mode_e
 *  \return error code
 */
int hmpdf_get_Cell(hmpdf_obj *d,
                   int Nell,
                   double ell[Nell],
                   double Cell[Nell],
                   hmpdf_Cell_mode_e mode);

/*! Returns the angular correlation function.
 *
 *  \param[in,out] d    hmpdf_init() must have been called on d
 *  \param[in] Nphi     number of phi-values the correlation function is to be output at
 *  \param[in] phi      array of length Nphi (in arcmin)
 *  \param[out] Cphi    output array, at least Nphi long
 *  \param[in] mode     one of #hmpdf_Cell_mode_e
 *  \return error code
 */
int hmpdf_get_Cphi(hmpdf_obj *d,
                   int Nphi,
                   double phi[Nphi],
                   double Cphi[Nphi],
                   hmpdf_Cell_mode_e mode);


#endif
