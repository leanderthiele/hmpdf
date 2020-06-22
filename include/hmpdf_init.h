/*! \file hmpdf_init.h */
#ifndef HMPDF_INIT_H
#define HMPDF_INIT_H

/*! Initializes the #hmpdf_obj and computes data that is needed for all outputs.
 *
 *  \param d            created with hmpdf_new()
 *  \param class_ini    path to a CLASS .ini file
 *  \param stype        signal type (either #hmpdf_kappa or #hmpdf_tsz)
 *  \param ...          variable argument list for optional arguments
 *  \return void
 *  
 *  \attention if stype=#hmpdf_kappa, the first entry in ... must be the source
 *             redshift (double)
 *  \attention the last argument in ... must be #hmpdf_end_configs, regardless
 *             of whether ... is empty otherwise.
 *
 *  The syntax to pass additional non-default settings through ... is as follows:
 *  A setting is passed as a pair <name>, <val>, where <name> is one of
 *  #hmpdf_configs_e, and <val> must have the type given in the documentation
 *  for that specific <name>.
 *  
 *  For example, to perform a weak lensing calculation with source redshift 1
 *  and all configurations at default except for the number of threads and the
 *  pixel sidelength, you would call
 *  \code
 *  hmpdf_init(d, "explanatory.ini", kappa, 1.0,
 *             hmpdf_N_threads, 4,    // use 4 threads
 *             hmpdf_pixel_side, 0.4, // use 0.4 arcmin pixels
 *             hmpdf_end_configs);
 *  \endcode
 *
 *  \attention perhaps counter-intuitively, successive hmpdf_init() calls are not
 *             cumulative (i.e. d retains no internal state between them).
 *             Thus, you need to pass all configuration options in a single call
 *             to hmpdf_init().
 */
int hmpdf_init(hmpdf_obj *d,
               char *class_ini,
               hmpdf_signaltype_e stype,
               ...);

#endif
