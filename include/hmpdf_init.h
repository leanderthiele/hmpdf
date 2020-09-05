/*! \file hmpdf_init.h */
#ifndef HMPDF_INIT_H
#define HMPDF_INIT_H

#include "hmpdf_object.h"
#include "hmpdf_configs.h"

/*! Initializes the #hmpdf_obj and computes data that is needed for all outputs.
 *
 *  \param[in,out] d        created with hmpdf_new()
 *  \param[in] class_ini    path to a CLASS .ini file
 *  \param[in] stype        signal type (either #hmpdf_kappa or #hmpdf_tsz)
 *  \param[in] ...          variable argument list for optional arguments
 *  \return error code
 *  
 *  \attention if stype=#hmpdf_kappa, the first entry in ... must be the source
 *             redshift (double)
 *  \remark    this function performs some basic sanity checks on the user inputs.
 *             If inputs fall out of recommended bounds,
 *             an error message will be printed.
 *             However, execution will continue (for the unlikely case the input was
 *                                               what you intended it to be).
 *
 *  The syntax to pass additional non-default settings through ... is as follows:
 *  A setting is passed as a pair \<name\>, \<value\>, where \<name\> is one of
 *  #hmpdf_configs_e, and \<value\> must have the type given in the documentation
 *  for that specific \<name\>.
 *  \attention In this context, be careful with literals passed as \<value\>,
 *             as the compiler will be unable to figure out any necessary casts.
 *             For example, if passing a literal as #hmpdf_N_signal,
 *             it needs to be an explicit long.
 *  
 *  For example, to perform a weak lensing calculation with source redshift 1
 *  and all configurations at default except for the number of threads and the
 *  pixel sidelength, you would call
 *  \code
 *  hmpdf_init(d, "example.ini", hmpdf_kappa, 1.0,
 *             hmpdf_N_threads, 4,
 *             hmpdf_pixel_side, 0.4);
 *  \endcode
 *
 *  \attention perhaps counter-intuitively, successive hmpdf_init() calls are not
 *             cumulative (i.e. d retains no internal state between them).
 *             Thus, you need to pass all configuration options in a single call
 *             to hmpdf_init().
 */
#define hmpdf_init(...) hmpdf_init_fct(__VA_ARGS__, hmpdf_end_configs)

/*! Underlying function. For convenience, it should be called through the hmpdf_init macro.
 */
int hmpdf_init_fct(hmpdf_obj *d,
                   char *class_ini,
                   hmpdf_signaltype_e stype,
                   ...);

#endif
