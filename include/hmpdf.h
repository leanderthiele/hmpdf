/*! \file hmpdf.h
 *
 *  For convenience, you can just include this single header file.
 */
#ifndef HMPDF_H
#define HMPDF_H

/*! \mainpage hmpdf Documentation
 *
 *  \section intro Introduction
 *  
 *  The hmpdf code computes one- and two-point PDFs, the covariance matrix
 *  of the one-point PDF, and the angular power spectrum/correlation function
 *  of cosmological fields.
 *  
 *  The formalism is based on the halo model, thus, only fields that are mainly
 *  sourced by virialized matter make sense.
 *  
 *  Currently, two fields are supported:
 *      + thermal Sunyaev-Zel'dovich effect (Compton-y)
 *      + weak lensing convergence
 *
 *  It code is interfaced both through C header files
 *  as well as through a convenient python wrapper.
 *
 *  \section build Building
 *
 *  Dependencies:
 *      + <a href="http://www.fftw.org/">FFTW3</a>
 *      + <a href="https://www.gnu.org/software/gsl/">GSL</a> (>= 2.4)
 *      + <a href="https://lesgourg.github.io/class_public/class.html">CLASS</a>
 *        (as of writing, the latest version has a
 *         <a href="https://github.com/lesgourg/class_public/issues/329">bug</a>
 *         that can occasionally cause segfaults; the fix is simple)
 *  
 *  Compilation should be straightforward.
 *  You may have to edit the following parts of the Makefile:
 *      + adapt PATHTOCLASS to your CLASS installation
 *      + add FFTW3/GSL locations to LINKER if those are non-standard
 *      + remove the OpenMP flag if you do not wish parallel execution
 *
 *  If you intend to use the python wrapper, you need to adapt the variable PATHTOHMPDF
 *  in hmpdf.py.
 *
 *  It may also be convenient to put libhmpdf.so into one of the default locations
 *  searched by the linker
 *  (in which case you can simply delete the PATHTOHMPDF line).
 *
 *  \section use Usage
 *
 *  The C interface can be found in hmpdf.h.
 *  All exposed objects' names start with hmpdf_ to help you keep your name space clean.
 *
 *  Workflow is as follows:
 *      1. allocate a new #hmpdf_obj with hmpdf_new().
 *      2. set all options (required and optional) with hmpdf_init(),
 *         this will also compute the data needed for all outputs.
 *         See #hmpdf_configs_e for optional inputs.
 *      3. get your output [hmpdf_get_op(), hmpdf_get_tp(), hmpdf_get_cov(),
 *                          hmpdf_get_Cell(), hmpdf_get_Cphi()].
 *      4. go to (3.) if you require any other outputs;
 *         go to (2.) if you want to re-run the code with different options.
 *      5. free the memory associated with the #hmpdf_obj with hmpdf_delete().
 *
 *  There is also a Python wrapper (hmpdf.py), with analogous workflow and
 *  simpler syntax.
 *
 *
 *  \section opt Options
 *
 *  There is a number of optional arguments to hmpdf_init(),
 *  the most important ones to get you started are
 *      + pixelization: #hmpdf_pixel_side
 *      + ell space filters: #hmpdf_tophat_radius, #hmpdf_gaussian_fwhm,
 *                           #hmpdf_custom_ell_filter (and #hmpdf_custom_ell_filter_params)
 *      + pixel-wise Gaussian noise: #hmpdf_noise
 *      + multithreading: #hmpdf_N_threads
 *
 *  See the documentation of #hmpdf_configs_e for all options.
 *
 *  \section time Runtime
 *
 *  The following are empirically found with standard settings (and only approximate!).
 *  All functions scale as #hmpdf_N_M x #hmpdf_N_z, this is omitted in the following:
 *      + hmpdf_init(): ~3 seconds on a single thread.
 *                      \par
 *                      scales as #hmpdf_N_theta to #hmpdf_N_theta^2.
 *                      Different if CLASS runtime starts to dominate.
 *                      #hmpdf_tsz takes longer than #hmpdf_kappa.
 *      + hmpdf_get_op(): fast if #hmpdf_monotonize is set to non-zero (default),
 *                        otherwise ~1 min (?)
 *                        \par
 *                        scales as #hmpdf_N_signal
 *      + hmpdf_get_tp(): ~1 min on a single thread.
 *                        \par
 *                        scales as #hmpdf_N_signal^2.
 *      + hmpdf_get_cov(): ~20 min on 40 threads.
 *                         \par
 *                         scales as #hmpdf_N_signal^2 x #hmpdf_N_phi.
 *
 *  .
 *
 *  Other functions are fast in comparison to hmpdf_init().
 *  hmpdf_init() and hmpdf_get_cov() are parallelized in critical parts,
 *  while hmpdf_get_tp() does not get faster if #hmpdf_N_threads is increased.
 *
 *  \section errrors Error handling
 *
 *  All functions [except hmpdf_new()] return a non-zero int if an error occured.
 *  A traceback will be printed to stderr.
 *
 *  \section threads Thread safety
 *
 *  All functions are *not* threadsafe: making two calls on the same #hmpdf_obj
 *  concurrently results in undefined behaviour.
 *
 *  \section examples Examples
 *
 *  Some examples are collected in examples/example.c and examples/example.py.
 *  You can compile the C code with
 *  \snippet example.c compile
 *
 *  A simple calculation of the weak lensing convergence one-point PDF
 *  with all settings at default would look like this:
 *  \snippet example.c example_kappa_onepoint
 *
 *  The same thing using the python wrapper is really easy:
 *  \snippet example.py example_kappa_onepoint
 *
 *  Including the effect of an ell-space filter would look like this:
 *  \snippet example.c example_ell_filter_use
 *  Here, we defined the function example_ell_filter conforming to the
 *  typedef #hmpdf_ell_filter_f:
 *  \snippet example.c example_ell_filter
 *
 *  We can also do this in python:
 *  \snippet example.py example_ell_filter_use
 */

#include "hmpdf_data.h"
#include "hmpdf_configs.h"
#include "hmpdf_init.h"
#include "hmpdf_onepoint.h"
#include "hmpdf_twopoint.h"
#include "hmpdf_covariance.h"
#include "hmpdf_powerspectrum.h"

#endif
