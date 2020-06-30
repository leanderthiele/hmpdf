/*! \file hmpdf_configs.h */
#ifndef HMPDF_CONFIGS_H
#define HMPDF_CONFIGS_H

/*! Halo mass definitions.
 *  
 *  Presently only used to specify in terms of which
 *  definition the radial cut-off of the halo profiles
 *  is given.
 */
typedef enum
{
    hmpdf_mdef_c, /*!< M200c (200x critial density) */
    hmpdf_mdef_v, /*!< Mvir (virial mass according to Bryan+Norman 1998) */
    hmpdf_mdef_m, /*!< M200m (200x mean matter density) */
} hmpdf_mdef_e;

/*! Signal types.
 *
 *  To specify which cosmological field the PDF should be computed for.
 */
typedef enum
{
    hmpdf_kappa, /*!< weak lensing convergence */
    hmpdf_tsz,   /*!< tSZ effect (Compton-y) */
} hmpdf_signaltype_e;

/*! Fixed point integration modes.
 *
 *  The integrals over halo mass and redshift are performed using fixed point quadratures,
 *  with modes taken from this enum.
 *  See the
 *  <a href="https://www.gnu.org/software/gsl/doc/html/integration.html#fixed-point-quadratures">
 *  GSL documentation</a>
 *  for details.
 *
 *  \attention Some of these integration modes are for infinite intervals and do not make sense
 *             for our application.
 */
typedef enum
{
    hmpdf_legendre, /*!< .*/
    hmpdf_chebyshev, /*!< .*/
    hmpdf_gegenbauer, /*!< .*/
    hmpdf_jacobi, /*!< .*/
    hmpdf_laguerre, /*!< .*/
    hmpdf_hermite, /*!< .*/
    hmpdf_exponential, /*!< .*/
    hmpdf_rational, /*!< .*/
    hmpdf_chebyshev2, /*!< .*/
} hmpdf_integr_mode_e;

/*! Options to hmpdf_init().
 *
 *  The variable argument list in hmpdf_init() can be used to pass non-default options.
 *  [syntax is explained in the documentation for hmpdf_init()].
 *
 *  There is a large number of options, many of which the typical user will not need to use.
 *  Here is a brief synopsis:
 *
 *  Always used: #hmpdf_end_configs
 *
 *  Frequently used options:
 *      + pixelization: #hmpdf_pixel_side
 *      + ell space filters: #hmpdf_tophat_radius, #hmpdf_gaussian_fwhm,
 *                           #hmpdf_custom_ell_filter (and #hmpdf_custom_ell_filter_params)
 *      + pixel-wise Gaussian noise: #hmpdf_noise
 *      + multithreading: #hmpdf_N_threads
 *  
 *  Less frequently used options:
 *      + verbosity: #hmpdf_verbosity
 *      + halo model fit parameters: #hmpdf_Duffy08_conc_params,
 *                                   #hmpdf_Tinker10_hmf_params,
 *                                   #hmpdf_Battaglia12_tsz_params
 *      + k space filter: #hmpdf_custom_k_filter (and #hmpdf_custom_k_filter_params)
 *      + monotonization of halo profiles: #hmpdf_monotonize
 *      + PDF internal sampling points: #hmpdf_N_signal, #hmpdf_signal_min, #hmpdf_signal_max
 *  
 *  Integration grids:
 *      + redshift integration: #hmpdf_N_z, #hmpdf_z_min, #hmpdf_z_max,
 *                              #hmpdf_zintegr_type, #hmpdf_zintegr_alpha, #hmpdf_zintegr_beta
 *      + halo mass integration: #hmpdf_N_M, #hmpdf_M_min, #hmpdf_M_max,
 *                               #hmpdf_Mintegr_type, #hmpdf_Mintegr_alpha, #hmpdf_Mintegr_beta
 *      + halo profile angular integration: #hmpdf_N_theta, #hmpdf_rout_scale, #hmpdf_rout_rdef
 *
 *  Covariance matrix calculation:
 *      + useful to improve numerical stability: #hmpdf_N_phi
 *      + integration/summation grid: #hmpdf_phi_max, #hmpdf_pixelexact_max, #hmpdf_phi_jitter,
 *                                    #hmpdf_phi_pwr
 */
typedef enum
{
    hmpdf_N_threads, /*!< number of threads to use in multithreaded parts of the code.
                    *   \par
                    *   Type: int. Default: 1.
                    *   \remark only applicable if code compiled with OpenMP.
                    *   \remark this does *not* control the multithreading behaviour of CLASS
                    */
    hmpdf_verbosity, /*!< larger values yield more detailed print output.
                      *   If 0 (default), only error messages are printed.
                      *   \par
                      *   Type: int. Default: 0.
                      */
    hmpdf_class_pre, /*!< optional CLASS precision file.
                      *   \par
                      *   Type: char *. Default: None (use CLASS's default precision settings).
                      */
    hmpdf_N_z, /*!< number of sample points in redshift integration.
                *   \par
                *   Type: int. Default: 65.
                *   \remark the default value is conservative.
                *           By playing with the #hmpdf_zintegr_type
                *           (as well as #hmpdf_zintegr_alpha and #hmpdf_zintegr_beta)
                *           you can gain some speed here.
                */
    hmpdf_z_min, /*!< minimum redshift.
                  *   \par
                  *   Type: double. Default: 0.
                  *   \remark non-zero values do not make much sense
                  */
    hmpdf_z_max, /*!< maximum redshift.
                  *   \par
                  *   Type: double. Default: 6 for tSZ, source redshift for weak lensing.
                  *   \warning If set to values larger than source redshift in weak lensing,
                  *            behaviour is undefined.
                  */
    hmpdf_N_M, /*!< number of sample points in halo mass integration.
                *   \par
                *   Type: int. Default: 65.
                *   \remark the default value is conservative.
                *           By playing with the #hmpdf_Mintegr_type
                *           (as well as #hmpdf_Mintegr_alpha and #hmpdf_Mintegr_beta)
                *           you can gain some speed here.
                */
    hmpdf_M_min, /*!< minimum halo mass, in M200m definition and Msun units.
                  *   \par
                  *   Type: double. Default: 1e11.
                  */
    hmpdf_M_max, /*!< maximum halo mass, in M200m definition and Msun units.
                  *   \par
                  *   Type: double. Default: 1e16.
                  */
    hmpdf_N_signal, /*!< number of points on which the one-point PDF is sampled internally
                     *   \par
                     *   Type: int. Default: 1024.
                     *   \remark two-point PDF will be sampled on the square of this.
                     *   \remark Should be chosen as a power of 2 for optimal speed.
                     *   \warning Setting this to large values can trigger numerical instability
                     *            in small-phi two-point PDF, and covariance matrix computations.
                     */
    hmpdf_signal_min, /*!< minimum signal value at which PDF is sampled internally.
                       *   \par
                       *   Type: double. Default: 0.
                       *   \warning non-zero values are not tested and will likely give wrong results.
                       */
    hmpdf_signal_max, /*!< maximum signal value at which PDF is sampled internally.
                       *   \par
                       *   Type: double. Default: 2e-4 (tSZ), 1 (weak lensing).
                       *   \remark Should be chosen such that PDF is well converged to zero to avoid ringing.
                       *   \remark Should be at least 2x larger than maximum signal value you are interested in.
                       *   \remark For weak lensing, you should change this option depending on source redshift.
                       */
    hmpdf_N_theta, /*!< number of points on which the signal profiles are sampled.
                    *   \par
                    *   Type: int. Default: 500.
                    *   \remark reasonable precision can be reached with as few as 100 sample points.
                    *           Runtime of hmpdf_init() is in many cases quite dominated by this setting.
                    */
    hmpdf_rout_scale, /*!< radial cut-off.
                       *   \par
                       *   Type: double. Default: 2.
                       */
    hmpdf_rout_rdef, /*!< mass/radius definition in terms of which the radial cut-off is specified.
                      *   \par
                      *   Type: #hmpdf_mdef_e. Default: #hmpdf_mdef_v (virial radius).
                      */
    hmpdf_pixel_side, /*!< pixel sidelength, in arcmin.
                       *   \par
                       *   Type: double. Default: None.
                       *   \remark required setting for covariance matrix calculation.
                       *   \remark negative values have no effect (and will not trigger a warning).
                       *   \warning PDFs computed with more than 5 arcmin in this setting should be treated
                       *            with caution.
                       */
    hmpdf_tophat_radius, /*!< includes the effect of smoothing the map with a tophat of given radius (in arcmin).
                          *   \par
                          *   Type: double. Default: None.
                          *   \remark negative values have no effect (and will not trigger a warning).
                          *   \warning PDFs computed with more than 5 arcmin in this setting should be treated
                          *            with caution.
                          */
    hmpdf_gaussian_fwhm, /*!< includes the effect of smoothing the map with a Gaussian of given FWHM (in arcmin).
                          *   \par
                          *   Type: double. Default: None.
                          *   \remark negative values have no effect (and will not trigger a warning).
                          *   \warning PDFs computed with more than 5 arcmin in this setting should be treated
                          *            with caution.
                          */
    hmpdf_custom_ell_filter, /*!< pass a user-defined ell-space filter
                              *   (for example, to include effect of a Wiener filter).
                              *   Signature of this function pointer has to conform
                              *   to the typedef #hmpdf_ell_filter_f.
                              *   \par
                              *   Type: #hmpdf_ell_filter_f. Default: None.
                              */
    hmpdf_custom_ell_filter_params, /*!< pass parameters to the ell-space filter (as its last argument).
                                     *   \par
                                     *   Type: void *. Default: None.
                                     *   \attention not supported in the python wrapper.
                                     */
    hmpdf_custom_k_filter, /*!< pass a user-defined k-space filter, with possible redshift-dependence
                            *   (for example, to emulate small-scale simulation resolution issues).
                            *   Signature of this function pointer has to conform
                            *   to the typedef #hmpdf_k_filter_f.
                            *   \par
                            *   Type: #hmpdf_k_filter_f. Default: None.
                            */
    hmpdf_custom_k_filter_params, /*!< pass parameters to the k-space filter (as its last argument).
                                   *   \par
                                   *   Type: void *. Default: None.
                                   *   \attention not supported in the python wrapper.
                                   */
    hmpdf_N_phi, /*!< Number of pixel-separation sample points in covariance matrix calculation.
                  *   \par
                  *   Type: int. Default: 1000.
                  *   \remark as long as no numerical instability at small pixel separations is encountered,
                  *           the default value is conservative.
                  *           In that case, if you want to decrease run-time, it is recommended you decrease
                  *           this value to a few hundred, and decrease #hmpdf_pixelexact_max as well.
                  *   \remark if, however, numerical instability is encountered,
                  *           results can get more accurate by increasing this value.
                  */
    hmpdf_phi_max, /*!< maximum pixel separation in covariance matrix calculation (in arcmin).
                    *   \par
                    *   Type: double. Default: 150.
                    *   \remark should be set to approximately the largest halo radius on the sky.
                    */
    hmpdf_pixelexact_max, /*!< up to which pixel separation (in units of pixel sidelength) the summation
                           *   is performed exactly (integration afterwards).
                           *   \par
                           *   Type: int. Default: 20.
                           *   \remark the default value is conservative, can be decreased to 10 without
                           *           losing much accuracy.
                           */
    hmpdf_phi_jitter, /*!< technical (for covariance matrix).
                       *   For numerical stability, at small pixel separations an average is used,
                       *   around the desired pixel separation.
                       *   This option sets how wide (in phi) the averaging is done.
                       *   \par
                       *   Type: double. Default: 0.02.
                       */
    hmpdf_phi_pwr, /*!< technical (for covariance matrix).
                    *   Increasing this number increases the sample points at small pixel separation.
                    *   \par
                    *   Type: double. Default: 2.
                    *   \remark If you encounter numerical instability (e.g., wild covariance matrix entries),
                    *           increasing this number can potentially be helpful.
                    */
    hmpdf_monotonize, /*!< If non-zero (default), non-monotonic sections in the signal profiles will be patched.
                       *   \par
                       *   Type: int. Default: 1.
                       *   \remark as long as you do not apply any custom ell- or k-space filter to the signal
                       *           profiles, leaving this as 1 is very safe.
                       *   \remark if you do use a custom filter, the rule of thumb is that as long as the filter
                       *           is close to 1 for small ell/k, it should be fine as well.
                       *   \remark this has to be set to 1 for two-point PDF and covariance matrix calculations
                       *           (otherwise the integrations are impossibly slow).
                       */
    hmpdf_zintegr_type, /*!< Fixed point integration mode for redshift integration.
                         *   \par
                         *   Type: #hmpdf_integr_mode_e. Default: #hmpdf_legendre.
                         */
    hmpdf_zintegr_alpha, /*!< Fixed point integration alpha for redshift integration
                          *   \par
                          *   Type: double. Default: 0.
                          */
    hmpdf_zintegr_beta, /*!< Fixed point integration beta for redshift integration.
                         *   \par
                         *   Type: double. Default: 0.
                         */
    hmpdf_Mintegr_type, /*!< Fixed point integration mode for halo mass integration.
                         *   \par
                         *   Type: #hmpdf_integr_mode_e. Default: #hmpdf_legendre.
                         */
    hmpdf_Mintegr_alpha, /*!< Fixed point integration alpha for halo mass integration.
                          *   \par
                          *   Type: double. Default: 0.
                          */
    hmpdf_Mintegr_beta, /*!< Fixed point integration beta for halo mass integration.
                         *   \par
                         *   Type: double. Default: 0.
                         */
    hmpdf_Duffy08_conc_params, /*!< Fit parameters in the
                                *   <a href="https://arxiv.org/abs/0804.2486">Duffy+2008</a>
                                *   concentration model.
                                *   \par
                                *   Type: double[9]. Default: see src/configs.c.
                                *   \remark in the python wrapper, pass a 1d numpy array
                                */
    hmpdf_Tinker10_hmf_params, /*!< Fit parameters in the
                                *   <a href="https://arxiv.org/abs/1001.3162">Tinker+2010</a>
                                *   halo mass function.
                                *   \par
                                *   Type: double[10]. Default: see src/configs.c.
                                *   \remark in the python wrapper, pass a 1d numpy array
                                */
    hmpdf_Battaglia12_tsz_params, /*!< Fit parameters in the
                                   *   <a href="https://arxiv.org/abs/1109.3711">Battaglia+2012</a>
                                   *   pressure profile model.
                                   *   \par
                                   *   Type: double[15]. Default: see src/configs.c.
                                   *   \remark in the python wrapper, pass a 1d numpy array
                                   */
    hmpdf_noise, /*!< Option to add pixel-wise Gaussian noise of this standard deviation.
                  *   \par
                  *   Type: double. Default: None.
                  *   \remark negative values have no effect (and will not trigger a warning).
                  */
    hmpdf_end_configs, /*!< \attention Required last argument in hmpdf_init(). */
    // keep this last
} hmpdf_configs_e;

/*! Function pointer typedef for user-defined ell-space filter.
 *  Passed to hmpdf_init() as #hmpdf_custom_ell_filter.
 *  \param ell      angular wavenumber
 *  \param p        pointer that allows user to pass other parameters.
 *                  Passed to hmpdf_init() as #hmpdf_custom_ell_filter_params.
 *  \return W(ell)  window function at ell
 */
typedef double (*hmpdf_ell_filter_f)(double,
                                     void *);

/*! Function pointer typedef for user-defined k-space filter.
 *  Passed to hmpdf_init() as #hmpdf_custom_k_filter.
 *  \param k        comoving wavenumber in 1/Mpc
 *  \param z        redshift
 *  \param p        pointer that allows user to pass other parameters.
 *                  Passed to hmpdf_init() as #hmpdf_custom_k_filter_params.
 *  \return W(k,z)  window function at k and z
 */
typedef double (*hmpdf_k_filter_f)(double,
                                   double,
                                   void *);

#endif
