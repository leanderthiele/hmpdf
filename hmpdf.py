## \file
from ctypes import *
from os.path import join
import numpy as np
from numpy.ctypeslib import ndpointer, as_ctypes

## \cond
PATHTOHMPDF = '/home/leander/Perimeter/Onepoint/C_implementation'

class _E(object) : # enumeration class
#{{{
    def __init__(self, names) :
        self.names = names
    def __call__(self, name) :
        return c_int(self.names.index(name))
#}}}

class _D(object) : # double pointer from numpy array
#{{{
    def __init__(self, l) :
        self.l = l
    def __call__(self, arr) :
        if not len(arr.shape) == 1 :
            raise TypeError('Expected numpy array to be 1-dimensional, \
                             but has dimension %d'%len(arr.shape))
        if self.l is not None :
            if not arr.shape[0] == self.l :
                raise TypeError('Expected numpy array to have length %d, \
                                 but has length %d'%(self.l, arr.shape[1]))
        return c_void_p(np.ascontiguousarray(arr).ctypes.data)
#}}}

class _F(object) : # function pointer
#{{{
    def __init__(self, nargs) :
        args = [c_double, ] * nargs
        args.append(ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1))
        self.fctptr = CFUNCTYPE(c_double, *args)
    def __call__(self, f) :
        return self.fctptr(f)
#}}}

# Enumerations
#{{{
_integr = ['legendre', 'chebyshev', 'gegenbauer', 'jacobi', 'laguerre',
          'hermite', 'exponential', 'rational', 'chebyshev2']
_stypes = ['kappa', 'tsz']
_mdefs = ['mdef_c', 'mdef_v', 'mdef_m', ]
_corr_types = ['onehalo', 'twohalo', 'total']
#}}}

class _Configs(object) :
#{{{
    configs = ['N_threads', c_int, 'verbosity', c_int, 'class_pre', c_char_p,
               'N_z', c_int, 'z_min', c_double, 'z_max', c_double,
               'N_M', c_int, 'M_min', c_double, 'M_max', c_double,
               'N_signal', c_int, 'signal_min', c_double, 'signal_max', c_double,
               'N_theta', c_int, 'rout_scale', c_double, 'rout_rdef', _E(_mdefs),
               'pixel_side', c_double, 'tophat_radius', c_double, 'gaussian_fwhm', c_double,
               'custom_ell_filter', _F(1), 'custom_ell_filter_params', _D(None),
               'custom_k_filter', _F(2), 'custom_k_filter_params', _D(None),
               'N_phi', c_int, 'phi_max', c_double, 'pixelexact_max', c_int,
               'phi_jitter', c_double, 'phi_pwr', c_double,
               'monotonize', c_int,
               'zintegr_type', _E(_integr), 'zintegr_alpha', c_double, 'zintegr_beta', c_double,
               'Mintegr_type', _E(_integr), 'Mintegr_alpha', c_double, 'Mintegr_beta', c_double,
               'Duffy08_conc_params', _D(9),
               'Tinker10_hmf_params', _D(10),
               'Battaglia12_tsz_params', _D(15),
               'noise', c_double, ]
    def __init__(self) :
        self.t = [] # keep references to the types
        self.l = []
    def append(self, key, value) :
        if key not in _Configs.configs :
            raise KeyError('Invalid configuration option %s.'%key)
        if _Configs.configs[_Configs.configs.index(key)+1] is None :
            raise NotImplementedError('Configuration option %s not \
                                       supported in the python wrapper'%key)
        pos = _Configs.configs.index(key)
        self.l.append(c_int(pos/2))
        self.t.append(_Configs.configs[pos+1])
        self.l.append(self.t[-1](value))
    def __call__(self) :
        self.l.append(c_int(len(_Configs.configs)/2)) # corresponds to end_configs
        return self.l
#}}}
## \endcond

## Python wrapper around hmpdf.h
#
#  The general interface is very similar to the C one,
#  so please read the documentation for this.
#  The differences are:
#       + all names have the hmpdf_ prefix removed
#       + the function signatures are a bit different
#       + the argument list to init() is not ended with #hmpdf_end_configs
#       + you have the option to let the object do the error handling
#       + the function signatures for the options #hmpdf_custom_ell_filter
#         and #hmpdf_custom_k_filter are restricted, see the \ref examples.
#
#  Best used in a context manager.
class HMPDF(object) :
#{{{
    # interaction with the DLL
    #{{{
    __libhmpdf = CDLL(join(PATHTOHMPDF, 'libhmpdf.so'))
    __new = __libhmpdf.hmpdf_new
    __new.restype = POINTER(c_int)
    __delete = __libhmpdf.hmpdf_delete
    __delete.restype = c_int
    __delete.argtypes = [c_void_p, ]
    __init = __libhmpdf.hmpdf_init
    __init.restype = c_int
    __get_op = __libhmpdf.hmpdf_get_op
    __get_op.restype = c_int
    __get_op.argtypes = [c_void_p, c_int,
                         ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1),
                         ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1),
                         c_int, c_int, ]
    __get_tp = __libhmpdf.hmpdf_get_tp
    __get_tp.restype = c_int
    __get_tp.argtypes = [c_void_p, c_double, c_int,
                         ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1),
                         ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1),
                         c_int, ]
    __get_cov = __libhmpdf.hmpdf_get_cov
    __get_cov.restype = c_int
    __get_cov.argtypes = [c_void_p, c_int,
                          ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1),
                          ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1),
                          c_int, ]
    __get_Cell = __libhmpdf.hmpdf_get_Cell
    __get_Cell.restype = c_int
    __get_Cell.argtypes = [c_void_p, c_int,
                           ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1),
                           ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1),
                           c_int, ]
    __get_Cphi = __libhmpdf.hmpdf_get_Cphi
    __get_Cphi.restype = c_int
    __get_Cphi.argtypes = [c_void_p, c_int,
                           ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1),
                           ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1),
                           c_int, ]
    __errmsg = '%s completed with non-zero exit code'
    #}}}
    def __ret(self, err, callname, *args) :
    #{{{
        if not self.catch_errors :
            return (err, ) + args
        else :
            if err :
                raise RuntimeError(HMPDF.__errmsg%callname)
            else :
                if len(args) == 0 :
                    return None
                elif len(args) == 1 :
                    return args[0]
                else :
                    return args
    #}}}
    
    ## calls hmpdf_new().
    def __init__(self) :
    #{{{
        self.d = HMPDF.__new()
    #}}}
    
    ## Initializes the object [calls hmpdf_init()].
    #
    #  \param class_ini     string, the CLASS .ini file
    #  \param stype         string, the signal type (either "kappa" or "tsz")
    #  \param *args         used if stype="kappa" to pass the source redshift
    #  \param **kwargs      optional settings, see the documentation of #hmpdf_configs_e.
    #                       \note there is one additional setting, "catch_errors".
    #                             If set to false, the python code does not perform error
    #                             checking and outputs the error code as the first entry
    #                             in any output tuple.
    #                             Default is True, i.e. the python code will raise errors
    #                             upon encountering a non-zero error code.
    #  \return the error code if catch_errors=True, else None
    def init(self, class_ini, stype, *args, **kwargs) :
    #{{{
        if stype == 'kappa' :
            if len(args) != 1 :
                raise SyntaxError('*args not of appropriate length for stype=kappa, '\
                                  'you need to pass source redshift.')
        elif stype == 'tsz' :
            if len(args) != 0 :
                raise SyntaxError('*args is not empty, which it should be for stype=tsz')
        else :  
            raise NotImplementedError('Unknown stype.')
        self.catch_errors = True
        self.c = _Configs()
        # add kwargs to configs
        for k,v in kwargs.items() :
            if k == 'catch_errors' :
                self.catch_errors = v
            else :
                self.c.append(k,v)
        arglist = self.c()
        if stype == 'kappa' :
            arglist.insert(0, c_double(args[0]))
        err = HMPDF.__init(self.d, c_char_p(class_ini),
                           _E(_stypes)(stype), *arglist)
        return self.__ret(err, 'init()')
    #}}}

    ## calls hmpdf_delete()
    def _delete(self) :
    #{{{
        err = HMPDF.__delete(self.d)
        return self.__ret(err, 'delete()')
    #}}}

    ## returns object
    def __enter__(self) :
    #{{{
        return self
    #}}}

    ## calls hmpdf_delete()
    def __exit__(self, exc_type, exc_value, exc_traceback) :
    #{{{
        self._delete()
    #}}}

    ## Get the one-point PDF [calls hmpdf_get_op()]
    #
    #  \param binedges      a 1d numpy array
    #  \param **kwargs      to pass optional arguments:
    #                           + incl_2h: default True
    #                           + noisy: default False
    #  \return the binned PDF (1d numpy array)
    def get_op(self, binedges, **kwargs) :
    #{{{
        incl_2h = kwargs['incl_2h'] if 'incl_2h' in kwargs else True
        noisy = kwargs['noisy'] if 'noisy' in kwargs else False
        out = np.empty(len(binedges)-1)
        err = HMPDF.__get_op(self.d, len(binedges)-1, binedges, out,
                             incl_2h, noisy)
        return self.__ret(err, 'get_op()', out)
    #}}}

    ## Get the two-point PDF [calls hmpdf_get_tp()]
    #
    #  \param binedges      a 1d numpy array
    #  \param phi           a float
    #  \param **kwargs      to pass optional arguments:
    #                           + noisy: default False
    #  \return the binned PDF (2d numpy array)
    def get_tp(self, phi, binedges, **kwargs) :
    #{{{
        noisy = kwargs['noisy'] if 'noisy' in kwargs else False
        out = np.empty((len(binedges)-1)*(len(binedges)-1))
        HMPDF.__get_tp(self.d, phi, len(binedges)-1, binedges, out, noisy)
        out = out.reshape((len(binedges)-1, len(binedges)-1))
        return self.__ret(err, 'get_tp()', out)
    #}}}

    ## Get the covariance matrix of the one-point PDF [calls hmpdf_get_cov()]
    #
    #  \param binedges      a 1d numpy array
    #  \param **kwargs      to pass optional arguments:
    #                           + noisy: default False
    def get_cov(self, binedges, **kwargs) :
    #{{{
        noisy = kwargs['noisy'] if 'noisy' in kwargs else False
        Nbins = len(binedges) - 1
        out = np.empty(Nbins*Nbins)
        err = HMPDF.__get_cov(self.d, Nbins, binedges, out, fname)
        out = out.reshape((Nbins, Nbins))
        return self.__ret(err, 'get_cov()', out)
    #}}}

    ## Get the angular power spectrum [calls hmpdf_get_Cell()]
    #
    #  \param ell           a 1d numpy array
    #  \param **kwargs      to pass optional arguments
    #                           + mode: default total
    def get_Cell(self, ell, **kwargs) :
    #{{{
        mode = kwargs['mode'] if 'mode' in kwargs else 'total'
        out = np.empty(len(ell))
        err = HMPDF.__get_Cell(self.d, len(ell), ell, out,
                               _E(_corr_types)(mode))
        return self.__ret(err, 'get_Cell()', out)
    #}}}

    ## Get the angular correlation function [calls hmpdf_get_Cphi()]
    #
    #  \param phi           a 1d numpy array
    #  \param **kwargs      to pass optional arguments
    #                           + mode: default total
    def get_Cphi(self, phi, **kwargs) :
    #{{{
        mode = kwargs['mode'] if 'mode' in kwargs else 'total'
        out = np.empty(len(phi))
        err = HMPDF.__get_Cphi(self.d, len(phi), phi, out,
                               _E(_corr_types)(mode))
        return self.__ret(err, 'get_Cphi()', out)
    #}}}
#}}}
