## \file
from ctypes import *
from os.path import join
import numpy as np
from numpy.ctypeslib import ndpointer, as_ctypes
from typing import Optional, Tuple
from pkg_resources import resource_filename

## \cond
class _C(object) : # char pointer
#{{{
    def __init__(self) :
        pass
    def __call__(self, name) :
        if isinstance(name, str) :
            return c_char_p(name.encode('utf-8'))
        elif isinstance(name, bytes) :
            return c_char_p(name)
        else :
            raise TypeError('Not a string.')
#}}}

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
            raise TypeError('Expected numpy array to be 1-dimensional, '\
                            'but has dimension %d'%len(arr.shape))
        if self.l is not None :
            if not arr.shape[0] == self.l :
                raise TypeError('Expected numpy array to have length %d, '\
                                'but has length %d'%(self.l, arr.shape[1]))
        return c_void_p(np.ascontiguousarray(arr).ctypes.data)
#}}}

class _F(object) : # function pointer
#{{{
    def __init__(self, nargs) :
        args = [c_double, ] * nargs
#        args.append(ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1))
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
    configs = ['N_threads', c_int, 'verbosity', c_int, 'warn_is_err', c_int,
               'class_pre', _C(), 
               'N_z', c_int, 'z_min', c_double, 'z_max', c_double,
               'N_M', c_int, 'M_min', c_double, 'M_max', c_double,
               'N_signal', c_int, 'signal_min', c_double, 'signal_max', c_double,
               'N_theta', c_int, 'rout_scale', c_double, 'rout_rdef', _E(_mdefs),
               'pixel_side', c_double, 'tophat_radius', c_double, 'gaussian_fwhm', c_double,
               'custom_ell_filter', _F(1), 'custom_ell_filter_params', None,
               'custom_k_filter', _F(2), 'custom_k_filter_params', None,
               'N_phi', c_int, 'phi_max', c_double, 'pixelexact_max', c_int,
               'phi_jitter', c_double, 'phi_pwr', c_double,
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
            raise NotImplementedError('Configuration option %s not '\
                                      'supported in the python wrapper'%key)
        pos = _Configs.configs.index(key)
        self.l.append(c_int(pos//2))
        self.t.append(_Configs.configs[pos+1])
        self.l.append(self.t[-1](value))
    def __call__(self) :
        self.l.append(c_int(len(_Configs.configs)//2)) # corresponds to end_configs
        return self.l
#}}}
## \endcond

class HMPDF(object) :
    """! Python wrapper around hmpdf.h

    The general interface is very similar to the C one,
    so please read the documentation for this.
    The differences are:
         + all names have the hmpdf_ prefix removed
         + enums are replaced by strings
         + the function signatures are a bit different
         + the argument list to init() is implemented with the **kwargs syntax
           and not ended with #hmpdf_end_configs
         + passing custom ell- and k-space filters works differently,
           see the \ref examples. The options #hmpdf_custom_ell_filter_params
           and #hmpdf_custom_k_filter_params are not supported.

    Best used in a context manager.
    """
#{{{
    # interaction with the DLL
    #{{{
    # locate the shared library
    __pathtohmpdf_name = resource_filename('hmpdf', 'PATHTOHMPDF.txt')
    with open(__pathtohmpdf_name, 'r') as f :
        PATHTOHMPDF = f.readline().rstrip()
    try :
        __libhmpdf = CDLL(join(PATHTOHMPDF, 'libhmpdf.so'))
    except OSError : # try to read the shared library from LD_LINKER_PATH
        __libhmpdf = CDLL('libhmpdf.so')
    __new = __libhmpdf.hmpdf_new
    __new.restype = POINTER(c_int)
    __delete = __libhmpdf.hmpdf_delete
    __delete.restype = c_int
    __delete.argtypes = [c_void_p, ]
    __init = __libhmpdf.hmpdf_init_fct
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
    __get_Nphi = __libhmpdf._get_Nphi
    __get_Nphi.restype = c_int
    __get_Nphi.argtypes = [c_void_p, POINTER(c_int), ]
    __get_cov_diagnostics = __libhmpdf.hmpdf_get_cov_diagnostics1
    __get_cov_diagnostics.restype = c_int
    __get_cov_diagnostics.argtypes = [c_void_p,
                                      ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1),
                                      ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1),
                                      ndpointer(c_double, flags='C_CONTIGUOUS', ndim=1), ]
    __errmsg = '%s completed with non-zero exit code'
    #}}}
    def __ret(self, err, callname, *args) :
    #{{{
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
    
    def __init__(self) -> None :
        """! calls hmpdf_new()."""
    #{{{
        self.__d = HMPDF.__new()
        if not self.__d :
            raise MemoryError('Could not allocate hmpdf_obj.')
    #}}}

    def __del__(self) -> None :
        """! calls hmpdf_delete()"""
    #{{{
        err = HMPDF.__delete(self.__d)
        if err :
            raise RuntimeError('delete failed.')
    #}}}

    def __enter__(self) -> 'HMPDF' :
        """! returns object"""
    #{{{
        return self
    #}}}

    def __exit__(self,
                 exc_type,
                 exc_value,
                 exc_traceback) -> None :
        """! calls hmpdf_delete()"""
    #{{{
        del self
    #}}}

    def init(self,
             class_ini: str,
             stype: str,
             zsource: Optional[float]=None,
             **kwargs) -> None :
        """! Initializes the object [calls hmpdf_init()].
        
        \param class_ini     CLASS .ini file
        \param stype         signal type (either "kappa" or "tsz")
        \param zsource       source redshift. Use only if stype="kappa"
        \param **kwargs      optional settings, see the documentation of #hmpdf_configs_e.
        """
    #{{{
        if stype == 'kappa' :
            if zsource is None :
                raise SyntaxError('you need to pass source redshift for stype=kappa.')
        elif stype == 'tsz' :
            if zsource is not None :
                raise SyntaxError('zsource should not be passed for stype=tsz.')
        else :  
            raise NotImplementedError('Unknown stype.')
        self.__c = _Configs()
        # add kwargs to configs
        for k,v in kwargs.items() :
            self.__c.append(k,v)
        arglist = self.__c()
        if stype == 'kappa' :
            arglist.insert(0, c_double(zsource))
        err = HMPDF.__init(self.__d, _C()(class_ini),
                           _E(_stypes)(stype), *arglist)
        return self.__ret(err, 'init()')
    #}}}

    def get_op(self,
               binedges: np.ndarray,
               incl_2h: bool=True,
               noisy: bool=False) -> np.ndarray :
        """! Get the one-point PDF [calls hmpdf_get_op()]
        
        \param binedges      1d, defines how the PDF is binned
        \param incl_2h       whether to include the two-halo term
        \param noisy         whether to include pixel-wise Gaussian noise
        \\return the binned PDF (1d)
        """
    #{{{
        out = np.empty(len(binedges)-1)
        err = HMPDF.__get_op(self.__d, len(binedges)-1, binedges, out,
                             incl_2h, noisy)
        return self.__ret(err, 'get_op()', out)
    #}}}

    def get_tp(self,
               phi: float,
               binedges: np.ndarray,
               noisy: bool=False) -> np.ndarray :
        """! Get the two-point PDF [calls hmpdf_get_tp()]
        
        \param phi           angular separation of the two sky locations (in arcmin)
        \param binedges      1d, defines how the PDF is binned
        \param noisy         whether to include pixel-wise Gaussian noise
        \\return the binned PDF (2d)
        """
    #{{{
        out = np.empty((len(binedges)-1)*(len(binedges)-1))
        err = HMPDF.__get_tp(self.__d, phi, len(binedges)-1, binedges, out, noisy)
        out = out.reshape((len(binedges)-1, len(binedges)-1))
        return self.__ret(err, 'get_tp()', out)
    #}}}

    def get_cov(self,
                binedges: np.ndarray,
                noisy: bool=False) -> np.ndarray :
        """! Get the covariance matrix of the one-point PDF [calls hmpdf_get_cov()]
        
        \param binedges      1d, defines how the covariance matrix is binned
        \param noisy         whether to include pixel-wise Gaussian noise
        \return the binned covariance matrix (2d)
        """
    #{{{
        Nbins = len(binedges) - 1
        out = np.empty(Nbins*Nbins)
        err = HMPDF.__get_cov(self.__d, Nbins, binedges, out, noisy)
        out = out.reshape((Nbins, Nbins))
        return self.__ret(err, 'get_cov()', out)
    #}}}

    def get_Cell(self,
                 ell: np.ndarray,
                 mode: str='total') -> np.ndarray :
        """! Get the angular power spectrum [calls hmpdf_get_Cell()]
        
        \param ell           1d, the angular wavenumbers
        \param mode          one of "onehalo", "twohalo", "total"
        \\return the power spectrum at ell (1d)
        """
    #{{{
        out = np.empty(len(ell))
        err = HMPDF.__get_Cell(self.__d, len(ell), ell, out,
                               _E(_corr_types)(mode))
        return self.__ret(err, 'get_Cell()', out)
    #}}}

    def get_Cphi(self,
                 phi: np.ndarray,
                 mode: str='total') -> np.ndarray :
        """! Get the angular correlation function [calls hmpdf_get_Cphi()]
        
        \param phi           1d, the angular separations (in arcmin)
        \param mode          one of "onehalo", "twohalo", "total"
        \\return the correlation function at phi (1d)
        """
    #{{{
        out = np.empty(len(phi))
        err = HMPDF.__get_Cphi(self.__d, len(phi), phi, out,
                               _E(_corr_types)(mode))
        return self.__ret(err, 'get_Cphi()', out)
    #}}}

    def get_cov_diagnostics(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """! Get the covariance diagnostics [calls hmpdf_get_cov_diagnostics()]
        
        \\return (phi, phiweights, corr_diagn)
        """
    #{{{
        Nphi = c_int(0)
        err = HMPDF.__get_Nphi(self.__d, byref(Nphi))
        if err :
            return self.__ret(err, 'get_cov_diagnostics()',
                              None, None, None)
        phi = np.empty(Nphi)
        phiweights = np.empty(Nphi)
        corr_diagn = np.empty(Nphi)
        err = HMPDF.__get_cov_diagnostics(self.__d, phi, phiweights, corr_diagn)
        return self.__ret(err, 'get_cov_diagnostics',
                          phi, phiweights, corr_diagn)
    #}}}
