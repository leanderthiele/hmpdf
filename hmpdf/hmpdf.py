## \file
import ctypes as ct
from os.path import join
import numpy as np
from numpy.ctypeslib import ndpointer, as_ctypes
from typing import Optional, Tuple, Union, Sequence
from pkg_resources import resource_filename

## \cond
class _C(object) : # char pointer
#{{{
    def __init__(self) :
        pass
    def __call__(self, name) :
        if isinstance(name, str) :
            return ct.c_char_p(name.encode('utf-8'))
        elif isinstance(name, bytes) :
            return ct.c_char_p(name)
        else :
            raise TypeError('Not a string.')
#}}}

class _E(object) : # enumeration class
#{{{
    integr = ['legendre', 'chebyshev', 'gegenbauer', 'jacobi', 'laguerre',
              'hermite', 'exponential', 'rational', 'chebyshev2']
    stypes = ['kappa', 'tsz']
    mdefs = ['mdef_c', 'mdef_v', 'mdef_m', ]
    corr_types = ['onehalo', 'twohalo', 'total']
    def __init__(self, names) :
        self.names = names
    def __call__(self, name) :
        return ct.c_int(self.names.index(name))
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
        return ct.c_void_p(np.ascontiguousarray(arr).ctypes.data)
#}}}

class _F(object) : # function pointer
#{{{
    def __init__(self, nargs) :
        args = [ct.c_double, ] * nargs
        self.fctptr = ct.CFUNCTYPE(ct.c_double, *args)
    def __call__(self, f) :
        return self.fctptr(f)
#}}}

class _Configs(object) :
#{{{
    configs = ['N_threads', ct.c_int, 'verbosity', ct.c_int, 'warn_is_err', ct.c_int,
               'class_pre', _C(), 
               'N_z', ct.c_int, 'z_min', ct.c_double, 'z_max', ct.c_double,
               'N_M', ct.c_int, 'M_min', ct.c_double, 'M_max', ct.c_double,
               'N_signal', ct.c_long, 'signal_min', ct.c_double, 'signal_max', ct.c_double,
               'N_theta', ct.c_int, 'rout_scale', ct.c_double, 'rout_rdef', _E(_E.mdefs),
               'pixel_side', ct.c_double, 'tophat_radius', ct.c_double, 'gaussian_fwhm', ct.c_double,
               'custom_ell_filter', _F(1), 'custom_ell_filter_params', None,
               'custom_k_filter', _F(2), 'custom_k_filter_params', None,
               'N_phi', ct.c_int, 'phi_max', ct.c_double, 'pixelexact_max', ct.c_int,
               'phi_jitter', ct.c_double, 'phi_pwr', ct.c_double,
               'zintegr_type', _E(_E.integr), 'zintegr_alpha', ct.c_double, 'zintegr_beta', ct.c_double,
               'Mintegr_type', _E(_E.integr), 'Mintegr_alpha', ct.c_double, 'Mintegr_beta', ct.c_double,
               'Duffy08_conc_params', _D(9), 'Tinker10_hmf_params', _D(10), 'Battaglia12_tsz_params', _D(15),
               'noise_pwr', _F(1), 'noise_pwr_params', None, ]
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
        self.l.append(ct.c_int(pos//2))
        self.t.append(_Configs.configs[pos+1])
        self.l.append(self.t[-1](value))
    def __call__(self) :
        self.l.append(ct.c_int(len(_Configs.configs)//2)) # corresponds to end_configs
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
           see the @ref examples. The options #hmpdf_custom_ell_filter_params
           and #hmpdf_custom_k_filter_params are not supported.
           Analogous for the noise power spectrum (and #hmpdf_noise_pwr_params).

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
        __libhmpdf = ct.CDLL(join(PATHTOHMPDF, 'libhmpdf.so'))
    except OSError : # try to read the shared library from LD_LIBRARY_PATH
        __libhmpdf = ct.CDLL('libhmpdf.so')
    __new = __libhmpdf.hmpdf_new
    __new.restype = ct.POINTER(ct.c_int)
    __delete = __libhmpdf.hmpdf_delete
    __delete.restype = ct.c_int
    __delete.argtypes = [ct.c_void_p, ]
    __init = __libhmpdf.hmpdf_init_fct
    __init.restype = ct.c_int
    __get_op = __libhmpdf.hmpdf_get_op
    __get_op.restype = ct.c_int
    __get_op.argtypes = [ct.c_void_p, ct.c_int,
                         ndpointer(ct.c_double, flags='C_CONTIGUOUS', ndim=1),
                         ndpointer(ct.c_double, flags='C_CONTIGUOUS', ndim=1),
                         ct.c_int, ct.c_int, ]
    __get_tp = __libhmpdf.hmpdf_get_tp
    __get_tp.restype = ct.c_int
    __get_tp.argtypes = [ct.c_void_p, ct.c_double, ct.c_int,
                         ndpointer(ct.c_double, flags='C_CONTIGUOUS', ndim=1),
                         ndpointer(ct.c_double, flags='C_CONTIGUOUS', ndim=1),
                         ct.c_int, ]
    __get_cov = __libhmpdf.hmpdf_get_cov
    __get_cov.restype = ct.c_int
    __get_cov.argtypes = [ct.c_void_p, ct.c_int,
                          ndpointer(ct.c_double, flags='C_CONTIGUOUS', ndim=1),
                          ndpointer(ct.c_double, flags='C_CONTIGUOUS', ndim=1),
                          ct.c_int, ]
    __get_Cell = __libhmpdf.hmpdf_get_Cell
    __get_Cell.restype = ct.c_int
    __get_Cell.argtypes = [ct.c_void_p, ct.c_int,
                           ndpointer(ct.c_double, flags='C_CONTIGUOUS', ndim=1),
                           ndpointer(ct.c_double, flags='C_CONTIGUOUS', ndim=1),
                           ct.c_int, ]
    __get_Cphi = __libhmpdf.hmpdf_get_Cphi
    __get_Cphi.restype = ct.c_int
    __get_Cphi.argtypes = [ct.c_void_p, ct.c_int,
                           ndpointer(ct.c_double, flags='C_CONTIGUOUS', ndim=1),
                           ndpointer(ct.c_double, flags='C_CONTIGUOUS', ndim=1),
                           ct.c_int, ]
    __get_Nphi = __libhmpdf._get_Nphi
    __get_Nphi.restype = ct.c_int
    __get_Nphi.argtypes = [ct.c_void_p, ct.POINTER(ct.c_int), ]
    __get_cov_diagnostics = __libhmpdf.hmpdf_get_cov_diagnostics1
    __get_cov_diagnostics.restype = ct.c_int
    __get_cov_diagnostics.argtypes = [ct.c_void_p,
                                      ndpointer(ct.c_double, flags='C_CONTIGUOUS', ndim=1),
                                      ndpointer(ct.c_double, flags='C_CONTIGUOUS', ndim=1),
                                      ndpointer(ct.c_double, flags='C_CONTIGUOUS', ndim=1), ]
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
        if self.__d is not None :
            err = HMPDF.__delete(self.__d)
            if err :
                raise RuntimeError('delete failed.')
            self.__d = None
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
        self.__del__()
    #}}}

    def init(self,
             class_ini: str,
             stype: str,
             zsource: Optional[float]=None,
             **kwargs) -> None :
        """! Initializes the object [calls hmpdf_init()].
        
        @param class_ini     CLASS .ini file
        @param stype         signal type (either "kappa" or "tsz")
        @param zsource       source redshift. Use only if stype="kappa"
        @param **kwargs      optional settings, see the documentation of #hmpdf_configs_e.
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
            arglist.insert(0, ct.c_double(zsource))
        err = HMPDF.__init(self.__d, _C()(class_ini),
                           _E(_E.stypes)(stype), *arglist)
        return self.__ret(err, 'init()')
    #}}}

    def get_op(self,
               binedges: Union[Sequence[float], np.ndarray],
               incl_2h: bool=True,
               noisy: bool=False) -> np.ndarray :
        """! Get the one-point PDF [calls hmpdf_get_op()]
        
        @param binedges      1d, defines how the PDF is binned
        @param incl_2h       whether to include the two-halo term
        @param noisy         whether to include pixel-wise Gaussian noise
        @return the binned PDF (1d)
        """
    #{{{
        out = np.empty(len(binedges)-1)
        err = HMPDF.__get_op(self.__d, len(binedges)-1, binedges, out,
                             incl_2h, noisy)
        return self.__ret(err, 'get_op()', out)
    #}}}

    def get_tp(self,
               phi: float,
               binedges: Union[Sequence[float], np.ndarray],
               noisy: bool=False) -> np.ndarray :
        """! Get the two-point PDF [calls hmpdf_get_tp()]
        
        @param phi           angular separation of the two sky locations (in arcmin)
        @param binedges      1d, defines how the PDF is binned
        @param noisy         whether to include pixel-wise Gaussian noise
        @return the binned PDF (2d)
        """
    #{{{
        out = np.empty((len(binedges)-1)*(len(binedges)-1))
        err = HMPDF.__get_tp(self.__d, phi, len(binedges)-1, binedges, out, noisy)
        out = out.reshape((len(binedges)-1, len(binedges)-1))
        return self.__ret(err, 'get_tp()', out)
    #}}}

    def get_cov(self,
                binedges: Union[Sequence[float], np.ndarray],
                noisy: bool=False) -> np.ndarray :
        """! Get the covariance matrix of the one-point PDF [calls hmpdf_get_cov()]
        
        @param binedges      1d, defines how the covariance matrix is binned
        @param noisy         whether to include pixel-wise Gaussian noise
        @return the binned covariance matrix (2d)
        """
    #{{{
        Nbins = len(binedges) - 1
        out = np.empty(Nbins*Nbins)
        err = HMPDF.__get_cov(self.__d, Nbins, binedges, out, noisy)
        out = out.reshape((Nbins, Nbins))
        return self.__ret(err, 'get_cov()', out)
    #}}}

    def get_Cell(self,
                 ell: Union[Sequence[float], np.ndarray],
                 mode: str='total') -> np.ndarray :
        """! Get the angular power spectrum [calls hmpdf_get_Cell()]
        
        @param ell           1d, the angular wavenumbers
        @param mode          one of "onehalo", "twohalo", "total"
        @return the power spectrum at ell (1d)
        """
    #{{{
        out = np.empty(len(ell))
        err = HMPDF.__get_Cell(self.__d, len(ell), ell, out,
                               _E(_E.corr_types)(mode))
        return self.__ret(err, 'get_Cell()', out)
    #}}}

    def get_Cphi(self,
                 phi: Union[Sequence[float], np.ndarray],
                 mode: str='total') -> np.ndarray :
        """! Get the angular correlation function [calls hmpdf_get_Cphi()]
        
        @param phi           1d, the angular separations (in arcmin)
        @param mode          one of "onehalo", "twohalo", "total"
        @return the correlation function at phi (1d)
        """
    #{{{
        out = np.empty(len(phi))
        err = HMPDF.__get_Cphi(self.__d, len(phi), phi, out,
                               _E(_E.corr_types)(mode))
        return self.__ret(err, 'get_Cphi()', out)
    #}}}

    def get_cov_diagnostics(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """! Get the covariance diagnostics [calls hmpdf_get_cov_diagnostics()]
        
        @return (phi, phiweights, corr_diagn)
        """
    #{{{
        Nphi = ct.c_int(0)
        err = HMPDF.__get_Nphi(self.__d, ct.byref(Nphi))
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
