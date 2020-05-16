from ctypes import *
import numpy as np
from numpy.ctypeslib import ndpointer, as_ctypes

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
        if not arr.shape[0] == self.l :
            raise TypeError('Expected numpy array to have length %d, \
                             but has length %d'%(self.l, arr.shape[1]))
        return c_void_p(np.ascontiguousarray(arr).ctypes.data)
#}}}

# Enumerations
#{{{
_integr = ['legendre', 'chebyshev', 'gegenbauer', 'jacobi', 'laguerre',
          'hermite', 'exponential', 'rational', 'chebyshev2']
_cl_uncl = ['uncl', 'cl']
_stypes = ['kappa', 'tsz']
_mdefs = ['mdef_c', 'mdef_v', 'mdef_m', ]
_corr_types = ['onehalo', 'twohalo', 'total']
#}}}

class _Configs(object) :
#{{{
    configs = ['N_cores', c_int, 'class_pre', c_char_p,
               'N_z', c_int, 'z_min', c_double, 'z_max', c_double, 'z_source', c_double,
               'N_M', c_int, 'M_min', c_double, 'M_max', c_double,
               'N_signal', c_int, 'signal_min', c_double, 'signal_max', c_double,
               'N_theta', c_int, 'rout_scale', c_double, 'rout_rdef', _E(_mdefs),
               'pixel_side', c_double, 'tophat_radius', c_double, 'gaussian_fwhm', c_double,
               'custom_ell_filter', None, 'custom_ell_filter_params', None,
               'custom_k_filter', None, 'custom_k_filter_params', None,
               'N_phi', c_int, 'phi_max', c_double, 'pixelexact_max', c_int,
               'phi_jitter', c_double, 'phi_pwr', c_double, 'regularize_tp', c_int,
               'monotonize', c_int,
               'zintegr_type', _E(_integr), 'zintegr_alpha', c_double, 'zintegr_beta', c_double,
               'Mintegr_type', _E(_integr), 'Mintegr_alpha', c_double, 'Mintegr_beta', c_double,
               'Duffy08_conc_params', _D(9),
               'Tinker10_hmf_params', _D(10),
               'Battaglia12_tsz_params', _D(15), ]
    def __init__(self) :
        assert len(_Configs.configs)/2 == 74-36, len(_Configs.configs)/2
        self.l = []
    def append(self, key, value) :
        if key not in _Configs.configs :
            raise KeyError('Invalid configuration option %s.'%key)
        if _Configs.configs[_Configs.configs.index(key)+1] is None :
            raise NotImplementedError('Configuration option %s not \
                                       supported in the python wrapper'%key)
        pos = _Configs.configs.index(key)
        self.l.append(c_int(pos/2))
        self.l.append(_Configs.configs[pos+1](value))
    def __call__(self) :
        self.l.append(c_int(len(_Configs.configs)/2)) # corresponds to end_configs
        return self.l
#}}}

class HMPDF(object) :
#{{{
    # interaction with the DLL
    #{{{
    __libhmpdf = CDLL("./libhmpdf.so")
    __new_data = __libhmpdf.new_data
    __new_data.restype = POINTER(c_int)
    __delete = __libhmpdf.delete_data
    __delete.restype = None
    __delete.argtypes = [c_void_p, ]
    __init = __libhmpdf.init
    __init.restype = None
    __get_op = __libhmpdf.get_op
    __get_op.restype = None
    __get_op.argtypes = [c_void_p, c_int,
                         ndpointer(c_double, flags="C_CONTIGUOUS", ndim=1),
                         ndpointer(c_double, flags="C_CONTIGUOUS", ndim=1),
                         c_int, ]
    __get_tp = __libhmpdf.get_tp
    __get_tp.restype = None
    __get_tp.argtypes = [c_void_p, c_double, c_int,
                         ndpointer(c_double, flags="C_CONTIGUOUS", ndim=1),
                         ndpointer(c_double, flags="C_CONTIGUOUS", ndim=1), ]
    __get_cov = __libhmpdf.get_cov
    __get_cov.restype = None
    __get_cov.argtypes = [c_void_p, c_int,
                          ndpointer(c_double, flags="C_CONTIGUOUS", ndim=1),
                          ndpointer(c_double, flags="C_CONTIGUOUS", ndim=1),
                          c_char_p, ]
    __get_Cell = __libhmpdf.get_Cell
    __get_Cell.restype = None
    __get_Cell.argtypes = [c_void_p, c_int,
                           ndpointer(c_double, flags="C_CONTIGUOUS", ndim=1),
                           ndpointer(c_double, flags="C_CONTIGUOUS", ndim=1),
                           c_int, ]
    __get_Cphi = __libhmpdf.get_Cphi
    __get_Cphi.restype = None
    __get_Cphi.argtypes = [c_void_p, c_int,
                           ndpointer(c_double, flags="C_CONTIGUOUS", ndim=1),
                           ndpointer(c_double, flags="C_CONTIGUOUS", ndim=1),
                           c_int, ]
    #}}}
    def __init__(self) :
        self.d = HMPDF.__new_data()
    def init(self, class_ini, stype, **kwargs) :
        c = _Configs()
        # add kwargs to configs
        for k,v in kwargs.items() :
            c.append(k,v)
        HMPDF.__init(self.d, c_char_p(class_ini),
                     _E(_stypes)(stype),
                     *c())
    def delete(self) :
        HMPDF.__delete(self.d)
    def __enter__(self) :
        return self
    def __exit__(self, exc_type, exc_value, exc_traceback) :
        self.delete()
    def get_op(self, binedges, mode='cl') :
        out = np.empty(len(binedges)-1)
        HMPDF.__get_op(self.d, len(binedges)-1, binedges, out,
                       _E(_cl_uncl)(mode))
        return out
    def get_tp(self, phi, binedges) :
        out = np.empty((len(binedges)-1)*(len(binedges)-1))
        binedges = np.ascontiguousarray(binedges)
        HMPDF.__get_tp(self.d, phi, len(binedges)-1, binedges, out)
        return out.reshape((len(binedges)-1, len(binedges)-1))
    def get_cov(self, binedges, fname) :
        Nbins = 0 if binedges is None \
                else len(binedges)-1
        out = None if binedges is None \
              else np.empty(Nbins*Nbins)
        HMPDF.__get_cov(self.d, Nbins, binedges, out, fname)
        if Nbins > 0 :
            return out.reshape((Nbins, Nbins))
    def get_Cell(self, ell, mode='total') :
        out = np.empty(len(ell))
        HMPDF.__get_Cell(self.d, len(ell), ell, out,
                         _E(_corr_types)(mode))
        return out
    def get_Cphi(self, phi, mode='total') :
        out = np.empty(len(phi))
        HMPDF.__get_Cphi(self.d, len(phi), phi, out,
                         _E(_corr_types)(mode))
        return out
#}}}
