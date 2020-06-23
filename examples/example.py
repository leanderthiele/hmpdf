import sys
import numpy as np

sys.path.append('../')
from hmpdf import HMPDF

## [example_kappa_onepoint]
def example_kappa_onepoint() :
    binedges = np.linspace(0.0, 0.3, num=101)
    with HMPDF() as d :
        d.init('example.ini', 'kappa', 1.0)
        op = d.get_op(binedges)
    # do something with the one-point PDF ...
## [example_kappa_onepoint]

## [example_ell_filter]
def example_ell_filter(ell, a) :
    return float(ell < a[0])
## [example_ell_filter]

## [example_ell_filter_use]
def example_ell_filter_use() :
    ell_max = 5000.0
    with HMPDF() as d :
        d.init('example.ini', 'kappa', 1.0,
               custom_ell_filter=example_ell_filter,
               custom_ell_filter_params=np.array([ell_max, ]))
        # get your results ...
## [example_ell_filter_use]

def main() :
    example_ell_filter_use()

if __name__ == '__main__' :
    main()
