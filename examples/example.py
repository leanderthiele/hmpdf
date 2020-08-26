import numpy as np

from hmpdf import HMPDF

## [example_kappa_onepoint]
def example_kappa_onepoint() :
    binedges = np.linspace(0.0, 0.3, num=101)
    with HMPDF() as d :
        d.init('example.ini', 'kappa', 1.0)
        op = d.get_op(binedges)
    # do something with the one-point PDF ...
## [example_kappa_onepoint]

## [example_ell_filter_use]
def example_ell_filter_use() :
    example_ell_filter = lambda ell : float(ell < 5000.0)
    with HMPDF() as d :
        d.init('example.ini', 'kappa', 1.0,
               custom_ell_filter=example_ell_filter)
        # get your results ...
## [example_ell_filter_use]

def main() :
    example_kappa_onepoint()

if __name__ == '__main__' :
    main()
