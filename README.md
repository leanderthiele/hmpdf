The hmpdf code computes one- and two-point PDFs,
the covariance matrix of the one-point PDF,
and the angular power spectrum/correlation function of cosmological fields.

The formalism is based on the halo model,
thus, only fields that are mainly sourced by virialized matter make sense.

Currently, two fields are supported:
* thermal Sunyaev-Zel'dovich effect (Compton-y)
* weak lensing convergence

The code is interfaced both through C header files as well as through a convenient python wrapper.

Using the python wrapper, a calculation of the weak lensing convergence one-point PDF
for sources at redshift 1 in a cosmology given by a CLASS input file named "example.ini"
would go like this:
```python
from hmpdf import HMPDF
binedges = np.linspace(0.0, 0.3, num=101)
with HMPDF() as d :
    d.init('example.ini', 'kappa', 1.0)
    op = d.get_op(binedges)
# do something with the one-point PDF ...
```

For further information on building the code and using it, please read the
[Documentation](https://leanderthiele.github.io/hmpdf/html/).
