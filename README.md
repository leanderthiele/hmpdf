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
binedges = np.linspace(0.0, 0.3)
with HMPDF() as d :
    d.init('example.ini', 'kappa', 1.0)
    op = d.get_op(binedges)
# do something with the one-point PDF ...
```

Applications of the code can be found in
* [Thiele, Hill, Smith 2019](https://ui.adsabs.harvard.edu/abs/2019PhRvD..99j3511T")
  (for the tSZ field)
* [Thiele, Hill, Smith 2020](https://ui.adsabs.harvard.edu/abs/2020arXiv200906547T")
  (for the weak lensing convergence field)

which we also request you cite if you use the code.

For further information on building the code and using it, please read the
[Documentation](https://leanderthiele.github.io/hmpdf/html/).
