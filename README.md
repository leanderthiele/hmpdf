Using the python wrapper, a calculation of the weak lensing convergence one-point PDF
for sources at redshift 1 in a cosmology given by a CLASS input file named "example.ini"
would go like this:
```python
binedges = np.linspace(0.0, 0.3, num=101)
with HMPDF() as d :
	d.init('example.ini', 'kappa', 1.0)
	op = d.get_op(binedges)
# do something with the one-point PDF ...
```

For further information, read the
[Documentation](https://leanderthiele.github.io/hmpdf/html/).
