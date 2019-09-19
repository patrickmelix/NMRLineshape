[![DOI](https://zenodo.org/badge/203989389.svg)](https://zenodo.org/badge/latestdoi/203989389)

# NMRLineshape

Implements the calculation of NMR lineshapes for anisotropic, linear molecules (e.g. CO<sub>2</sub>) as presented in Bon et al. ([10.1016/j.micromeso.2015.02.042](http://dx.doi.org/10.1016/j.micromeso.2015.02.042)). The original code and theory by E. Eisbein are published in his [PhD Thesis](http://d-nb.info/1109807481) (only in German and not available online). *Working on my own publication including a description of the methodology at the moment.*

## Requirements:

* Python Version >3.6 (any should be fine)
* Numpy (any recent version should be fine)

## Citing:
You should cite the [original paper](http://dx.doi.org/10.1016/j.micromeso.2015.02.042) and my paper (*working on it, use the Zenodo Link [10.5281/zenodo.3450684](https://doi.org/10.5281/zenodo.3450684) until then*) if you're publishing data obtained using this code.

## Contributing:

Please do contact me with questions or use the GitHub features to contribute.

## Usage:
```python
from nmrlineshape import NMRLineshape
import matplotlib.pyplot as plt

#Frozen, all vectors aligned
frozen = NMRLineshape([np.ones(3) for i in range(10)], nBins=1000, nIntersections=256)

#Gaseous, homogeneous random orientation
gas = NMRLineshape([np.random.normal(size=3) for i in range(1000000)], nBins=1000, nIntersections=256)

print("Eigenvalues of gaseous matrix:")
print(gas.eigenvalues)

print("Eigenvalues of frozen matrix:")
print(frozen.eigenvalues)

#Plot results
ax = plt.figure().gca()
ax.set_xlim([1.0,0.0])
ax.set_ylim([0.0,5.0])
plt.xlabel(r"Relative Chemical Shift $\sigma / \Delta$")
plt.ylabel(r"Intensity [arb. unit]")
plt.plot(frozen.shifts, frozen.intensities,marker=None,color='blue',label='Frozen',linestyle='-')
plt.plot(gas.shifts, gas.intensities,marker=None,color='red',label='Gas',linestyle='-')
plt.legend()
plt.tight_layout()
plt.savefig('lineshape-frozen-gas.png')
plt.show()
```
