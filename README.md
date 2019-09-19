# NMRLineshape

Implements the calculation of NMR lineshapes for linear molecules (e.g. CO2) as presented in Bon et al. ([10.1016/j.micromeso.2015.02.042](http://dx.doi.org/10.1016/j.micromeso.2015.02.042), code and theory by E. Eisbein, check also his PhD Thesis). You should cite the paper and this repo (DOI via Zenodo available soon) if you're publishing data obtained using this code.

Please do contact me with questions or use the GitHub features to contribute.

## Requirements:

* Python Version >3.6 (any should be fine)
* Numpy (any recent version should be fine)

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
