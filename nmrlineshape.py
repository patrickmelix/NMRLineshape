#Calculate NMR Lineshapes for Linear Molecules
#
#
# Based on the work of Eisbein et al.
# http://dx.doi.org/10.1016/j.micromeso.2015.02.042
#
# You might also want to read the PhD Thesis of E. Eisbein (2015), Section 4.2.6
#
# Please cite at least the paper and this repository when using this.
#

import numpy as np

__all__ = ['NMRLineshape']

class NMRLineshape(object):
    """NMRLineshape class

    Parameters for initialization:
    vectors: list of numpy arrays
        All vectors representing the orientation of the molecules. If you want a time-average
        just pass all vectors as one list. Will be deepcopied into the instance.
    sigma_parallel: float
        Value of chemical shift tensor when parallel to field.
    sigma_perpendicular: float
        Value of chemical shift tensor when perpendicular to field.
    nBins: int
        Number of interavals to divide the xAxis in. The more points, the finer
        the grid between sigma_parallel and sigma_perpendicular.
    nIntersections: int
        Number of intersection of the triangular grid along the edge of the octahedron
        in the Alderman-Powder-Averaging Routine (N in Alderman et al.).
    """

    def __init__(self, vectors, sigma_parallel=0.0, sigma_perpendicular=1.0,
                 nBins=1000, nIntersections=32):
        from copy import deepcopy

        self._cache = {}

        assert isinstance(vectors, list), ("Vectors must be of type list!")
        for item in vectors:
            assert isinstance(item, np.ndarray), ("Vector must be of type numpy.array")
        self.vectors = deepcopy(vectors)

        assert isinstance(sigma_parallel, float), ("sigma_parallel must be of type float!")
        assert isinstance(sigma_perpendicular, float), ("sigma_perpendicular must be of type float!")
        assert isinstance(nBins, int), ("nBins must be of type int!")
        assert isinstance(nIntersections, int), ("nIntersections must be of type int!")
        self.sigma = ( sigma_parallel, sigma_perpendicular )
        self.nBins = nBins
        self.nIntersections = nIntersections
        xAxis = np.arange(min(self.sigma), max(self.sigma), abs(self.sigma[1]-self.sigma[0])/self.nBins)
        self._cache['shifts'] = xAxis


    @property
    def eigenvalues(self):
        """Eigenvalues of the O Matrix

        Only getter, no setter or deleter.
        """
        if not 'eigenvalues' in self._cache:
            self._calc_eigenvalues()

        return self._cache['eigenvalues']


    @property
    def intensities(self):
        """Intensities

        Only getter, no setter or deleter.
        """
        if not 'intensities' in self._cache:
            self.calc()

        return self._cache['intensities']


    @property
    def shifts(self):
        """Chemical Shift Axis

        Only getter, no setter or deleter.
        """
        return self._cache['shifts']


    def calc(self):
        """Calculate the NMR Lineshape"""

        #build sigma matrix
        sigma = np.zeros((3,3))
        np.fill_diagonal(sigma, self.eigenvalues)

        #powder average
        self._cache['intensities'] = _alderman_interpolation(sigma, self._cache['shifts'],
                                                    N=self.nIntersections)
        return

    def _calc_eigenvalues(self):
        """Calculate the Eigenvalues of O"""

        yValues = np.zeros(self.nBins)
        nVec = len(self.vectors)
        tmpMatrix = np.zeros((3,3))

        ### PhD Thesis E. Eisbein equation 4.23
        #iterate over all vectors
        for i in range(nVec):
            v = self.vectors[i]
            tmpMatrix += np.divide(np.outer(v,v), np.dot(v,v))

        #devide by number of vectors
        tmpMatrix /= nVec
        #1-matrix
        tmpMatrix = np.identity(3) - tmpMatrix

        #get eigenvalues
        eigenVals = np.linalg.eigvals(tmpMatrix)
        #sort them
        idx = eigenVals.argsort()
        eigenVals = eigenVals[idx]
        self._cache['eigenvalues'] = eigenVals

        return

def _alderman_interpolation(sigma, bins, N=32):
    """Interpolation following Alderman et al.

    Method described in Alderman, Solum and Grant
    (10.1063/1.450211)[https://doi.org/10.1063/1.450211].
    This code follows their variable naming

    Parameters:
    sigma: numpy array of shape (3,3)
        Sigma Matrix.
    bins: list or numpy array
        Frequency Axis, can be oriented either way, can be unsorted.
        Must contain lower and upper limit, so len(bins) = len(intensities)+1.
    N: integer
        Number of triangles along one octahedron edge,
        results in N**2 triangles per octahedron face.
    """

    #accuracy is not very good, catch this
    small = 1.0e-8

    def _ten_cases(flow, fhigh, fmin, fmid, fmax):
        """Ten cases of Alderman et al. 1a to 1j"""
        case = None
        if (flow <= fmin) and (fmax < fhigh):
            a = 1.0
        elif (fhigh <= fmin):
            a = 0.0
        elif (flow <= fmin < fhigh <= fmid):
            a = (fhigh-fmin)**2/((fmax-fmin)*(fmid-fmin))
        elif (fmin < flow) and (fhigh <= fmid):
            a = (fhigh-flow)*(fhigh+flow-2*fmin)/((fmax-fmin)*(fmid-fmin))
        elif (flow <= fmin) and (fmid < fhigh <= fmax):
            a = (fmid-fmin)/(fmax-fmin)
            a += (fhigh-fmid)*(2*fmax-fhigh-fmid)/((fmax-fmin)*(fmax-fmid))
        elif (fmin < flow <= fmid < fhigh <= fmax):
            a = (fmid-flow)*(fmid+flow-2*fmin)/((fmax-fmin)*(fmid-fmin))
            a += (fhigh-fmid)*(2*fmax-fhigh-fmid)/((fmax-fmin)*(fmax-fmid))
        elif (fmin < flow <= fmid) and (fmax < fhigh):
            a = (fmid-flow)*(fmid+flow-2*fmin)/((fmax-fmin)*(fmid-fmin))
            a += (fmax-fmid)/(fmax-fmin)
        elif (fmid < flow) and (fhigh <= fmax):
            a = (fhigh-flow)*(2*fmax-fhigh-flow)/((fmax-fmin)*(fmax-fmid))
        elif (fmid < flow <= fmax < fhigh):
            a = (fmax-flow)**2/((fmax-fmin)*(fmax-fmid))
        elif (fmax < flow):
            a = 0
        #Error Check
        else:
            m = "None of the ten cases in the 2D interpolation came up!" +\
                 "This seems like an implementation error?!" +\
                  str([flow, fhigh, fmin, fmid, fmax])
            raise RuntimeError(m)

        #Error Check
        if not (1.0+small >= a >= 0-small):
            m = "This points to an implementation error?!" +\
                  str([a, flow, fhigh, fmin, fmid, fmax])
            raise RuntimeError(m)
        return a



    #build matrix of lmn and distances with indices i and j
    lmn = np.zeros((N+1,N+1), dtype=object)
    #R = np.zeros((N+1,N+1)) #not needed, saving the memory
    weights = np.zeros((N+1,N+1))
    #only upper triangle of matrices matters
    for i in range(N+1):
        for j in range(N-i):
            k = N - i - j
            assert k >= 0
            #append
            tmp = np.sqrt( i**2 + j**2 + k**2 )
            lmn[i,j] = np.array([i,j,k])/tmp
            R = tmp/N
            weights[i,j] = 1/R**3

    #calculate (lmn) sigma (lmn)^T
    #equation 4 in Alderman et al.
    delta = np.zeros((N+1,N+1))
    for i in range(N+1):
        for j in range(N-i):
            delta[i,j] = np.dot(np.dot(lmn[i,j].T, sigma), lmn[i,j])

    ######
    #2D Interpolation
    #####

    #build corner information for all triangles
    #see Fig. 3 b) of Alderman et al.
    cornerIndices = []
    for i in range(N+1):
        for j in range(N-i):
            #corners of the triangle are indices of the lmnList
            firstIdx = (i, j)
            secondIdx = (i, j+1)
            thirdIdx = (i+1, j)
            cornerIndices.append((firstIdx, secondIdx, thirdIdx))
            if i-1 >= 0:
                #opposite direction
                thirdIdx = (i-1, j+1)
                cornerIndices.append((firstIdx, secondIdx, thirdIdx))
    #Must have N**2 triangles defined
    assert len(cornerIndices) == N**2

    #calculate intensities
    intensities = np.zeros(len(bins))
    sortedBins = np.sort(bins)
    #iterate over triangles
    for iTriangle in range(N**2):
        tentF = [ delta[idx] for idx in cornerIndices[iTriangle] ]
        #sort frequencies
        tentF.sort()

        #go through all relevant bins
        for iBin, flow in enumerate(sortedBins):
            if iBin+1 < len(bins):
                fhigh = bins[iBin+1]
            #dirty, but should cover most cases
            else:
                fhigh = flow + 10000
            #cycle if bins not relevant, cases j and b
            if (flow > tentF[2]) or (fhigh <= tentF[0]):
                continue

            #now the ten cases
            a = _ten_cases(flow, fhigh, *tentF)
            #weights following eq. 6 of Alderman et al.
            weight = np.array([ weights[idx] for idx in cornerIndices[iTriangle] ]).mean()

            intensities[iBin] += a * weight

    #normalize intensities
    intensities /= (np.sum(intensities) / len(intensities))

    return intensities
