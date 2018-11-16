import numpy as np
import numpy.polynomial.chebyshev as cheb
import schrodingerutils


#initializes an array for wavefunction values, where X is the number of x positions on the grid,
# and T is the number of time steps
X = 10000
T = 1000
psi = np.array((T,X))
