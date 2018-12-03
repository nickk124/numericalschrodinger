import numpy as np
import scipy as sp
import schrodingerutils as ut

hbar = ut.hbar #sp.constants.hbar
m = ut.m

# Crank-Nicholson integrator.
# Returns the full solution array (including
# initial conditions at t=0). Array should be
# of shape (J,N), with J the spatial and N
# the temporal support points.
# Uses numpy matrix inversion to solve
def cranknicholson(x,t,potential,delt,delx,fBNC,psi_0):
    V = ut.initPotential(potential, x)
    J        = len(x)
    N        = len(t)
    q = hbar*delt/(4*m*delx**2)
    r = delt/(2*hbar)
    psi = np.zeros((J+2,N), dtype=np.complex_)
    psi[:, 0] = psi_0

    A = np.zeros((J,J),dtype=np.complex_) # Matrix operator
    V = ut.initPotential(potential, x)
    for i in range(len(A)):
        if i == 0:
            A[i][i] += (1 + 3*1.j*q + 1.j*r*V[i])
            A[i][i+1] += (-1.j*q)
        elif i == len(A)-1:
            A[i][i] += (1 + 3*1.j*q + 1.j*r*V[i])
            A[i][i-1] += (-1.j*q)
        else:
            A[i][i] += (1 + 2*1.j*q + 1.j*r*V[i])
            A[i][i+1] += (-1.j*q)
            A[i][i-1] += (-1.j*q)
    Ainv = np.linalg.inv(A)

    rhs = np.zeros(J,dtype=np.complex_) # matrix to store each new RHS term in matrix equation
    # this comes from the mixture of explicit and implicit methods (we need more terms to calculate RHS array vals)
    for n in range(1,N):
        for l in range(J): # fill in RHS values for use in tridiag for current iteration
            rhs[l] = (1.j*q)*(psi[l,n-1] + psi[l+2,n-1]) + (1. - (2.*1.j*q) - (1.j*r*V[l]))*psi[l+1,n-1] # deleted factor of r
        psi[1:-1,n] = np.dot(Ainv,rhs)
        psi[:,n] = fBNC(potential, psi[:,n], psi[:,n-1])

    return psi[1:-1,:]
