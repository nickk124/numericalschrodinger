import numpy as np
import scipy as sp
from scipy import constants
from scipy import special
import schrodingerutils as ut

hbar = ut.hbar #sp.constants.hbar
m = ut.m
e = ut.e

class CFFT:
    def __init__(self, V, dt, dx):
        self.V = V
        self.dt = dt
        self.dx = dx

        Vmin = np.amin(self.V)
        Vmax = np.amax(self.V)

        self.Emin = Vmin
        Emax = ((hbar**2.0 * np.pi**2.0) / (2.0 * m * dx**2.0) ) + Vmax
        self.deltaE = (Emax - self.Emin)/2.0


def D(n):
        if n == 0:
            return 1
        elif n >= 1:
            return 2
def JBess(n,cfft):
    deltaE = cfft.deltaE
    dt = cfft.dt

    return sp.special.jv(n, ((deltaE*dt)/hbar))

def Phi(x,n,H,psi,cfft): #recursion relation for Phi (see equation 9 in paper)
    # H and psi are most recently computed; H is a function pointer, psi is an array
    deltaE = cfft.deltaE
    Emin = cfft.Emin
    dt = cfft.dt

    if n == 0:
        return psi
    if n == 1:
        return -1.j * (1.0 / deltaE) * H(Phi(x,0,H,psi,cfft),x,cfft) + 1.j*(1.0 + Emin/deltaE) * Phi(x,0,H,psi,cfft)
    elif n > 1:
        return -2.j * (1.0 / deltaE) * H(Phi(x,n-1,H,psi,cfft),x,cfft) + 2.j * (1.0 + Emin/deltaE)*Phi(x,n-1,H,psi,cfft) + Phi(x,n-2,H,psi,cfft)

def DDchi(chi,x): #uses fast fourier transforms to evaluate the second derivative of some function chi(x) (discretized in an array)
    Chi = np.fft.fft(chi)
    ddchi = np.fft.ifft(-4.0 * np.pi**2.0 * x**2.0 * Chi) # see top of page 5
    return ddchi #returns array

def Hamiltonian(chi,x,cfft): #for some discretized function chi (array), returns an array of the hamiltonian acting on chi
    V = cfft.V

    return (-1.0 * hbar**2.0 / (2.0*m) )*DDchi(chi,x) + V*chi
    #numpy array element-wise arithmetic is wonderful

# Stability criteria outlined at the end of the VT paper
def isStable():
    M = 1 # This is introduced as the number of complex roots of the Pade approximant
    # Is this just equivalent to sumcount??? Or is it 1?  Trials suggest that it should scale
    h = scipy.constants.h
    v = hbar*cfft.dt/(2*m*cfft.dc**2) # Equation 4.6



# name      type                explanation
# potential string              name of potential energy function
# psi_0     array               initial wavefunction (1D array of size J)
# solver    function pointer    type of scheme (chebfft or cranknicholson)
# J         int                 number of spatial points
# xbounds   tuple of floats     boundary points for x
# dt        float               time step
# fBNC      function pointer    put in array of length J+2, and apply boundary conditions to it

#schrodinger_solve(potential, psi_0, solver, J, xbounds, dt, FBNC):
def chebyshev_fft(x,t,potential,fBNC, psi_0,**kwargs):
    print("running CFFT solver")

    J = len(x)
    N = len(t)
    h = x[1] - x[0]
    #psi = psi_0
    # need to apply boundary conditions for psi_0!!! dependent on potential
    Psi = np.zeros((J+2,N), dtype=np.complex_)
    Psi[:,0] = psi_0

    sumcount = kwargs.pop('sumcount', 10) # pre-chosen amount of terms to do in summation (equation 8 in the paper)

    dx = np.abs(x[1] - x[0])
    dt = np.abs(t[1] - t[0])

    V = ut.initPotential(potential, x)
    cfft = CFFT(V, dt, dx)

    a = np.zeros(sumcount, dtype=np.complex_) #list of a_n (see equation 8)

    deltaE = cfft.deltaE
    Emin = cfft.Emin
    dt = cfft.dt

    for n in range(sumcount):
        a[n] = np.exp( (-1.j*(deltaE + Emin)*dt) / hbar) * D(n) * JBess(n,cfft)

    for time in range(1,N):
        ''' # Nick's original code.  Didn't want to delete but can be done in 3 lines (below)
        psi_old = Psi[:,time-1] #these psi vectors are of size psi
        # need to apply boundary
        psi_new = fBNC(potential, psi_old)

        psi_new_inner = np.zeros(J, dtype=np.complex_)

        for n in range(sumcount):
            psi_new_inner += a[n]*Phi(x,n,Hamiltonian,psi_old[1:J+1], cfft)

        psi_new[1:J+1] = psi_new_inner
        # now psi_new has boundaries generated from last iteration, and inner points generated for current iteration

        Psi[:,time] = psi_new
        '''

        for n in range(sumcount):
            Psi[1:-1,time] += a[n]*Phi(x,n,Hamiltonian,Psi[1:J+1,time-1], cfft)

        Psi[:,time] = fBNC(potential, Psi[:,time])


    return Psi[1:J+1,:]
