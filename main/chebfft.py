import numpy as np
import scipy as sp
from scipy import constants
import schrodingerutils as ut

hbar = sp.constants.hbar

def D(n):
        if n == 0:
            return 1
        elif n >= 1:
            return 2
def J(n):
    return sp.special.jv(n, ((deltaE*dt)/hbar))

def Phi(x,n,H,psi): #recursion relation for Phi (see equation 9 in paper)
    # H and psi are most recently computed; H is a function pointer, psi is an array
    if n == 0:
        return psi
    if n == 1:
        return -1.j * (1.0 / deltaE) * H( Phi(x,0,H,psi) ) + 1.j*(1.0 + Emin/deltaE) * Phi(x,0,H,psi)
    elif n > 1:
        return -2.j * (1.0 / deltaE) * H( Phi(x,n-1,H,psi) ) + 2.j * (1.0 + Emin/deltaE)*Phi(x,n-1,H,psi) + Phi(x,n-2,H,psi)

def DDchi(chi,x): #uses fast fourier transforms to evaluate the second derivative of some function chi(x) (discretized in an array)
    Chi = np.fft.fft(chi)
    ddchi = np.fft.ift(-4.0 * np.pi**2.0 * x**2.0 * Chi) # see top of page 5
    return ddchi #returns array

def Hamiltonian(chi,V,x,m): #for some discretized function chi (array), returns an array of the hamiltonian acting on chi
    return (-1.0 * hbar**2.0 / (2.0*m) )*DDchi(chi,x) + V*chi
    #numpy array element-wise arithmetic is wonderful


# name      type                explanation
# potential string              name of potential energy function
# psi_0     array               initial wavefunction (1D array of size J)
# solver    function pointer    type of scheme (chebfft or cranknicholson)
# J         int                 number of spatial points
# xbounds   tuple of floats     boundary points for x
# dt        float               time step
# fBNC      function pointer    put in array of length J+2, and apply boundary conditions to it

#schrodinger_solve(potential, psi_0, solver, J, xbounds, dt, FBNC):
def chebyshev_fft(x,t,potential,psi_0,m,fBNC,**kwargs):
    J = len(x)
    N = len(t)
    #psi = psi_0
    # need to apply boundary conditions for psi_0!!! dependent on potential
    Psi = np.zeros((J+2,N), dtype=np.complex_)
    Psi[1:J+1,0] = psi_0

    sumcount = 10 # pre-chosen amount of terms to do in summation (equation 8 in the paper)
    for key in kwargs:
        if (key=='sumcount'):
            sumcount = kwargs[key]

    dx = np.abs(x[1] - x[0])
    dt = np.abs(t[1] - t[0])

    V = ut.initPotential(potential, J, h, x[0])
    Vmin = np.amin(V)
    Vmax = np.amax(V)

    Emin = Vmin
    Emax = ((hbar**2.0 * np.pi**2.0) / (2.0 * m * dx**2.0) ) + Vmax
    deltaE = (Emax - Emin)/2.0

    a = np.zeros(sumcount, dtype=np.complex_) #list of a_n (see equation 8)

    for n in range(sumcount):
        a[n] = np.exp( (-1.j*(deltaE + Emin)*dt) / hbar) * D(n) * J(n)


    for time in range(1,N):
        psi_old = Psi[:,time-1] #these psi vectors are of size psi
        # need to apply boundary conditions
        psi_new = fBNC(potential, psi_old)

        psi_new_inner = np.zeros(J, dtype=np.complex_)

        for n in range(sumcount):
            psi_new_inner += a[n]*Phi(x,n,Hamiltonian,psi_old[1:J+1])

        psi_new[1:J+1] = psi_new_inner 
        # now psi_new has boundaries generated from last iteration, and inner points generated for current iteration

        Psi[:,time] = psi_new

    return Psi[1:J+1,:]
