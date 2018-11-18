import numpy as np
import scipy as sp
import schrodingerutils as ut


# name      type                explanation
# potential string              name of potential energy function
# psi_0     array               initial wavefunction (1D array of size J)
# solver    function pointer    type of scheme (chebfft or cranknicholson)
# J         int                 number of spatial points
# xbounds   tuple of floats     boundary points for x
# dt        float               time step
# fBNC      function pointer    put in array of length J+2, and apply boundary conditions to it

#schrodinger_solve(potential, psi_0, solver, J, xbounds, dt, FBNC):
def chebyshev_fft(x, t, potential, psi_0, J, m):
    hbar = sp.constants.hbar
    
    sumcount = 10 # pre-chosen amount of terms to do in summation (equation 8 in the paper)
    
    dx = np.abs(x[1] - x[0])
    dt = np.abs(t[1] - t[0])
    V = initPotential(potential, J, h, x[0])

    Vmin = np.amin(V)
    Vmax = np.amax(V)

    Emin = Vmin
    Emax = ((hbar**2.0 * np.pi**2.0) / (2.0 * m * dx**2.0) ) + Vmax
    deltaE = (Emax - Emin)/2.0

    def D(n):
        if n == 0:
            return 1
        if n >= 1:
            return 2
    def J(n):
        return sp.special.jv(n, ((deltaE*dt)/hbar))

    def Phi(x,n,H): #recursion relation for Phi

    a = np.zeros(sumcount, dtype=np.complex_) #list of a_n

    for n in range(sumcount):
        a[n] = np.exp( (-1.j*(deltaE + Emin)*dt) / hbar) * D(n) * J(n)


    while True:




    return