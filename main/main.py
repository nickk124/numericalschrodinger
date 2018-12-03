import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import cranknicholson as cn
import chebfft as cf
import schrodingerutils as ut

hbar = ut.hbar
m = ut.m
dx = 1000
dt = 1000

def boundary(potential,u,u_last):
    J = len(u)-2
    deltak = 1e-21
    c = (np.sqrt((ut.k+deltak)/deltak) - np.sqrt(ut.k))/deltak
    # For now, all of the boundaries are the same and simple
    if potential == 'free' or potential == 'harmonic':
        # Initially we wanted to use open boundary conditions but these require we
        # know the speed of the wave, which would have been difficult to implement
        # Instead, these are periodic
        u[0] = u[-2]
        u[-1] = u[1]
    elif potential == 'infwell' or potential == 'barrier':
        # Reflective conditions: ghost cell values are simply those of the nearest real cell
        # such that du/dx = 0 at ends
        # http://hplgit.github.io/INF5620/doc/pub/sphinx-wave/._main_wave003.html
        u[0] = u[1]
        u[-1] = u[-2]
    elif potential == 'hydrogen':
        u[0] = u[-2]
        u[-1] = u[1]
    return u

# name      type                explanation
# potential string              name of potential energy function
# psi_0     array               initial wavefunction (1D array of size J)
# solver    function pointer    type of scheme (chebfft or cranknicholson)
# J         int                 number of spatial points
# dt        float               time step
# fBNC      function pointer    put in array of length J+2, and apply boundary conditions to it

def schrodinger_solve(potential,solver,psi0_name,J,N,xbounds,dt,fBNC, **kwargs):

    n = kwargs.pop('n', 3)

    n = int(n)

    x = np.linspace(xbounds[0], xbounds[-1], J) #array of x coordinates
    t = np.arange(0, N, dt)
    global dx
    dx = (xbounds[-1]-xbounds[0])/J

    psi_0 = np.zeros(J+2) # Initial guess for
    psi_0[1:-1] = ut.fINC(potential,psi0_name, x, n)
    psi_0 = fBNC(potential, psi_0, psi_0)

    if solver == 'CN':
        psi = cn.cranknicholson(x,t,potential,dt,dx,fBNC,psi_0)
    elif solver == 'CFFT':
        psi = cf.chebyshev_fft(x,t,potential,fBNC,psi_0,sumcount=5)
    return psi # returned psi is a J by N array of the wavefunction

def main():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument("J",type=int,
                        help="number of spatial support points (including boundaries)")
    parser.add_argument("dt",type=float,
                        help="timestep")
    parser.add_argument("solver",type=str,
                        help="schrodinger equation solver:\n"
                             "    CN      : Crank-Nicholson\n"
                             "    CFFT    : Chebyshev-FFT")
    parser.add_argument("potential",type=str,
                        help="potentials:\n"
                             "    free    : constant potential\n"
                             "    infwell : infinite square well\n"
                             "    finwell : finite square well\n"
                             "    barrier : well with barrier at center\n"
                             "    harmonic : harmonic oscillator\n"
                             "    hydrogen : radial component of hydrgoen atom")

    parser.add_argument("initial",type=str,
                    help="initial wavefunction:\n"
                            "    stationarystate   : stationarystate\n"
                            "    boundstate    : bound state\n"
                            "    wavepacket    : Chebyshev-FFT")
    parser.add_argument("-j","--n",type=int,default='3',help="energy level of stationary state (for infinite well)")

    # -----------------------------------------------------
    global dt
    args         = parser.parse_args()
    J            = args.J
    dt           = args.dt
    solver       = args.solver
    potential    = args.potential
    psi0_name    = args.initial
    n            = args.n

    N = 1e3 # Use 1000 time support points
    if potential == 'free': # This one gets boring after the wave dissipates
        N = 1e2

    xbounds = [0,1] # Say we're looking only at the interval [0,1]
    if potential == 'hydrogen':
        xbounds = [0, 10]
    elif potential == 'free':
        xbounds = [-10, 10]
    elif potential == 'harmonic':
        xbounds = [-1, 1]
    Jorig = J # The number of support points in the interval [0,1]
    J *= (xbounds[-1] - xbounds[0]) # The total number of support points (over entire range)

    psi = schrodinger_solve(potential,solver,psi0_name,J,N,xbounds,dt,boundary, n=n)
    if (potential == 'hydrogen'):
        psi = psi[0:Jorig,:]
        xbounds = [0,1]
    elif potential == 'free' or potential == 'harmonic':
        psi = psi[J//2 - Jorig//2 : J//2 + Jorig//2, :]
        if potential == 'free':
            xbounds = [-5,5]

    x = np.linspace(xbounds[0], xbounds[-1], Jorig) # x values in interval [0,1]
    t = np.arange(0, N, dt)
    V = ut.initPotential(potential, x)

    ut._3DPlot(psi,x,t,xbounds,V)
    ut.animPlot(psi,x,t,xbounds,V)

# --------------------------------------------------
main()
