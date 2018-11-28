import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import CrankDat as cn
import chebfft as cf
import schrodingerutils as ut

hbar = ut.hbar
m = ut.m

def boundary(potential, u):
    J = len(u)-2
    # For now, all of the boundaries are the same and simple
    if potential == 'free':
        # Periodic boundary condition: let the wave travel through the boundary unchanged
        # We should stop the simulation when the wave reaches the boundary
        # This is way more confusing than intended, so leaving for now
        # https://www.asc.tuwien.ac.at/~arnold/pdf/graz/graz.pdf
        u[0] = u[-2]
        u[-1] = u[1]
    elif potential == 'infwell' or potential == 'barrier':
        # Reflective conditions: ghost cell values are simply those of the nearest real cell
        # such that du/dx = 0 at ends
        # http://hplgit.github.io/INF5620/doc/pub/sphinx-wave/._main_wave003.html
        u[0] = u[1]
        u[-1] = u[-2]
    return u

# name      type                explanation
# potential string              name of potential energy function
# psi_0     array               initial wavefunction (1D array of size J)
# solver    function pointer    type of scheme (chebfft or cranknicholson)
# J         int                 number of spatial points
# dt        float               time step
# fBNC      function pointer    put in array of length J+2, and apply boundary conditions to it

def schrodinger_solve(potential,solver,J,N,xbounds,dt,fBNC):
    dx = (xbounds[-1]-xbounds[0])/J
    x = np.arange(xbounds[0], J*dx, dx) #array of x coordinates
    t = np.arange(0, N, dt)

    psi_0 = np.zeros(J+2) # Initial guess for psi
    psi_0[1:-1] = ut.fINC(potential, x)
    psi_0 = fBNC(potential, psi_0)

    if solver == 'CN':
        psi = cn.cranknicholson(x,t,potential,dt,dx,fBNC,psi_0)
    elif solver == 'CFFT':
        if psi_0.size > J:
            psi_0 = psi_0[1:J+1].copy()
        psi = cf.chebyshev_fft(x,t,potential,psi_0,fBNC,sumcount = 10)
    return psi, x, t # returned psi is a J by N array of the wavefunction

def main():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument("J",type=int,
                        help="number of spatial support points (including boundaries)")
    parser.add_argument("dt",type=float,
                        help="timestep")
    parser.add_argument("solver",type=str,
                        help="diffusion equation solver:\n"
                             "    CN      : Crank-Nicholson\n"
                             "    CFFT    : Chebyshev-FFT")
    parser.add_argument("potential",type=str,
                        help="potentials:\n"
                             "    free    : constant potential\n"
                             "    infwell : infinite square well\n"
                             "    finwell : finite square well\n"
                             "    barrier : well with barrier at center\n"
                             "    harmonic : harmonic oscillator")

    # -----------------------------------------------------
    args         = parser.parse_args()
    J            = args.J
    dt           = args.dt
    solver       = args.solver
    potential    = args.potential

    N = 1e3 # Use 1000 time support points
    xbounds = [0,1] # Say we're looking only at the interval [0,1]
    psi, x, t = schrodinger_solve(potential,solver,J,N,xbounds,dt,boundary)

    V = ut.initPotential(potential, x)

    #ut._3DPlot(psi, x, t, V)
    ut.animPlot(psi, x, t, V)

# --------------------------------------------------
main()
