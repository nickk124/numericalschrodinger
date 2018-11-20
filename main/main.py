import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import matplotlib.pyplot as plt

import CrankDat as cn
import chebfft as cf
import schrodingerutils as ut


m = 1.0     #Define mass as a global variable

def fillerBoundary(potential, u):
    J = u.shape[0]-2
    # For now, all of the boundaries are the same and simple
    if potential == 'free':
        u[0:J+2,0    ] = -u[0:J+2,1]
        u[0:J+2,J+1  ] = -u[0:J+2,J]
        u[0    ,0:J+2] = -u[1,0:J+2]
        u[J+1  ,0:J+2] = -u[J,0:J+2]
    elif potential == 'infwell':
        u[0:J+2,0    ] = -u[0:J+2,1]
        u[0:J+2,J+1  ] = -u[0:J+2,J]
        u[0    ,0:J+2] = -u[1,0:J+2]
        u[J+1  ,0:J+2] = -u[J,0:J+2]

    elif potential == 'barrier':
        u[0:J+2,0    ] = -u[0:J+2,1]
        u[0:J+2,J+1  ] = -u[0:J+2,J]
        u[0    ,0:J+2] = -u[1,0:J+2]
        u[J+1  ,0:J+2] = -u[J,0:J+2]
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
    x = np.arange(xbounds[0], (J+1)*dx, dx) #array of x coordinates
    t = np.arange(0, N, dt)
    psi_0 = np.zeros(J) # Initial guess for psi

    if solver == 'CN':
        V_0 = np.zeros(J)
        psi, t = cn.cranknicholson(x,t,potential,dt,dx,fBNC,V_0,psi_0,m)
    elif solver == 'CFFT':
        psi, t = cf.chebyshev_fft(x,t,potential,psi_0,m)
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
                             "    barrier : well with barrier at center")

    # -----------------------------------------------------
    args         = parser.parse_args()
    J            = args.J
    dt           = args.dt
    solver       = args.solver
    potential    = args.potential

    N = 1e3 # Use 1000 time support points
    xbounds = [0,1] # Say we're looking only at the interval [0,1]
    psi, x, t = schrodinger_solve(potential,solver,J,N,xbounds,dt,fillerBoundary)

    #ut._3DPlot()
    #ut.animPlot()

# --------------------------------------------------
main()
