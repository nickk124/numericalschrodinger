import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import matplotlib.pyplot as plt

import cranknicholson as cn
import chebfft as cf
import schrodingerutils as ut

# name      type                explanation
# potential string              name of potential energy function
# psi_0     array               initial wavefunction (1D array of size J)
# solver    function pointer    type of scheme (chebfft or cranknicholson)
# J         int                 number of spatial points
# xbounds   tuple of floats     boundary points for x
# dt        float               time step
# fBNC      function pointer    put in array of length J+2, and apply boundary conditions to it

def schrodinger_solve(potential, psi_0, solver, J, xbounds, dt, FBNC):

    
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
                             "    infwell : infinite square well")

    # -----------------------------------------------------
    args         = parser.parse_args()
    J            = args.J
    dt           = args.dt
    solver       = args.solver
    potential    = args.potential

    #psi_0 = 
    #xbounds = 

    psi, x, t = schrodinger_solve(potential, psi_0, solver, J, xbounds, dt, fBNC)

    #ut._3DPlot()
    #ut.animPlot()

# --------------------------------------------------
main()