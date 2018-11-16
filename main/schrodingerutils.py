# general utilities file for the numerical schrodinger equation project
import numpy as np

def initPotential(name, X, h, x0): #initialzes a vector corresponding to the potential, evaluated at each X. h is the x step size, x0 is the lowest value of x
    x = np.arange(x0, (X+1)*h, h) #array of x coordinates
    if name == "free": #free particle
        V = np.zeros(X)
        return V
    if name == "infwell":
        V = np.zeros(X)
        V[0] = np.inf
        V[-1] = np.inf
        return V


def _3DPlot(): #plotting function that will plot, in 3D, Psi vs. x vs. t
    return
#IDEA: add support for colored display of wavefunction phase

def animPlot(): #plotting function that creates time-animated function of Psi vs. x
    return