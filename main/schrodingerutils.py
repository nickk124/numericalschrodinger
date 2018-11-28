# general utilities file for the numerical schrodinger equation project
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import scipy as sp
from scipy import constants

m = 1.0e-31 # approx mass of an e-
hbar = sp.constants.hbar # use hbar constant value from scipy
omega = 100

# begin with a Gaussian wave-packet
def fINC(x):
    width = 1/50 # Variance
    a = 1/(width*np.sqrt(2*np.pi))
    mu = np.mean(x)
    f = a*np.exp(-.5*pow(((x-mu)/width), 2))
    return f

def initPotential(name, x): #initialzes a vector corresponding to the potential, evaluated at each X. h is the x step size, x0 is the lowest value of x
    x = x.copy() #array of x coordinates
    J = len(x)
    V = np.zeros(J)
    if name == 'free' or name == 'infwell': #free particle
        pass
    elif name == 'barrier':
        V[J//2] = 100 # Place a big ass barrier in the center of the well
    elif name == 'finwell':
        V[0:J//4] = 2
        V[3*J//4:] = 2
    elif name == 'harmonic':
        V = 0.5*m*omega*(x-np.mean(x))**2
    return V

def analytical(name, x, t):
    if name=='free':
        return np.exp(1j*(x-t)) # How to add wavenumber and frequency?
    elif name=='infwell':
        a = X/h # Width of the box
        return np.sqrt(2/a)*np.sin(np.pi*x/a)
    else:
         return 0 # Analytical solution unknown


def _3DPlot(psi,x,t,V,analytical=None): #plotting function that will plot, in 3D, Psi vs. x vs. time
    h = x[1]-x[0] # Distance between support points
    X = len(x) # Number of support points

    plotnum = 121
    if analytical != None:
        plotnum = 221

    solutions = plt.figure(num=1,figsize=(8,8),dpi=100,facecolor='white')
    solutions.suptitle('Numerical and Analytical Solutions for ')
    Jmesh, Nmesh = np.meshgrid(t, x)
    print('len Jmesh = ' + str(len(Jmesh)))
    print('len Nmesh = ' + str(len(Nmesh)))

    approx = solutions.add_subplot(plotnum,projection='3d')
    approx.plot_surface(Jmesh[1], Nmesh[0],psi,cmap='rainbow')
    approx.set_zlabel('$\Psi$')
    approx.set_ylabel('t')
    approx.set_xlabel('x')


    if analytical != None:
        true = fig.add_subplot(122, projection='3d')
        true.plot_surface(Jmesh[1], Nmesh[0],anal,cmap='rainbow')
        true.set_zlabel('$\Psi$')
        true.set_ylabel('t')
        true.set_xlabel('x')


    plt.show()
#IDEA: add support for colored display of wavefunction phase

def animPlot(psi,x,t,V,analytical=None): #plotting function that creates time-animated function of Psi vs. x
    plotnum = 211
    potPlotnum = 212
    if analytical != None:
        plotnum = 221
        potPlotnum = 222

    fig = plt.figure(num=1,figsize=(8,8),dpi=100,facecolor='white')
    numPlot = fig.add_subplot(plotnum)

    def animateNumerical(i): # Takes the interval number as input: used to index values
        psiVal = psi[:, i]
        time = t[i]
        numPlot.clear()
        numPlot.plot(x, psiVal)
        numPlot.set_title('t = ' + str(round(t[i], 2)))
        numPlot.set_xlabel('x')
        numPlot.set_ylabel('Numerical $\Psi$')

    ani = animation.FuncAnimation(fig, animateNumerical, interval=1)

    # If known, show analytical plot
    if analytical != None:
        truePlot = fig.add_subplot(122)

        def animateTrue(i): # Takes the interval number as input: used to index values
            psiVal = analytical[:, i]
            time = t[i]
            truePlot.clear()
            truePlot.plot(x, psiVal)
            truePlot.set_title('t = ' + str(round(t[i], 2)))
            truePlot.set_xlabel('x')
            truePlot.set_ylabel('$\Psi$')

        plt.xlabel('x')
        plt.ylabel('Analytical $\Psi$')
        ani = animation.FuncAnimation(fig, animateTrue, interval=1)

    pot = fig.add_subplot(potPlotnum)
    pot.plot(x, V)
    pot.set_xlabel('x')
    pot.set_ylabel('V')

    plt.tight_layout()
    plt.show()
