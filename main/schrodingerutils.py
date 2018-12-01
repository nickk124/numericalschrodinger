# general utilities file for the numerical schrodinger equation project
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import scipy as sp
from scipy import constants

m = 1e-31 # approx mass of an e-
e = 1.60218e-19 # Charge of e-

e_0 = sp.constants.epsilon_0
hbar = sp.constants.hbar # use hbar constant value from scipy
k = 1e-20

omega = np.sqrt(k/m)

# Sets up the initial conditions for each potential configuration
def fINC(name, x):
    J = len(x)
    mu = np.mean(x)
    if name == 'free':
        width = 1/50 # Variance
        a = 1/(width*np.sqrt(2*np.pi))
        f = a*np.exp(-.5*pow(((x-mu)/width), 2))
    elif name == 'infwell':
        a = 1
        f = np.sqrt(2/a)*np.sin(3*np.pi*x/a) # Equation for third harmonic
    elif name == 'finwell':
        L = x[3*J//4] - x[J//4] # Size of the well
        A = 0.05 # Random guesses at parameters
        B = 0.5
        C = 1
        alpha = k*np.tan(k*L/2)
        # Outside the well, decaying exponential
        f = np.zeros(J)
        f[0:J//4] = A*np.exp(-alpha*(mu - x[0:J//4]))
        f[3*J//4] = B*np.exp(alpha*(mu - x[3*J//4]))
        f[J//4:3*J//4] = C*np.cos(k*(mu - x[J//4:3*J//4]))
    elif name == 'barrier': # Gaussian on the left side of the barrier
        width = 1/50
        a = 1/(width*np.sqrt(2*np.pi))
        mu = x[J//4]
        f = a*np.exp(-.5*pow(((x-mu)/width), 2))
    elif name == 'harmonic':
        f = pow((m*omega/(np.pi*hbar)),.25)*np.exp(-m*omega*pow((x-mu),omega)/(2*hbar))
    #plt.plot(x, f)
    #plt.show()
    return f

def initPotential(name, x): #initialzes a vector corresponding to the potential, evaluated at each X. h is the x step size, x0 is the lowest value of x
    x = x.copy() #array of x coordinates
    J = len(x)
    V = np.zeros(J)
    if name == 'free' or name == 'infwell': #free particle
        pass
    elif name == 'barrier':
        # Value of the barrier was chosen to allow a tiny bit of tunneling
        V[J//2-4:J//2+4] = 3e-34 # Place a big barrier in the center of the well
    elif name == 'finwell':
        V[0:J//4] = 1e-34
        V[3*J//4:] = 1e-34
    elif name == 'harmonic':
        V = 0.5*m*omega*(x-np.mean(x))**2
    elif name == 'hydrogen':
        V[0] = 10000
        V[1:] = -pow(e,2)/(4*np.pi*x[1:])
    return V

### Analytical solutions for the potentials that have them and a getter =============
def getAnalytical(potential):
    if potential == 'free':
        return freeParticle
    elif potential == 'infwell':
        return None
    elif potential == 'finwell':
        return None
    elif potential == 'harmonic':
        return None
    else:
        return None

def freeParticle(x, t):
    J = len(x); N = len(t)
    psi = np.zeros((J, N), dtype = np.complex_)
    for j in range(J):
        for n in range(N):
            psi[j, n] = np.exp(1j*(x[j]-t[n]))
    return psi

def infiniteWell(x, t):
    J = len(x); N = len(t)
    a = 1 # Not sure this constant matters
    psi = np.zeros((J, N), dtype = np.complex_)
    for j in range(J):
        for n in range(N):
            np.sqrt(2/a)*np.sin(np.pi*x[j]/a)

def finiteWell(x, t):
    return 0

def harmonic(x, t):
    return 0

#====================================================================================


def animPlot(psi,x,t,V,analytical=None): #plotting function that creates time-animated function of Psi vs. x
    plotnum = 211
    potPlotnum = 212
    if analytical != None:
        plotnum = 311
        potPlotnum = 313

    fig = plt.figure(num=1,figsize=(8,8),dpi=100,facecolor='white')

    psiReal = psi.real
    psiComplex = psi.imag
    numPlot = fig.add_subplot(plotnum)

    def animateNumerical(i): # Takes the interval number as input: used to index values
        psiRealVal = psiReal[:, i]
        psiComplexVal = psiComplex[:, i]
        numPlot.clear()
        numPlot.plot(x, psiRealVal, label='Real component of $\Psi$')
        numPlot.plot(x, psiComplexVal, label='Imaginary component of $\Psi$')
        numPlot.set_title('t = ' + str(round(t[i], 2)))
        numPlot.set_xlabel('x')
        numPlot.set_ylabel('Numerical $\Psi$')
        numPlot.set_ylim([np.min(psiReal),np.max(psiReal)])
        numPlot.set_xlim([x[0], x[-1]])
        numPlot.legend(loc=1)

    numAni = animation.FuncAnimation(fig, animateNumerical, interval=1)

    # If known, show analytical plot
    if analytical != None:
        truePlot = fig.add_subplot(312)
        psiTrue = analytical(x, t)
        realTrue = psiTrue.real
        complexTrue = psiTrue.imag

        def animateTrue(i): # Takes the interval number as input: used to index values
            realTrueVal = realTrue[:, i]
            complexTrueVal = complexTrue[:, i]
            truePlot.clear()
            truePlot.plot(x, realTrueVal,label='Real component of $\Psi$')
            truePlot.plot(x, complexTrueVal,label='Imaginary component of $\Psi$')
            truePlot.set_title('t = ' + str(round(t[i], 2)))
            truePlot.set_xlabel('x')
            truePlot.set_ylabel('Analytical $\Psi$')
            truePlot.set_ylim([np.min(psi),np.max(psi)])
            truePlot.set_xlim([x[0], x[-1]])
            truePlot.legend(loc=1)

        trueAni = animation.FuncAnimation(fig, animateTrue, interval=1)

    pot = fig.add_subplot(potPlotnum)
    pot.plot(x, V)
    pot.set_xlabel('x')
    pot.set_ylabel('V')
    pot.set_xlim([x[0], x[-1]])
    plt.show()

def _3DPlot(psi,x,t,V,analytical=None): #plotting function that will plot, in 3D, Psi vs. x vs. time
    h = x[1]-x[0] # Distance between support points
    X = len(x) # Number of support points

    plotnum = 121
    if analytical != None:
        plotnum = 221

    solutions = plt.figure(num=1,figsize=(8,8),dpi=100,facecolor='white')
    solutions.suptitle('Numerical and Analytical Solutions for ')
    Jmesh, Nmesh = np.meshgrid(t, x)

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
