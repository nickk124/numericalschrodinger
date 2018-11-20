# general utilities file for the numerical schrodinger equation project
import numpy as np
import matplotlib.pyplot as plot
import matplotlib.animation as animation

def initPotential(name, J, h, x0): #initialzes a vector corresponding to the potential, evaluated at each X. h is the x step size, x0 is the lowest value of x
    x = np.arange(x0, (J+1)*h, h) #array of x coordinates
    V = np.zeros(J)
    if name == 'free': #free particle
        pass
    elif name == 'infwell':
        maxFloat = np.finfo('d').max
        V[0] = maxFloat
        V[-1] = maxFloat
    elif name == 'barrier':
        V[J//2] = 100 # Place a big ass barrier in the center of the well
    return V


def analytical(name, x, t):
    if name=='free':
        return np.exp(1j*(x-t)) # How to add wavenumber and frequency?
    elif name=='infwell':
        a = X/h # Width of the box
        return np.sqrt(2/a)*np.sin(np.pi*x/a)


def _3DPlot(psi, x, t, analytical=None): #plotting function that will plot, in 3D, Psi vs. x vs. time
    h = x[1]-x[0] # Distance between support points
    X = len(x) # Number of support points

    plotnum = 111
    if analytical != None:
        plotnum = 121

    solutions = plt.figure(num=1,figsize=(8,8),dpi=100,facecolor='white')
    solutions.suptitle('Numerical and Analytical Solutions for ')
    mesh = np.meshgrid(x,t)

    approx = fig.add_subplot(plotnum,projection='3d')
    approx.plot_surface(mesh[0],mesh[1],psi,cmap='rainbow')
    approx.set_zlabel('$\Psi$')
    approx.set_ylabel('t')
    approx.set_xlabel('x')
    approx.tick_params(labelsize=ftsz)

    if analytical != None:
        true = fig.add_subplot(122, projection='3d')
        true.plot_surface(mesh[0],mesh[1],anal,cmap='rainbow')
        true.set_zlabel('$\Psi$')
        true.set_ylabel('t')
        true.set_xlabel('x')
        true.tick_params(labelsize=ftsz)


    plt.show()
#IDEA: add support for colored display of wavefunction phase

def animPlot(psi,x,t,analytical=None): #plotting function that creates time-animated function of Psi vs. x
    plotnum = 111
    if analytical == None:
        plotnum = 121

    fig = plt.figure(num=1,figsize=(8,8),dpi=100,facecolor='white')
    numPlot = fig.add_subplot(plotnum)

    def animateNumerical(i): # Takes the interval number as input: used to index values
        psiVal = psi[i]
        xVal = x[i]
        time = t[i]
        numPlot.clear()
        numPlot.plot(xVal, psiVal)

    numPlot.xlabel('x')
    numPlot.ylabel('Numerical $\Psi$')

    # Optional show analytical plot
    if analytical != None:
        truePlot = fig.add_subplot(plotnum)

        def animateTrue(i): # Takes the interval number as input: used to index values
            psiVal = analytical[i]
            xVal = x[i]
            time = t[i]
            truePlotPlot.clear()
            truePlot.plot(xVal, psiVal)

        plt.xlabel('x')
        plt.ylabel('Analytical $\Psi$')

        ani = animation.FuncAnimation(fig, animateTrue, interval=100)
    plt.show()
