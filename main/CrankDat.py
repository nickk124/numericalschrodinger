import numpy as np
import scipy as sp
import schrodingerutils as ut

# Crank-Nicholson integrator.
# Returns the full solution array (including
# initial conditions at t=0). Array should be
# of shape (J,N), with J the spatial and N
# the temporal support points.
# Uses tridiag to solve the tridiagonal matrix
def cranknicholson(x,t,potential,delt,delx,fBNC,psi_0,m,hbar):
    J        = len(x)
    N        = len(t)
    q = hbar*delt/(4*m*delx**2)
    r = delt/(2*hbar)
    psi = np.zeros((J+2,N), dtype=np.complex_)
    # psi_0 = ? Need to figure out what initial psi array will look like
    # check me on these definitions for the tridiag array part
    # I think this is how we should solve it with the CN fn given, but tridiag still confuses me
    a = np.zeros(J, dtype=np.complex_) + (-1.j*q)
    b = np.zeros(J, dtype=np.complex_)
    c = a.copy()

    psi[:, 0] = psi_0

    V = ut.initPotential(potential, x)

    for k in range(J):
        b[k] += (1 + 2*1.j*q + 1.j*r*V[k])
        if k == 0 or k == J-1:
            b[k] += 1.j*q # account for boundary conditions in middle diagonal end terms
    rhs = np.zeros(J,dtype=np.complex_) # matrix to store each new RHS term in matrix equation
    # this comes from the mixture of explicit and implicit methods (we need more terms to calculate RHS array vals)
    #for j in range(1, len(psi) - 2): # fill y with initial temperature array, leaving space for boundary conditions
    #    psi[j+1,0] = psi_0[j]

    for n in range(1,N):
        for l in range(J): # fill in RHS values for use in tridiag for current iteration
            rhs[l] = (1.j*q)*(psi[l,n] + psi[l+2,n]) + (1. - (2.*1.j*q) - (1.j*V[l]))*psi[l+1,n] # deleted factor of r
        psi[1:-1,n] = tridiag(a,b,c,rhs) # use tridiag to solve now-implicit equation
        psi[:,n] = fBNC(potential, psi[:,n-1])
#        for j in range(1,J+1,1): # fill y with CN-solved values
#            psi[j][n+1] = psi_next[j-1]
    # to here ??????
    return psi[1:-1,:]

# Solver for a tridiagonal matrix.
# a,b,c are the lower, center, and upper diagonals,
# r is the RHS vector.
def tridiag(a,b,c,r):
    n    = b.size
    gam  = np.zeros(n,dtype=np.complex_)
    u    = np.zeros(n,dtype=np.complex_)
    bet  = b[0]
    u[0] = r[0]/bet
    for j in range(1,n):
        gam[j] = c[j-1]/bet
        bet    = b[j]-a[j]*gam[j]
        if (bet == 0.0):
            print('[tridiag]: matrix not invertible.')
            exit()
        u[j]   = (r[j]-a[j]*u[j-1])/bet
    for j in range(n-2,-1,-1):
        u[j] = u[j]-gam[j+1]*u[j+1]
    return u

# Driver for the actual integrators. Sets the initial conditions
# and generates the support point arrays in space and time.
# input: J      : number of spatial support points
#        dt0    : timestep
#        minmaxx: 2-element array containing minimum and maximum of spatial domain
#        minmaxt: 2-element array, same for time domain
#        fINT   : integrator (one of ftcs, implicit, cranknicholson)
#        fBNC   : boundary condition function
#        fINC   : initial condition function
def diffusion_solve(J,minmaxx,dt0,minmaxt,fINT,fBNC,fINC,**kwargs):
    kappa   = 1.0
    for key in kwargs:
        if (key=='kappa'):
            kappa = kwargs[key]
    # time and space discretization
    N  = int((minmaxt[1]-minmaxt[0])/dt0)+1
    dt = (minmaxt[1]-minmaxt[0])/float(N-1) # recalculate, to make exact
    dx = (minmaxx[1]-minmaxx[0])/float(J)
    x  = minmaxx[0]+(np.arange(J)+0.5)*dx
    t  = minmaxt[0]+np.arange(N)*dt
    # alpha factor
    alpha    = kappa*dt/dx**2
    print('[diffusion_solve]: alpha = %13.5e' % (alpha))
    print('[diffusion_solve]: N     = %7i' % (N))
    y        = fINT(x,t,alpha,fBNC,fINC)
    return x,t,y
