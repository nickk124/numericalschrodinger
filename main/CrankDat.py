import numpy as np
import scipy as sp
import schrodingerutils as sch

# Crank-Nicholson integrator.
# Returns the full solution array (including
# initial conditions at t=0). Array should be
# of shape (J,N), with J the spatial and N
# the temporal support points.
# Uses tridiag to solve the tridiagonal matrix.
def cranknicholson(x,t,delt,delx,fBNC,V_0):
    J        = x.size
    N        = t.size
    m = 1.0 # I let m be 1 for now, but we can change it later if need be
    q = (sp.constants.hbar)*delt/(4*m*delx**2)
    r = delt/(2*sp.constants.hbar)
    phi = np.zeros((J+2,N),complex)
    # phi_0 = ? Need to figure out what initial phi array will look like
    # check me on these definitions for the tridiag array part
    # I think this is how we should solve it with the CN fn given, but tridiag still confuses me
    a = np.zeros(J, complex) + (-1.j*q)
    b = np.zeros(J, complex)
    c = a.copy()
    for k in range(J):
        b[k] += (1 + 2*1.j*q + 1.j*r*V[k])
        if k == 0 or k == J-1:
            b[k] += 2*1.j*q # account for boundary conditions in middle diagonal end terms
    V = V_0.copy()
    r = np.zeros(len(V),complex) # matrix to store each new RHS term in matrix equation
    # this comes from the mixture of explicit and implicit methods (we need more terms to calculate RHS array vals)
    for j in range(J): # fill y with initial temperature array, leaving space for boundary conditions
        phi[j+1][0] = phi_0[j]
    for n in range(N-1):
        phi[0][n] = fBNC(0,y[:,n]) # update left bound
        phi[J+1][n] = fBNC(1,y[:,n]) # update right bound
        for l in range(J): # fill in RHS values for use in tridiag for current iteration
            r[l] = (1.j*q)*(phi[l][n] + phi[l+2][n]) + (1. - (2.*1.j*q) - (1.j*r*V[l]))*phi[l+1][n]
        phi_next = tridiag(a,b,c,r) # use tridiag to solve now-implicit equation
        for j in range(1,J+1,1): # fill y with CN-solved values
            phi[j][n+1] = phi_next[j-1]
    # to here ??????
    return phi[1:J+1,:]

# Solver for a tridiagonal matrix.
# a,b,c are the lower, center, and upper diagonals,
# r is the RHS vector.
def tridiag(a,b,c,r):
    n    = b.size
    gam  = np.zeros(n,complex)
    u    = np.zeros(n,complex)
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
