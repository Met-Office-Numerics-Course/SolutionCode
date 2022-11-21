import numpy as np
import matplotlib.pylab as plt


def smooth_ic(x):
    """ Initial Condition for Shallow water equations """
    nx = len(x)
    zeros = np.zeros(nx)
    ones = np.ones(nx)
    ic = np.where(np.logical_or(x <= 0.25, x >= 0.75), zeros,
                  0.5*(1+np.cos(4*np.pi*(x-0.5)))*ones)
    return ic


def wave_packet(x):
    """ Initial Condition for Shallow water equations """
    nx = len(x)
    zeros = np.zeros(nx)
    ones = np.ones(nx)
    wp = np.where(np.logical_or(x <= 0.25, x >= 0.75), zeros,
                  0.5*(1+np.cos(4*np.pi*(x-0.5)))*ones)
    wp[::2] *= -1
    return wp


# In the code below we solve the linear shallow water equations

phiref = 1.0  # The mean geopotential

# Grid points
nx = 40
dx = 1./float(nx)
twodx = 2.0*dx

# Grid for phi values
xphi = np.arange(nx)*dx

# Grid for u values - these are the same for unstaggered grids
xu = np.arange(nx)*dx

# Set initial condition
phi = wave_packet(xphi)
u = wave_packet(xu)*np.sqrt(phiref)

# Set initial variables
lfirst = True
phim = np.copy(phi)
um = np.copy(u)
phip = np.copy(phi)
up = np.copy(u)

# Set up plots
plt.plot(xphi, phi,'g+')
plt.plot(xu, u, 'b-')
plt.show()

# timestep
dt = .01

# number of steps to take
nstep = 100

# Take FTCS step on the first step
for istep in range(nstep):
    twodt = 2.*dt
    if lfirst:
        twodt = dt
        lfirst = False
    phip = phim - (np.roll(u, -1) - np.roll(u, 1))*(phiref*twodt/twodx)
    up = um - (np.roll(phi, -1) - np.roll(phi, 1))*(twodt/twodx)
    um = np.copy(u)
    u = np.copy(up)
    phim = np.copy(phi)
    phi = np.copy(phip)

        
plt.plot(xphi, phi, 'g+')
plt.plot(xu, u, 'b-')
plt.show()



