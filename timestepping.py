import matplotlib.pyplot as plt
from cmath import exp, sqrt

def F(y, lamda):
    #  F(y) = -lamda y
    return -lamda * y

def G(y, omega):
    #  G(y) = i omega y
    return 1j*omega*y

def forward_euler(f, yn, dt, *args):
    # Computes the forward Euler approximation to y^{n+1} given the
    # value y^n, the right hand side f, and the timestep dt
    # Any arguments required by f can be passed using *args
    return yn + dt * f(yn, *args)

def leapfrog(f, yn, ynm1, dt, *args):
    # Computes the leapfrog approximation to y^{n+1} given the
    # values y^n and y^{n-1}, the right hand side f, and the timestep dt
    # Any arguments required by f can be passed using *args
    return ynm1 + 2 * dt *f(yn, *args)

def matsuno(f, yn, dt, *args):
    # Computes the Matsuno approximation to y^{n+1} given the
    # value y^n, the right hand side f, and the timestep dt
    # Any arguments required by f can be passed using *args
    ystar = yn + dt * f(yn, *args)
    return yn + dt * f(ystar, *args)

def exact_F(t, lamda):
    # The exact solution of the ODE with right hand side F(y) at time t
    return exp(-lamda*t)

def exact_G(t, omega):
    # The exact solution of the ODE with right hand side G(y) at time t
    return exp(1j*omega*t)


# initial time, t=0
t0 = 0

# timestep
dt = 0.1

# max time to integrate until
tmax = 5

# initial condition y(0)=1
y0 = 1

# some lists to store the time levels and solutions
tvals = [t0]
fe_sol = [y0]
lf_sol = [y0]
ms_sol = [y0]
exact = [y0]

# lambda required for right have side function F(y) = -lamda y
lamda = 1

# initial values going in to timeloop
t = t0
fe_yn = y0
lf_yn = y0
ms_yn = y0

# flag for first timestep - needed for leapfrog since we don't have y^{n-1}
first_step = True

# timestepping loop:
while t < tmax:

    #  increment t and store
    t += dt
    tvals.append(t)

    # forward Euler approximation
    fe_ynp1 = forward_euler(F, fe_yn, dt, lamda)
    fe_sol.append(fe_ynp1)
    fe_yn = fe_ynp1

    # leapfrog approximation
    if first_step:
        lf_ynp1 = forward_euler(F, lf_yn, dt, lamda)
        first_step = False
    else:
        lf_ynp1 = leapfrog(F, lf_yn, lf_ynm1, dt, lamda)
    lf_sol.append(lf_ynp1)
    lf_ynm1 = lf_yn
    lf_yn = lf_ynp1

    # Matsuno approximation
    ms_ynp1 = matsuno(F, ms_yn, dt, lamda)
    ms_sol.append(ms_ynp1)
    ms_yn = ms_ynp1
    exact.append(exact_F(t, lamda))

# plot exact solution and approximations
plt.plot(tvals, fe_sol, 'b')
plt.plot(tvals, lf_sol, 'r')
plt.plot(tvals, ms_sol, 'g')
plt.plot(tvals, exact, 'k')
plt.show()


# some lists to store the time levels and solutions
tvals = [t0]
fe_sol = [y0]
lf_sol = [y0]
ms_sol = [y0]
exact = [y0]

# omega required for right have side function F(y) = -lamda y
omega = 1

# initial values going in to timeloop
t = t0
fe_yn = y0
lf_yn = y0
ms_yn = y0

# flag for first timestep - needed for leapfrog since we don't have y^{n-1}
first_step = True

# timestepping loop:
while t < tmax:

    #  increment t and store
    t += dt
    tvals.append(t)

    # forward Euler approximation
    fe_ynp1 = forward_euler(G, fe_yn, dt, omega)
    fe_sol.append(fe_ynp1.real)
    fe_yn = fe_ynp1

    # leapfrog approximation
    if first_step:
        lf_ynp1 = forward_euler(G, lf_yn, dt, omega)
        first_step = False
    else:
        lf_ynp1 = leapfrog(G, lf_yn, lf_ynm1, dt, omega)
    lf_sol.append(lf_ynp1.real)
    lf_ynm1 = lf_yn
    lf_yn = lf_ynp1

    # Matsuno approximation
    ms_ynp1 = matsuno(G, ms_yn, dt, omega)
    ms_sol.append(ms_ynp1.real)
    ms_yn = ms_ynp1
    exact.append(exact_G(t, omega))

# plot exact solution and approximations
plt.plot(tvals, fe_sol, 'b')
plt.plot(tvals, lf_sol, 'r')
plt.plot(tvals, ms_sol, 'g')
plt.plot(tvals, exact, 'k')
plt.show()
