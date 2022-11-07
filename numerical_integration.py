from numpy import pi, sin
import matplotlib.pyplot as plt

def f(x):
    # This function returns sin(x)
    return sin(x)

def trap(f, a, b, n):
    # This function calculates the trapezoidal approximation to the
    # integral of f(x) between points x=a and x=b using n intervals

    # calculate h
    h = (b - a) / n

    # initialise integral with the values from the endpoints
    integral = 0.5 * h * (f(a) + f(b))

    # loop over interior points, adding their contributions
    for i in range(1, n):
        integral += h * f(a + i * h)

    return integral

def simp(f, a, b, n):
    # This function calculates Simpson's approximation to the integral
    # of f(x) between points x=a and x=b using n intervals by applying
    # the trapezoidal approximation twice
    return (4/3) * trap(f, a, b, n) - (1/3) * trap(f, a, b, int(n/2))

# initialise some empty lists for storing values to plot:
nvals = []
trap_err = []
simp_err = []

# loop over values of n:
for n in range(10, 100, 10):
    nvals.append(n)

    # compute trapezoidal approximation and its error:
    ans = trap(f, 0, 5*pi, n)
    trap_err.append(abs(2 - ans))

    # compute Simpson's approximation and its error:
    ans = simp(f, 0, 5*pi, n)
    simp_err.append(abs(2 - ans))

# log log plots:
plt.loglog(nvals, trap_err, 'rx')
plt.loglog(nvals, simp_err, 'bx')

plt.show()
