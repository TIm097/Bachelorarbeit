import numpy as np
import functools
import matplotlib.pyplot as plt
import sympy as sp
from sympy import Symbol
from sympy.matrices import Matrix, zeros
from scipy.integrate import complex_ode

def f(t,x):
    #A = np.array(((0,1),(1,0)))
    A = np.zeros((2,2))*1j
    A[1,0] = 1
    A[0,1] = 1
    print(1j*A)
    return 1j*A.dot(x)

# x1' = ix2
# x2' = ix1
# -> x1'' = -x1

# Startwerte:
N = 100

x = np.zeros((2,N)) *1j
x0 = np.array((1,1))
x[:,0] = x0

t0 = 0
t1 = 2*np.pi
t = np.linspace(t0,t1,N)

r = complex_ode(f).set_integrator('dopri5')
r.set_initial_value(x0, t0)

i = 1
while r.successful() and r.t < t1:
    r.integrate(t[i])
    x[:,i] = r.y
    i += 1


print(x[:,0], x[:,25], x[:,50], x[:,75])
