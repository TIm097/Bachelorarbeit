import numpy as np
import functools
import matplotlib.pyplot as plt
import sympy as sp
from sympy import Symbol
from sympy.matrices import Matrix, zeros
from scipy.integrate import ode

# Zustände des 36 D Hilbertraums:
Z = np.genfromtxt('Octave/Hubb_Ham/Hubb_Zust.txt', unpack = 'True')

# Hamiltonian:
H_j = np.genfromtxt('Octave/Hubb_Ham/Hubb_Ham_j.txt', unpack = 'True')
H_d = np.genfromtxt('Octave/Hubb_Ham/Hubb_Ham_d.txt', unpack = 'True')

H_0 = 0.3* H_j + 3* H_d # Hubbard-Hamiltonian ohne zeitabhängiges Potential

E,v = np.linalg.eig(H_0) # Eigs

# Grundzustand bei t=0:
psi_0 = v[:,0]
Egz_0 = E[0]

# Variablen:
aE_0 = Symbol('aE_0') # Gitterkonstante mal E-Feld-Amplitude
w = Symbol('w') # Frequenz

w = w.subs(w,1)
aE_0 = aE_0.subs(aE_0,1)

# Funktion für die Gitterpunkte:
def gp(k): # k ist der jeweilige Gitterplatz, Rückgegeben werden die Koordinaten
    if k == 1:
        return [-1,-1]
    if k == 2:
        return [1,-1]
    if k == 3:
        return [1,1]
    if k == 4:
        return [-1,1]

# Funktion für die rechte Seite der DGL:
def f(t,x):
    global H_0
    global aE_0
    global w
    # Berechnung von H_E:
    H_E = np.zeros((36,36))
    for i in range(0,36): # Alle möglichen Zustände
        for k in range(0,4): # Summe über alle Elektronen
            H_E[i,i] += gp(Z[i,k])[0]*sp.cos(w*t) + gp(Z[i,k])[1]*sp.sin(w*t) # x*cos(wt) + y*sin(wt)

    H = H_0 - aE_0 * H_E # 36 D zeitabhängige Hamiltonmatrix, stimmt
    return H.dot(x).T # stimmt

# Startwerte:
N = 5 # Anzahl Werte
x = np.zeros((36,N))*1j
y0 = np.array(psi_0)
x[:,0] = y0 # Imaginärteil ist 0

t0 = 0
t1 = 2*np.pi
t = np.linspace(t0,t1,N)

# Integration:
r = ode(f).set_integrator('dopri5')
r.set_initial_value(y0, t0)

n = 0 # Normierung
for i in range(0,36):
    n += x[i,0]*np.conj(x[i,0])

print(n)

i = 1
print('Fortschritt:')
while r.successful() and r.t < t1:
    r.integrate(t[i])
    print('--',i,'(',N,')','--')
    x[:,i] = r.y
    n = 0 # Normierung
    for p in range(0,36):
        n += x[p,i]*np.conj(x[p,i])
    print(n)
    x[:,i] *= 1/np.sqrt(n)
    print(x[35,i])
    i += 1



#plt.plot(t,x[1,:].real**2,'kx')
#plt.savefig('build/Hubb_Ham_Zeit_gzplot.pdf')
