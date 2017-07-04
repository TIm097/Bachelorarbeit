import numpy as np
import functools
import matplotlib.pyplot as plt
from scipy.integrate import complex_ode
from tqdm import tqdm

# Zustände des 36 D Hilbertraums:
Z = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Zust.txt', unpack = 'True')

# Hamiltonian:
H_j = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Ham_j.txt', unpack = 'True')
H_d = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Ham_d.txt', unpack = 'True')

H_0 = 0.3* H_j + 3* H_d # Hubbard-Hamiltonian ohne zeitabhängiges Potential

E,v = np.linalg.eig(H_0) # Eigs

# Grundzustand bei t=0:
psi_0 = v[:,0]
Egz_0 = E[0]

# Variablen:
E_0_n = 1 # E-Feld-Amplitude in Einheiten von a in eV
w = 1 # Frequenz

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
    global E_0_n
    global w
    # Berechnung von H_E:
    H_E = np.zeros((36,36))
    for i in range(0,36): # Alle möglichen Zustände
        for k in range(0,4): # Summe über alle Elektronen
            H_E[i,i] += gp(Z[i,k])[0]*np.cos(w*t) + gp(Z[i,k])[1]*np.sin(w*t) # x*cos(wt) + y*sin(wt)

    H = E_0_n * H_E - H_0 # 36 D zeitabhängige Hamiltonmatrix, -H/hbar
    return 1j*H.dot(x) # stimmt

# Startwerte:
N = 1000 # Anzahl Werte
x = np.zeros((36,N))*1j
y0 = np.array(psi_0)
x[:,0] = y0 # Imaginärteil ist 0
# x ist ein complexes 36xN array

t0 = 0
t1 = 2*np.pi
t = np.linspace(t0,t1,N)

#Integration:
r = complex_ode(f).set_integrator('dopri5')
r.set_initial_value(y0, t0)

i = 1
#print('Fortschritt:')
with tqdm(total=N-1) as pbar:
    while r.successful() and r.t < t1:
        r.integrate(t[i])
        pbar.update(1)
        x[:,i] = r.y
        i += 1


# Normierung
#n = 0
#for m in range(0,36):
#    n += x[m,0].real**2 + x[m,0].imag**2
#
#print(n)

plt.plot(t,x[5,:].real**2 + x[5,:].imag**2,'k-')
plt.xlim(t0,t1)
plt.savefig('build/Hubb_Ham_Zeit_gzplot.pdf')
