import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import complex_ode
from tqdm import tqdm

# Zustände des 4D-Hilbertraums:
Z = np.array((1,2,3,4))

# Hamiltonian:
H_0 = np.genfromtxt('Ham1.txt', unpack = 'True')

# Grundzustand:
v = -0.5 * np.ones(4)

psi0 = v

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


def f(t,x):
    E_0_n = 1
    w = 0.5
    # Berechnung von H_E:
    H_E = np.zeros((4,4))
    for i in range(0,4): # Alle möglichen Zustände
        H_E[i,i] += gp(Z[i])[0]*np.cos(w*t) + gp(Z[i])[1]*np.sin(w*t) # x*cos(wt) + y*sin(wt)
    H = E_0_n * H_E - H_0 # 36 D zeitabhängige Hamiltonmatrix, -H/hbar
    return 1j*H.dot(x) # stimmt

# Startwerte:
N = 1001 # Anzahl Werte
x = np.zeros((4,N))*1j # für RK5
y0 = np.array(psi0)
x[:,0] = y0 # Imaginärteil ist 0
# x ist complexes 4xN array

t0 = 0
t1 = 20
t = np.linspace(t0,t1,N)

#Integration mit Runge-Kutta 5:
r = complex_ode(f).set_integrator('dopri5')
r.set_initial_value(y0, t0)

i = 1
print('Runge-Kutta(5) für',N-1,'Werte:')
with tqdm(total=N-1) as pbar:
    while r.successful() and r.t < t1:
        r.integrate(t[i])
        pbar.update(1)
        x[:,i] = r.y
        i += 1

J = np.zeros((4,4))
for i in range(0,4): # Alle Zeilen (jede i-te Zeile bestimmt alpha des i-ten End-Zustandes)
    for j in range(0,4): # Alle Spalten (Multiplikation von Jij mit alpha|j>)
        M = np.zeros((2,1)) # jede erste Komponente negativer Strom, jede zweite Komponente positiver Strom
        M[0,0] = Z[j]
        M[1,0] = Z[j]
        # Sprung gegen Stromrichtung:
        if M[0,0] == 1: # Randbedingung
            M[0,0] = 4
        else:
            M[0,0] -= 1
        # Sprung mit Stromrichtung:
        if M[1,0] == 4: # Randbedingung
            M[1,0] = 1
        else:
            M[1,0] += 1
        for p in range(0,2): # Summe über alle erzeugten Zustände (jeweils mit und gegen Stromrichtung)
            vzp = p*2-1 # -1 oder +1 abhängig von der Stromrichtung
            if M[p,0] == Z[i]: # Vergleich
                J[j,i] = vzp * 1

J = 1j*J
print(J)
Jpsi = np.zeros((N,4))*1j
for n in range(0,N): # alle Zeiten
    Jpsi[n,:] = J.dot(x[:,n])

# Berechnung von psi Jpsi:
psiJpsi = np.zeros((N,1))*1j

for n in range(0,N): # alle Zeiten
    psiJpsi[n] = (x[:,n].real).dot(Jpsi[n,:].real) + (x[:,n].imag).dot(Jpsi[n,:].imag) + 1j* ((x[:,n].real).dot(Jpsi[n,:].imag)-(x[:,n].imag).dot(Jpsi[n,:].real))

plt.plot(t,psiJpsi.real , 'b-', label=r'$Re \left[ \bra{\psi_0(t)} \hat{J} \ket{\psi_0(t)} \right]$')

plt.legend(loc='best')
plt.xlim(t0,t1)
plt.xlabel(r'$t/J^{-1}$')
plt.ylabel(r'$I/J \cdot x$')
plt.savefig('build/Ham_Strom.pdf')
