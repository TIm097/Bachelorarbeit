import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import complex_ode
from tqdm import tqdm

# (0) Aufstellen und Lösen der DGL (Schrödinger):
# Zustände des 36 D Hilbertraums:
Z = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Zust.txt', unpack = 'True')

# Hamiltonian:
H_j = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Ham_j.txt', unpack = 'True')
H_d = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Ham_d.txt', unpack = 'True')

# Grundzustand:
v0_oct = np.genfromtxt('Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_Ham_psi0.txt', unpack = 'True')

# J auf 1 normiert
U = 0.02
H_0 = - H_j + U* H_d # Hubbard-Hamiltonian ohne zeitabhängiges Potential

E,v = np.linalg.eig(H_0) # Eigs
i = np.argmin(E)
v0 = v[:,i]
E = np.sort(E)

# Grundzustand bei t=0:
psi_0 = v0
Egz_0 = E[0]
E_2 = E[2]
deltaE1 = E_2 - Egz_0
print('deltaE1:',deltaE1)

# Variablen:
E_0_n = 0.3 # E-Feld-Amplitude in Einheiten von a in eV
w = 0.01 # Frequenz

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
    global H_0
    # Berechnung von H_E:
    H_E = np.zeros((36,36))
    for i in range(0,36): # Alle möglichen Zustände
        for k in range(0,4): # Summe über alle Elektronen
            H_E[i,i] += gp(Z[i,k])[0]*np.cos(w*t) + gp(Z[i,k])[1]*np.sin(w*t) # x*cos(wt) + y*sin(wt)
            # H_E[i,i] += -gp(Z[i,k])[0]*np.sin(w*t) + gp(Z[i,k])[1]*np.cos(w*t) # -x*sin(wt) + y*cos(wt) Phasenverschiebung um pi/2
    H = H_0 - E_0_n * H_E # 36 D zeitabhängige Hamiltonmatrix, -H/hbar
    return -1j *H.dot(x) # stimmt

# Startwerte:
N = 2001 # Anzahl Werte
x = np.zeros((36,N))*1j # für RK5
y0 = np.array(psi_0)
x[:,0] = y0 # Imaginärteil ist 0
# x ist complexes 36xN array

t0 = 0
t1 = 30*np.pi
t = np.linspace(t0,t1,N)
np.savetxt('Hubb_Ham_Zeit_txt/Hubb_Ham_Zeit_Linspace.txt', t)

#Ausgabe:
print('E_0*a:', E_0_n, '// w:', w,
'// Intervall für t: [',t0,';',t1,']'
)

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

##Integration mit Runge-Kutta 8:
#z = np.zeros((36,N))*1j # für RK5
#y0 = np.array(psi_0)
#z[:,0] = y0 # Imaginärteil ist 0
## z ist complexes 36xN array
#
#r = complex_ode(f).set_integrator('dop853')
#r.set_initial_value(y0, t0)
#
#i = 1
#print('Runge-Kutta(5) für',N-1,'Werte:')
#with tqdm(total=N-1) as pbar:
#    while r.successful() and r.t < t1:
#        r.integrate(t[i])
#        pbar.update(1)
#        z[:,i] = r.y
#        i += 1
#
#u = z-x
#np.savetxt('Hubb_Ham_Zeit_txt/Runge_Kutta_vgl.txt', np.column_stack([u.real, u.imag]))

np.savetxt('Hubb_Ham_Zeit_txt/Hubb_Ham_Zeit_Lsg.txt', np.column_stack([x.real, x.imag]))
