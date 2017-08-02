import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import complex_ode
from tqdm import tqdm

# (0) Aufstellen und Lösen der DGL (Schrödinger):
# Zustände des 36 D Hilbertraums:
Z = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Zust.txt', unpack = 'True')

# Hamiltonian:
H_j = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Ham_j.txt', unpack = 'True')
H_j_lambda = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/lambda_Hubb_Ham_j.txt', unpack = 'True') # Diago-Hüpfen lambda
H_d = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Ham_d.txt', unpack = 'True')

# Grundzustand:
v0_oct = np.genfromtxt('Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_Ham_psi0.txt', unpack = 'True')

# J auf 1 normiert
U = 4

def Eigen(U, Hj, Hd):
    H_0 = - Hj + U* Hd # Hubbard-Hamiltonian ohne zeitabhängiges Potential
    E,v = np.linalg.eig(H_0) # Eigs
    i = np.argmin(E)
    v0 = v[:,i]
    E = np.sort(E)
    return E, v0

E,v0 = Eigen(U,H_j,H_d)
E_lambda,v0_lambda = Eigen(U,H_j_lambda, H_d)

def deltaE(E):
    return E[2] - E[0]

deltaE = deltaE(E)
#deltaE_lambda = deltaE(E_lambda)

print('deltaE1:',deltaE)
#print('deltaE1_lambda:',deltaE_lambda)

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
def f(t,x,E_0_n,w,H_0):
    # Berechnung von H_E:
    H_E = np.zeros((36,36))
    for i in range(0,36): # Alle möglichen Zustände
        for k in range(0,4): # Summe über alle Elektronen
            H_E[i,i] += gp(Z[i,k])[0]*np.cos(w*t) + gp(Z[i,k])[1]*np.sin(w*t) # x*cos(wt) + y*sin(wt)
            # H_E[i,i] += -gp(Z[i,k])[0]*np.sin(w*t) + gp(Z[i,k])[1]*np.cos(w*t) # -x*sin(wt) + y*cos(wt) Phasenverschiebung um pi/2
    H = H_0 - E_0_n * H_E # 36 D zeitabhängige Hamiltonmatrix, -H/hbar
    return -1j *H.dot(x) # stimmt

# Anzahl Werte:
N = 501

# Intervall
t0 = 0
t1 = 50*np.pi
t = np.linspace(t0,t1,N)
np.savetxt('Hubb_Ham_Zeit_txt/Hubb_Ham_Zeit_Linspace.txt', t)

def DGLloesen(E_0_n, w, U, Hj, Hd, v0, t0, t1, N):
    H_0 = - Hj + U* Hd

    # Startwerte:
    x = np.zeros((36,N))*1j # für RK5
    y0 = np.array(v0)

    x[:,0] = y0 # Imaginärteil ist 0
    # x ist complexes 36xN array

    #Integration mit Runge-Kutta 5:
    t = np.linspace(t0,t1,N)
    r = complex_ode(f).set_integrator('dopri5')
    r.set_initial_value(y0, t0).set_f_params(E_0_n,w,H_0)

    i = 1
    print('Runge-Kutta(5) für',N-1,'Werte:')
    with tqdm(total=N-1) as pbar:
        while r.successful() and r.t < t[N-1]:
            r.integrate(t[i])
            pbar.update(1)
            x[:,i] = r.y
            i += 1
    return x

E_0_n = 0.1
w_range = 51
u = np.zeros((36*w_range,N))*1j

for alpha in range(1,w_range):
    w = alpha/50
    print(alpha,'(',w_range-1,')')
    u[(alpha-1)*36:alpha*36] = DGLloesen(E_0_n, w, U, H_j_lambda, H_d, v0_lambda, t0, t1, N)

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

np.savetxt('Hubb_Ham_Zeit_txt/Hubb_Ham_Zeit_Lsg.txt', np.column_stack([u.real, u.imag]))
