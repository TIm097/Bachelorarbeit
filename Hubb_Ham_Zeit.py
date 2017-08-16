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
# stimmen um ein Minuszeichen überein

# J auf 1 normiert
U = 4

def Eigen(U, Hj, Hd):
    H_0 = - Hj + U* Hd # Hubbard-Hamiltonian ohne zeitabhängiges Potential
    E,v = np.linalg.eig(H_0) # Eigs
    i = np.argmin(E)
    v0 = v[:,i]
    #v2 = v[:,4]
    #v0 = v[:,1]
    #print(E[1])
    #v1 = v[:,2]
    #print(E[2])
    #v2 = v[:,5]
    #print(E[5])
    #v3 = v[:,6]
    #print(E[6])
    #E = np.sort(E)
    return E, v0 #, v1, v2, v3

#E,v0,v1,v2,v3 = Eigen(U,H_j,H_d)
E_lambda,v0_lambda = Eigen(U,H_j_lambda, H_d)

def deltaE(E):
    return E[2] - E[0]

#dE = deltaE(E)
#dE_lambda = deltaE(E_lambda)

#print('deltaE1:',dE)
#print('deltaE1_lambda:',dE_lambda)

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

# Startwerte:
H_0 = - H_j + U* H_d
A_0 = 0.1
A_0 = 0.5*A_0 # wegen gittervektoren
# omega immer von 0.5 bis 2

# Anzahl Werte:
N = 2501

# Intervall
t0 = 0
t1 = 80
t = np.linspace(t0,t1,N)
np.savetxt('Hubb_Ham_Zeit_Lösungen/Linspace.txt', t)

w = 1
# Funktion für die rechte Seite der DGL:
def f(t,x):
    global H_0
    # Berechnung von H_E:
    H_E = np.zeros((36,36))
    for i in range(0,36): # Alle möglichen Zustände
        for k in range(0,4): # Summe über alle Elektronen
            H_E[i,i] += gp(Z[i,k])[0]*np.cos(w*t) + gp(Z[i,k])[1]*np.sin(w*t) # x*cos(wt) + y*sin(wt)
            # H_E[i,i] += -gp(Z[i,k])[0]*np.sin(w*t) + gp(Z[i,k])[1]*np.cos(w*t) # -x*sin(wt) + y*cos(wt) Phasenverschiebung um pi/2
    H = H_0 + A_0 * H_E # 36 D zeitabhängige Hamiltonmatrix, -H/hbar
    return -1j *H.dot(x) # stimmt


def DGLloesen(A_0, w, v0, t0, t1, N):

    # Startwerte:
    x = np.zeros((36,N))*1j # für RK5
    y0 = np.array(v0)

    x[:,0] = y0 # Imaginärteil ist 0
    # x ist complexes 36xN array

    #Integration mit Runge-Kutta 5:
    t = np.linspace(t0,t1,N)
    r = complex_ode(f).set_integrator('dopri5')
    r.set_initial_value(y0, t0)

    i = 1
    print('Runge-Kutta(5) für',N-1,'Werte:')
    with tqdm(total=N-1) as pbar:
        while r.successful() and r.t < t[N-1]:
            r.integrate(t[i])
            pbar.update(1)
            x[:,i] = r.y
            i += 1
    return x

omega_anzahl = 2
u = np.zeros((36*omega_anzahl,N))*1j
for alpha in range(0,omega_anzahl):
    #w = 1.5*alpha/(omega_anzahl-1) + 0.5
    w = alpha +1
    #w = dE
    print(alpha+1,'(',omega_anzahl,')')
    u[(alpha)*36:(alpha+1)*36,:] = DGLloesen(A_0, w, v0_lambda, t0, t1, N)

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

#np.savetxt('Daten_für_Besetzungszahl/SG_Lsg_U8_E01.txt', np.column_stack([u.real, u.imag]))
np.savetxt('Hubb_Ham_Zeit_Lösungen/U4_E01_lambda.txt', np.column_stack([u.real, u.imag]))
