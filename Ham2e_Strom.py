import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import complex_ode
from tqdm import tqdm

# Zustände:
Z = np.array([(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)])

# Hamilton-Matrix, Tight-Binding-Anteil:
H_J = np.zeros((6,6))

for i in range(1,6): # Zeile
    for j in range(0,i): # Spalte
        for k in range(0,2): # beide Elektronen
            M = np.zeros((2,2))
            M[0,:] = Z[j,:]
            M[1,:] = Z[j,:]
            # Hochhüpfen:
            if M[0,k] == 4: # Randbedingung
                M[0,k] = 1
            else:
                M[0,k] += 1
            # Runterhüpfen:
            if M[1,k] == 1: # Randbedingung
                M[1,k] = 4
            else:
                M[1,k] -= 1
            for l in range(0,2): # Hoch und Runter
                if M[l,0] > M[l,1]: # Dann Permutieren!
                    M[l,:] = np.array((M[l,1],M[l,0]))
                    if (M[l,:]-Z[i,:]).dot(M[l,:]-Z[i,:]) == 0:
                        H_J[j,i] -=1
                        H_J[i,j] -=1
                else:
                    if (M[l,:]-Z[i,:]).dot(M[l,:]-Z[i,:]) == 0:
                        H_J[j,i] +=1
                        H_J[i,j] +=1

H_J *= -1

# Hamilton-Matrix,Potential(Bandisolator)-Anteil:
a = 1*10**(-20) 
H_a = np.zeros((6,6))
for i in range(0,6):
    if Z[i,0] == 1 or Z[i,0] == 3:
        H_a[i,i] += a
    else:
        H_a[i,i] -=a
    if Z[i,1] == 1 or Z[i,1] == 3:
        H_a[i,i] += a
    else:
        H_a[i,i] -=a

H_0 = -H_J + H_a

Eu,vu = np.linalg.eig(H_0)
Esort = np.argsort(Eu)
E = np.zeros((6,1))
v = np.zeros((6,6))
for i in range(0,6):
    E[i] = Eu[Esort[i]]
    v[i,:] = vu[Esort[i],:]

print('E=',E
,'v=',v)

psi_0 = v[1,:]

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
    H_E = np.zeros((6,6))
    for i in range(0,6): # Alle möglichen Zustände
        for k in range(0,2): # Summe über alle Elektronen
            H_E[i,i] += gp(Z[i,k])[0]*np.cos(w*t) + gp(Z[i,k])[1]*np.sin(w*t) # x*cos(wt) + y*sin(wt)
            # H_E[i,i] += -gp(Z[i,k])[0]*np.sin(w*t) + gp(Z[i,k])[1]*np.cos(w*t) # -x*sin(wt) + y*cos(wt) Phasenverschiebung um pi/2
    H = H_0 - E_0_n * H_E # 36 D zeitabhängige Hamiltonmatrix, -H/hbar
    return -1j *H.dot(x) # stimmt

# Startwerte:
N = 2001 # Anzahl Werte
x = np.zeros((6,N))*1j # für RK5
y0 = np.array(psi_0)
x[:,0] = y0 # Imaginärteil ist 0
# x ist complexes 36xN array

t0 = 0
t1 = 30*np.pi
t = np.linspace(t0,t1,N)

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

# Stromoperator Erwartungswert:
J = np.zeros((6,6))
for i in range(0,6): # Alle Zeilen (jede i-te Zeile bestimmt alpha des i-ten End-Zustandes)
    for j in range(0,6): # Alle Spalten (Multiplikation von J mit |j>)
        M = np.zeros((4,2)) # jede erste Komponente negativer Strom, jede zweite Komponente positiver Strom
        for s in range(0,2): # Alle Elektronen
            M[s*2,:] = Z[j,:]
            M[s*2+1,:] = Z[j,:]
            # Sprung gegen Stromrichtung:
            if M[s*2,s] == 1: # Randbedingung
                M[s*2,s] = 4
            else:
                M[s*2,s] -= 1
            # Sprung mit Stromrichtung:
            if M[s*2+1,s] == 4: # Randbedingung
                M[s*2+1,s] = 1
            else:
                M[s*2+1,s] += 1
            if M[s*2,0] != M[s*2,1] or M[s*2+1,0] != M[s*2+1,1]: # Pauliverbot
                for p in range(s*2,s*2+2): # Summe über alle erzeugten Zustände (jeweils mit und gegen Stromrichtung)
                    vzp = (p-2*s)*2-1 # -1 oder +1 abhängig von der Stromrichtung

                    R = np.zeros((2,2)) # Permutation

                    R[0,:] = M[p,:]
                    R[0,:] = np.flip(R[0,:],0) # Ein Flip

                    R[1,:] = M[p,:] # kein Flip

                    for l in range(0,2): # Vergleich mit allen Permutationen mit Zustand |j>
                        vzl = np.sign(l-0.5)# -1 oder +1 abhängig von der Anzahl an Permutationen
                        if (R[l,:] - Z[i,:]).T.dot(R[l,:] - Z[i,:]) == 0: # Vergleich etwas umständlich, aber funktioniert
                            J[j,i] = vzp * vzl * 1

print(J)
J = 1j*J

#Erwartungswert des Stromoperators in Abhängigkeit von der Zeit t:
EJ =  np.zeros((N,1))*1j
for n in range(0,N): # alle Zeiten
    EJ[n] = np.conj(x[:,n]).dot(J.dot(x[:,n]))

q = 1
EJ = q * EJ

# Mittelwert:
M = np.sum(EJ)/N
print('Mittelwert:', M)

plt.plot(t, EJ.real, 'r-', label=r'$\bra{\psi_0(t)} \hat{J} \ket{\psi_0(t)}$')
plt.legend(loc='best')
plt.xlim(t0,t1)
plt.xlabel(r'$t/J^{-1}$')
plt.ylabel(r'$I/J \cdot x$')
plt.tight_layout(pad=0, h_pad=1.10, w_pad=1.08)
plt.savefig('build/Ham2e_Strom.pdf')
