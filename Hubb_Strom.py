# Stromoperator Erwartungswert im Zeitmittel
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# Zustände des 36 D Hilbertraums:
Z = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Zust.txt', unpack = 'True')

# Ort:
u = np.genfromtxt('Hubb_Ham_Zeit_txt/Hubb_Ham_Zeit_Lsg.txt', unpack = 'True').T # Lösung der SG, real und imaginär untereinander
N = np.round(np.shape(u)[1]/2).astype(int) # Anzahl Iterationsschritte
x = u[:,0:N] + 1j*u[:,N:N*2]

# Zeit:
t = np.genfromtxt('Hubb_Ham_Zeit_txt/Hubb_Ham_Zeit_Linspace.txt', unpack = 'True')
t0 = t[0]
t1 = t[np.shape(t)[0]-1]

# Stromoperatormatrix:
J = np.genfromtxt('Hubb_Ham_Zeit_txt/Hubb_Strom_Matrix.txt', unpack = 'True')
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
plt.savefig('build/Hubb_Strom.pdf')
