# Besetzungszahloperator-Erwartungswert
import numpy as np
import matplotlib.pyplot as plt

u = np.genfromtxt('Hubb_Ham_Zeit_Lösungen/U4_E01_besetzung.txt', unpack = 'True').T # Lösung der SG, real und imaginär untereinander
N = np.round(np.shape(u)[1]/2).astype(int) # Anzahl Iterationsschritte
x = u[:,0:N] + 1j*u[:,N:N*2]

t = np.genfromtxt('Hubb_Ham_Zeit_Lösungen/Linspace.txt', unpack = 'True')
t0 = t[0]
t1 = t[np.shape(t)[0]-1]

M = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Besetzung.txt', unpack = 'True')  # Besetzungsmatrix

n = np.zeros((4,N)) # Besetzungszahl Für alle Gitterplätze

for w in range(0,N): #Alle Zeiten
    for i in range(0,36): #Alle Zustände
        for k in range(0,4): #Alle Gitterplätze
            n[k,w] += (x[i,w].real**2 + x[i,w].imag**2) * M[k,i]

plt.plot(t,n[0,:],label = r'$\langle n_1 \rangle$', linewidth = 0.9)
plt.plot(t,n[1,:],label = r'$\langle n_2 \rangle$', linewidth = 0.9)
plt.plot(t,n[2,:],label = r'$\langle n_3 \rangle$', linewidth = 0.9)
plt.plot(t,n[3,:],label = r'$\langle n_4 \rangle$', linewidth = 0.9)
#plt.plot(t,n[0,:]+n[1,:]+n[2,:]+n[3,:],'k-', label = r'$\langle n_\text{alle} \rangle$')
plt.xlim(t0,t1+0.001)
#plt.ylim(0.996,1.008)
plt.ylabel(r'$\langle n \rangle$')
plt.xlabel(r'$t/\frac{\hbar}{J}$')
plt.legend(loc = 'upper right')
#plt.xticks(np.arange(0, 7*np.pi, np.pi), (['0',r'$1\pi$',r'$2\pi$',r'$3\pi$',r'$4\pi$',r'$5\pi$',r'$6\pi$']))
#plt.tight_layout(pad=0, h_pad=1.18, w_pad=1.08)
plt.grid()
plt.savefig('build/Hubb_Ham_Zeit_Besetzung.pdf')
plt.savefig('Plots/U4_E01_besetzung.pdf')
