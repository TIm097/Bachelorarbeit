# Besetzungszahloperator-Erwartungswert
import numpy as np
import matplotlib.pyplot as plt

u = np.genfromtxt('Hubb_Ham_Zeit_txt/Hubb_Ham_Zeit_Lsg.txt', unpack = 'True').T # Lösung der SG, real und imaginär untereinander
N = np.round(np.shape(u)[1]/2).astype(int) # Anzahl Iterationsschritte
x = u[:,0:N] + 1j*u[:,N:N*2]

t = np.genfromtxt('Hubb_Ham_Zeit_txt/Hubb_Ham_Zeit_Linspace.txt', unpack = 'True')
t0 = t[0]
t1 = t[np.shape(t)[0]-1]

M = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Besetzung.txt', unpack = 'True')  # Besetzungsmatrix

n = np.zeros((4,N)) # Besetzungszahl Für alle Gitterplätze

for w in range(0,N): #Alle Zeiten
    for i in range(0,36): #Alle Zustände
        for k in range(0,4): #Alle Gitterplätze
            n[k,w] += (x[i,w].real**2 + x[i,w].imag**2) * M[k,i]

plt.plot(t,n[0,:],'y-', label = r'$\langle n_1 \rangle$')
plt.plot(t,n[1,:],'b-', label = r'$\langle n_2 \rangle$')
plt.plot(t,n[2,:],'r-', label = r'$\langle n_3 \rangle$')
plt.plot(t,n[3,:],'g-', label = r'$\langle n_4 \rangle$')
#plt.plot(t,n[0,:]+n[1,:]+n[2,:]+n[3,:],'k-', label = r'$\langle n_\text{alle} \rangle$')
plt.xlim(t0,t1)
plt.ylabel(r'$\langle n \rangle$')
plt.xlabel(r'$t/J^{-1}$')
plt.legend(loc = 'best')
#plt.xticks(np.arange(0, 16*np.pi+np.pi, 2*np.pi), (['0',r'$2\pi$',r'$4\pi$',r'$6\pi$',r'$8\pi$',r'$10\pi$',r'$12\pi$',r'$14\pi$',r'$16\pi$',r'$18\pi$']))
plt.tight_layout(pad=0, h_pad=1.18, w_pad=1.08)
plt.savefig('build/Hubb_Ham_Zeit_Besetzung.pdf')
