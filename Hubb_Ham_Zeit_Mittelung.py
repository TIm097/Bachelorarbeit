# Periodenzeitmittelung
import numpy as np
import matplotlib.pyplot as plt

N = np.genfromtxt('Hubb_Ham_Zeit_Anzahl.txt', unpack = 'True')[0] # Anzahl Iterationen bei Lösung der SG

u = np.genfromtxt('Hubb_Ham_Zeit_Lsg.txt', unpack = 'True').T # Lösung der SG
x = u[0:36,:] + 1j*u[36:72]

A = 10 # Anzahl Perioden für ein 20 pi intervall
V = (N-1)/A # Werte pro Periode

x_betrag = np.zeros((36,N))
for a in range(0,36):
    for b in range(0,N):
        x_betrag[a,b] = x[a,b].real**2 + x[a,b].imag**2


z = np.zeros((36,A))*1j # gemittelte Werte
for a in range(0,A):
    for i in range(a*50,(a+1)*50):
        for l in range(0,36):
            z[l,a] += x_betrag[l,i]

z = 1/V * z

T = np.linspace(0+np.pi,A*2*np.pi+np.pi,A)
plt.plot(T,z[0,:],'k-', label = r'Zustand $\ket{1212}$')
plt.plot(T,z[7,:],'r-', label = r'Zustand $\ket{1313}$')
plt.xlim(t0,t1)
plt.ylabel(r'$|\alpha|^2$')
plt.xlabel(r'$\omega t$')
plt.legend(loc = 'best')
#plt.xticks(np.arange(0, 4*np.pi+np.pi, np.pi), (['0',r'$\pi$',r'$2\pi$',r'$3\pi$',r'$4\pi$',r'$5\pi$']))
plt.savefig('build/Hubb_Ham_Zeit_Mittelung.pdf')
