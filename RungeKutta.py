import numpy as np
import matplotlib.pyplot as plt

# RungeKuttaVgl-Daten (Differenzen zwischen RK5 und RK8):
u = np.genfromtxt('Hubb_Ham_Zeit_txt/Runge_Kutta_vgl.txt', unpack = 'True').T
N = np.round(np.shape(u)[1]/2).astype(int) # Anzahl Iterationsschritte
x = u[:,0:N] + 1j*u[:,N:N*2]

t = np.genfromtxt('Hubb_Ham_Zeit_txt/Hubb_Ham_Zeit_Linspace.txt', unpack = 'True')
t0 = t[0]
t1 = t[np.shape(t)[0]-1]

y = np.zeros((N,1))*1j
for i in range(0,N):
    for j in range(0,36):
        y[i] += 1/36 * np.sum(x[j,i].real**2 + x[j,i].imag**2)

plt.plot(t, y.real**2 + y.imag**2, 'k-')

plt.savefig('build/RungeKuttaVgl.pdf')
