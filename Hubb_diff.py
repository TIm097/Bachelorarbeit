import numpy as np
import matplotlib.pyplot as plt

U, diff, anregung_1 = np.genfromtxt('Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_anr_diff.txt', unpack = 'True')

plt.plot(U,diff,'r-', label=r'Energieaufwand für einfache Spinanregung')
x = np.linspace(0.1,100,10000)
plt.plot(x,4/x,'k-', label=r'Approximation der Spinanregung t(U)')
plt.plot(x,x-4,'b-', label=r'Approximation der Ladungsanregung U-W')
# Maximum um U = 6.3

plt.legend(loc='best')
plt.xlim(0,30)
plt.ylim(0,3)
plt.xlabel(r'$U/J$')
plt.ylabel(r'$E_\text{diff}/J$')
plt.savefig('build/Hubb_diff.pdf')
