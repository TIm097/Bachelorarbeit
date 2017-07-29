import numpy as np
import matplotlib.pyplot as plt

U, dE = np.genfromtxt('Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_ediff.txt', unpack = 'True')

print('dE(U=4):',dE[39],
'dE(U=8):',dE[79])
plt.plot(U,dE,'k-',label = r'$\Delta E$')
plt.legend(loc = 'best')
plt.xlim(0,10)
plt.ylim(0,1.1)
plt.xlabel(r'$U/J^{-1}$')
plt.ylabel(r'$\Delta E/J^{-1}$')
plt.savefig('build/Hubb_ediff_plot.pdf')
