import numpy as np
import matplotlib.pyplot as plt

U, dE = np.genfromtxt('Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_ediff.txt', unpack = 'True')

plt.plot(U,dE,color = (1,0.5,0),label = r'$\Delta E_{0 \to 2}$')
plt.legend(loc = 'best')
plt.xlim(0,15)
plt.ylim(0,1.4)
plt.xlabel(r'$U/J$')
plt.ylabel(r'$E/J$')
plt.tight_layout()
plt.grid()
plt.savefig('build/Hubb_ediff_plot.pdf')
