import numpy as np
import matplotlib.pyplot as plt

U, diff, anregung_1 = np.genfromtxt('Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_anr_diff.txt', unpack = 'True')

plt.plot(U,diff,'r-', label=r'Differenz zwischen GZ und AZ1')
# Maximum um U = 6.3

plt.legend(loc='best')
plt.xlim(0,50)
plt.ylim(0,0.4)
plt.xlabel(r'$U/J$')
plt.ylabel(r'$E_\text{diff}/J$')
plt.savefig('build/Hubb_diff.pdf')
