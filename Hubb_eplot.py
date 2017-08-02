import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

# Eigenwerte-Plot:
U, E_0, E_1, E_2, E_3, E_4, E_5, E_6, E_7 = np.genfromtxt('Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_eplot.txt', unpack = 'True')

plt.plot(U, E_0, 'k-', label = r'$E_0$')
plt.plot(U, E_1, 'b-', label = r'$E_1$')
plt.plot(U, E_2, 'c-', label = r'$E_2$')
plt.plot(U, E_3, 'r-', label = r'$E_3$')
#plt.plot(U, E_5, 'k-', label = r'$E_5$')
#plt.plot(U, E_4, 'y--', label = r'$E_4$')
#plt.plot(U, np.append(E_3[0:245],E_5[245:]), 'y-', label = r'$E_3$')
#plt.plot(U, E_4, 'r-', label = r'$E_4, E_5$')
#plt.plot(U, np.append(E_5[0:245],E_3[245:]), 'c--', label = r'$E_5$')
#plt.plot(U, E_6, 'b--', label = r'$E_6$')
#plt.plot(U, E_7, 'g--', label = r'$E_7$')

plt.xlabel(r'$U/J$')
plt.ylabel(r'$E/J$')
plt.xlim(0.01,20)
plt.ylim(-4.1,0)
plt.legend(loc = 'best')
plt.grid()
plt.savefig('build/Hubb_eplot.pdf')
