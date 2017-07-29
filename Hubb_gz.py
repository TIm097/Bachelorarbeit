# Plot von Eigenvektorkomponenten im Hubbard-4-Elektronen-System
import numpy as np
import matplotlib.pyplot as plt

# Eigenvektor-Plot:
U, vgz_diag, vgz_para, vgz_rest = np.genfromtxt('Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_gzv.txt', unpack = 'True')
plt.plot(U, vgz_rest, 'k-', label = r'$v_\text{rest}^2$')
plt.plot(U, vgz_para, 'b-', label = r'$v_{6,16,21,31}^2$')
plt.plot(U, vgz_diag, 'r-', label = r'$v_{11,26}^2$')

plt.xlabel(r'$U/J$')
plt.ylabel(r'$v^2$')
plt.xlim(0.1,10)
plt.legend(loc = 'best')
plt.savefig('build/Hubb_gzv_plot.pdf')

plt.close()

# Eigenwerte-Plot:
U, E_0, E_1, E_2, E_3, E_4, E_5, E_6, E_7 = np.genfromtxt('Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_eplot.txt', unpack = 'True')
plt.plot(U, E_0, 'k-', label = r'$E_0$')
x = np.linspace(0.1,25,10000)
plt.plot(x,-12/x,'r-', label = r'erste Korrektur')

plt.xlabel(r'$U/J$')
plt.ylabel(r'$E/J$')
plt.xlim(0.01,20)
plt.ylim(-8,0)
plt.legend(loc = 'best')
plt.savefig('build/Hubb_gze_plot.pdf')

plt.close()
