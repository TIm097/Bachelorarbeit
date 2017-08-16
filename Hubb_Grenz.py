# Grenzfall U -> infty für das Hubbardmodell
# Plot für kleine und sehr große U und Approximationsgerade für große U

import numpy as np
import matplotlib.pyplot as plt

def linregress(x, y):
    assert len(x) == len(y)

    x, y = np.array(x), np.array(y)

    N = len(y)
    Delta = N * np.sum(x**2) - (np.sum(x))**2

    A = (N * np.sum(x * y) - np.sum(x) * np.sum(y)) / Delta
    B = (np.sum(x**2) * np.sum(y) - np.sum(x) * np.sum(x * y)) / Delta

    sigma_y = np.sqrt(np.sum((y - A * x - B)**2) / (N - 2))

    A_error = sigma_y * np.sqrt(N / Delta)
    B_error = sigma_y * np.sqrt(np.sum(x**2) / Delta)

    return A, A_error, B, B_error

U, E = np.genfromtxt('Stationäre_Systeme/Hubb_Grenz/Hubb_Grenz.txt', unpack = 'True')
U_lin, E_lin = np.genfromtxt('Stationäre_Systeme/Hubb_Grenz/Hubb_Grenz_lin.txt', unpack = 'True')

logE = np.log(np.abs(E))
logU = np.log(U)

logE_lin = np.log(np.abs(E_lin))
logU_lin = np.log(U_lin)

sl, sl_err, i, i_err = linregress(logU_lin, logE_lin)
print('Approximationsgerade:',sl, sl_err, i, i_err)

x = np.linspace(0,8)
plt.plot(x, sl*x + i, label = r'Fit', linewidth = 1.20)
plt.plot(logU, logE, label = r'$E_0$', linewidth = 1.20)
plt.xlim(0,6.1)
plt.ylim(-4,3.2)
plt.grid()
plt.xlabel(r'$\ln(U/J)$')
plt.ylabel(r'$\ln(\lvert E \rvert/J)$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('build/Hubb_Grenz_Plot.pdf')
