# Stromoperator Erwartungswert im Zeitmittel
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# Zustände des 36 D Hilbertraums:
Z = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Zust.txt', unpack = 'True')

# Orte:
u = np.genfromtxt('Hubb_Ham_Zeit_Lösungen/U0_E01_v0.txt', unpack = 'True').T # Lösung der SG, real und imaginär untereinander
N = np.round(np.shape(u)[1]/2).astype(int) # Anzahl Iterationsschritte pro Lösung
W = np.round(np.shape(u)[0]/36).astype(int)
print(W)
x = u[:,0:N] + 1j*u[:,N:N*2]

u1 = np.genfromtxt('Hubb_Ham_Zeit_Lösungen/U0_E01_v1.txt', unpack = 'True').T # Lösung der SG, real und imaginär untereinander
x1 = u1[:,0:N] + 1j*u1[:,N:N*2]

u2 = np.genfromtxt('Hubb_Ham_Zeit_Lösungen/U0_E01_v2.txt', unpack = 'True').T # Lösung der SG, real und imaginär untereinander
x2 = u2[:,0:N] + 1j*u2[:,N:N*2]

u3 = np.genfromtxt('Hubb_Ham_Zeit_Lösungen/U0_E01_v3.txt', unpack = 'True').T # Lösung der SG, real und imaginär untereinander
x3 = u3[:,0:N] + 1j*u3[:,N:N*2]

#u_lambda = np.genfromtxt('Hubb_Ham_Zeit_txt/Hubb_Ham_Zeit_Lsg_lambda.txt', unpack = 'True').T # Lösung der SG mit Diagohüpfen, real und imaginär untereinander
#x_lambda = u_lambda[:,0:N] + 1j*u_lambda[:,N:N*2]

# Zeit:
t = np.genfromtxt('Hubb_Ham_Zeit_Lösungen/Linspace.txt', unpack = 'True')
t0 = t[0]
t1 = t[np.shape(t)[0]-1]

# Stromoperatormatrix:
J = np.genfromtxt('Hubb_Ham_Zeit_Lösungen/Hubb_Strom_Matrix.txt', unpack = 'True').T
J = 1j*J

#Erwartungswert des Stromoperators in Abhängigkeit von der Zeit t:
def StromErwartungswert(z):
    global J
    EJ =  np.zeros((N,1))*1j
    for n in range(0,N): # alle Zeiten
        EJ[n] = np.conj(z[:,n]).dot(J.dot(z[:,n]))

    q = 1
    EJ = q * EJ
    return EJ

#EJ_lambda =  np.zeros((N,1))*1j

#for n in range(0,N): # alle Zeiten
#    EJ_lambda[n] = np.conj(x_lambda[:,n]).dot(J.dot(x_lambda[:,n]))

#EJ_lambda = q * EJ_lambda

SE1 = 0.25*(StromErwartungswert(x[0*36:(0+1)*36,0:N])+StromErwartungswert(x1[0*36:(0+1)*36,0:N])+StromErwartungswert(x2[0*36:(0+1)*36,0:N])+StromErwartungswert(x3[0*36:(0+1)*36,0:N]))
SE2 = 0.25*(StromErwartungswert(x[1*36:(1+1)*36,0:N])+StromErwartungswert(x1[1*36:(1+1)*36,0:N])+StromErwartungswert(x2[1*36:(1+1)*36,0:N])+StromErwartungswert(x3[1*36:(1+1)*36,0:N]))

# Mittelwert:
def StromMittelwert(z):
    J_t = StromErwartungswert(z)
    M = 0.25*np.sum(J_t)/N
    return M

#J_w = np.zeros((W,1)) # Strommittelwerte in Abhängigkeit von W
#w_lin = np.zeros((W,1))
#for w in range(0,W):
#    J_w[w] = StromMittelwert(x[w*36:(w+1)*36,0:N])
#    w_lin[w] = 1.5*w/(W-1) + 0.5
#    w_lin[w] = w +1
#
#print('Strommittelwerte (omega aufsteigend):', J_w)
#M_lambda = np.sum(EJ_lambda)/N
#print('Mittelwert:', M,
#'Mittelwert(lambda):', M_lambda)

plt.plot(t, SE1.real, label= r'$\omega = 1.0 \, J/\hbar$', linewidth = 1)
plt.plot(t, SE2.real, label= r'$\omega = 2.0 \, J/\hbar$', linewidth = 0.75)
#plt.plot((76,76),(-0.03,0.03),'k-.', label= r'$t_\text{max}$', linewidth = 0.85)
#plt.plot(t, StromErwartungswert(x[0*36:(0+1)*36,0:N]).real, label= r'$\omega = 1.0$')
#plt.plot(t, StromErwartungswert(x[1*36:(1+1)*36,0:N]).real, label= r'$\omega = 2.0$')
#plt.plot(t, StromErwartungswert(x[8*36:(8+1)*36,0:N]).real, label= r'$\omega = 1.4$')
#plt.plot(t, StromErwartungswert(x[10*36:(10+1)*36,0:N]).real, label= r'$\omega = 1.6$')
#plt.plot(t, StromErwartungswert(x[12*36:(12+1)*36,0:N]).real, label= r'$\omega = 1.8$')
#plt.plot(t, StromErwartungswert(x[14*36:(14+1)*36,0:N]).real, label= r'$\omega = 2.0$')

#plt.plot(w_lin, J_w, 'rx')
#, label=r'$\bra{\psi_0(t)} \hat{J} \ket{\psi_0(t)}$')
plt.legend(loc='best')
plt.xlim(t0,t1)
plt.ylim(-0.0075, 0.01)
plt.xlabel(r'$t/\tfrac{\hbar}{J}$')
plt.ylabel(r'$I/\tfrac{J \symup{e}}{\hbar}$')
plt.grid()
plt.tight_layout()
plt.savefig('Plots/U0_E01_schoen.pdf')
plt.savefig('build/Hubb_Strom.pdf')

#plt.close()

#plt.plot(w_lin, J_w, 'rx')
#plt.xlabel(r'$\omega/\frac{J}{\hbar}$')
#plt.ylabel(r'$I/J \cdot x$')
#plt.savefig('Plots/U0_E01_omega_schoen.pdf')
