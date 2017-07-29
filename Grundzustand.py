# Berechnung des zeitunabhängigen Grundzustandes für gegebenes U,J
import numpy as np

# Zustände des 36 D Hilbertraums:
Z = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Zust.txt', unpack = 'True')

# Hamiltonian:
H_j = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Ham_j.txt', unpack = 'True')
H_d = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Ham_d.txt', unpack = 'True')

J = 0.3 # in eV
U = 3 # in eV
H_0 = -J* H_j + U* H_d # Hubbard-Hamiltonian ohne zeitabhängiges Potential

E,v = np.linalg.eig(H_0) # Eigs

# Grundzustand bei t=0:
psi_0 = v[:,0]
Egz_0 = E[0]

print('J =',J,'|U =',U,'|Egz_0 = ', Egz_0)
np.savetxt('Psi_0.txt',psi_0) #Grundzustand wird abgespeichert
