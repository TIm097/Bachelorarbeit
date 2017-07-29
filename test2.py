import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import complex_ode
from tqdm import tqdm

# (0) Aufstellen und Lösen der DGL (Schrödinger):
# Zustände des 36 D Hilbertraums:
Z = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Zust.txt', unpack = 'True')

# Hamiltonian:
H_j = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Ham_j.txt', unpack = 'True')
H_d = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Ham_d.txt', unpack = 'True')

# Grundzustand:
v0 = np.genfromtxt('Stationäre_Systeme/Hubb_Eig_Ergebn/Hubb_Ham_psi0.txt', unpack = 'True')

# J auf 1 normiert
U = 4
H_0 = - H_j + U* H_d # Hubbard-Hamiltonian ohne zeitabhängiges Potential

E,v = np.linalg.eig(H_0) # Eigs
i = np.argmin(E)
print(E[i])
print(v[:,i]+v0)
