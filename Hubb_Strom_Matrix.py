import numpy as np

# Zustände des 36 D Hilbertraums:
Z = np.genfromtxt('Stationäre_Systeme/Hubb_Ham/Hubb_Zust.txt', unpack = 'True')

J = np.zeros((36,36))
for i in range(0,36): # Alle Zeilen (jede i-te Zeile bestimmt alpha des i-ten End-Zustandes)
    for j in range(0,36): # Alle Spalten (Multiplikation von J mit |j>)
        M = np.zeros((8,4)) # jede erste Komponente negativer Strom, jede zweite Komponente positiver Strom
        for s in range(0,4): # Alle Elektronen
            M[s*2,:] = Z[j,:]
            M[s*2+1,:] = Z[j,:]
            # Sprung gegen Stromrichtung:
            if M[s*2,s] == 1: # Randbedingung
                M[s*2,s] = 4
            else:
                M[s*2,s] -= 1
            # Sprung mit Stromrichtung:
            if M[s*2+1,s] == 4: # Randbedingung
                M[s*2+1,s] = 1
            else:
                M[s*2+1,s] += 1
            if M[s*2,0] != M[s*2,1] or M[s*2+1,2] != M[s*2,3]: # Pauliverbot
                for p in range(s*2,s*2+2): # Summe über alle erzeugten Zustände (jeweils mit und gegen Stromrichtung)
                    vzp = (p-2*s)*2-1 # -1 oder +1 abhängig von der Stromrichtung

                    R = np.zeros((4,4)) # alle 4 Permutationen

                    R[0,:] = M[p,:]
                    R[0,0:2] = np.flip(R[0,0:2],0) # Spin Up Flip

                    R[1,:] = M[p,:]
                    R[1,2:4] = np.flip(R[1,2:4],0) # Spin Down Flip

                    R[2,:] = M[p,:]
                    R[2,0:2] = np.flip(R[2,0:2],0) # Spin Up Flip
                    R[2,2:4] = np.flip(R[2,2:4],0) # Spin Down Flip

                    R[3,:] = M[p,:] # kein Flip

                    for l in range(0,4): # Vergleich mit allen Permutationen mit Zustand |j>
                        vzl = np.sign(l-1.5)# -1 oder +1 abhängig von der Anzahl an Permutationen
                        if (R[l,:] - Z[i,:]).T.dot(R[l,:] - Z[i,:]) == 0: # Vergleich etwas umständlich, aber funktioniert
                            J[j,i] = vzp * vzl * 1

print(J)
np.savetxt('Hubb_Ham_Zeit_Lösungen/Hubb_Strom_Matrix.txt', J)
