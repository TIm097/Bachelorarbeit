# Ein-Elektron-Vier-Niveau-System, Berechnung der Eigenwerte/-vektoren
H = -[0,1,0,1; 1,0,1,0; 0,1,0,1; 1,0,1,0]
save('Ham1.txt','H')
[v,lambda] = eig(H)
u = v(:,1)
save('GZ1.txt','u')
