# Ein-Elektron-Vier-Niveau-System, Berechnung der Eigenwerte/-vektoren mit Diagonalhüpfen l
l = 0.2
H = -[0,1,l,1; 1,0,1,l; l,1,0,1; 1,l,1,0]
save('lambda_Ham1.txt','H')
[v,lambda] = eig(H)
u = v(:,1)
save('lambda_GZ1.txt','u')