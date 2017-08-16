# Hubbard-4-Elektronen-System: Berechnung von <psi| S^2 |psi> für alle eigs 

# Zustände:
Z = cell2mat(struct2cell(load('Hubb_Zust.txt')));
z = linspace(1,36,36);

[z;Z]

# Hamiltonian:
H_diag = cell2mat(struct2cell(load('Hubb_Ham_d.txt')));
H_j = cell2mat(struct2cell(load('Hubb_Ham_j.txt')));

# Eigenwerte + Eigenvektoren:
H = -H_j + 1*H_diag;
[ev,lambda] = eig(H);
vgz = ev(:,1); 
v1az = ev(:,2);
ka = zeros(36,1);

# Diagonalelemente sind 0, da Sz |GZ> = 0  

# S-|GZ> Komponenten:
Zspinanr = cell2mat(struct2cell(load('Spinanr_Zustgroß.txt')));

# S- Hilbertraum:
sminus = zeros(16,36);

# Für alle Zustände:

# <GZ| S+S- |GZ> = 0, wie erwartet
# <AZ| S+S- |AZ> = 2, wie erwartet

for p = 1:36;
  v = ev(:,p);
  for i = 1:36;
    Mges = [Z(:,i),[flip(Z(1:2,i)); Z(3:4,i)]]; #negatives Vorzeichen für zweiten Spaltenvektor
    # (Jetzt S- auf die zweite Komponente anwenden) => zwei Zustände aus dem 16 D-Raum
    for k = 1:2;
      l = -(k*2 - 3); #für das negative Vorzeichen
      M = Mges(:,k);
      N = [M(1);flip(M(2:3));M(4)]; #flip 1
      O = [M(1);M(2);flip(M(3:4))]; #flip 2 
      P = [M(1);M(4); M(3); M(2)]; #flip 3
      Q = [M(1);N(2); flip(N(3:4))]; #doppelflip1
      R = [M(1);N(4); N(3); N(2)]; #doppelflip2
      for j=1:16;
        if M == Zspinanr(:,j);
          sminus(j,p) += l*v(i);
        end;
        if N == Zspinanr(:,j);
          sminus(j,p) -= l*v(i);
        end;
        if O == Zspinanr(:,j);
          sminus(j,p) -= l*v(i);
        end;
        if P == Zspinanr(:,j);
          sminus(j,p) -= l*v(i);
        end;
        if Q == Zspinanr(:,j);
          sminus(j,p) += l*v(i);
        end;
        if R == Zspinanr(:,j);
          sminus(j,p) += l*v(i);
        end;
      end;
    end;
  end;
end;  
 
Squad = zeros(36,1);
S = zeros(36,1);
for i = 1:36; 
  Squad(i) = transpose(sminus(:,i))*sminus(:,i); 
  # Lösen der quadratischen Gleichung a = S(S+1):
  S(i) = -0.5 + sqrt(0.25 + Squad(i));
end;

zeta = [transpose(z),Squad,S]
save('Hubb_S_quad.txt', 'zeta');

lambda(3,3)-lambda(1,1)