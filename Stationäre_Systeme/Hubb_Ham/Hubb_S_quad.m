# Hubbard-4-Elektronen-System: Berechnung von <GZ| S^2 |GZ> und <1AZ| S^2 |1AZ>

# Zustände:
Z_a = load('Hubb_Zust.txt');
Z = cell2mat(struct2cell(Z_a));
Z

# Hamiltonian:
H_diag = cell2mat(struct2cell(load('Hubb_Ham_d.txt')));
H_j = cell2mat(struct2cell(load('Hubb_Ham_j.txt')));

# Eigenwerte + Eigenvektoren:
H = H_diag + 10*H_j;
[v,lambda] = eig(H);
vgz = v(:,1); 
v1az = v(:,2);
ka = zeros(36,1);

# Diagonalelemente sind 0, da Sz |GZ> = 0  

# S-|GZ> Komponenten:
Zspinanr = cell2mat(struct2cell(load('Spinanr_Zust.txt')));

# S- Hilbertraum:
sminus = zeros(16,1);

# Für GZ:
v = vgz;
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
        sminus(j) += l*v(i);
      end;
      if N == Zspinanr(:,j);
        sminus(j) -= l*v(i);
      end;
      if O == Zspinanr(:,j);
        sminus(j) -= l*v(i);
      end;
      if P == Zspinanr(:,j);
        sminus(j) -= l*v(i);
      end;
      if Q == Zspinanr(:,j);
        sminus(j) += l*v(i);
      end;
      if R == Zspinanr(:,j);
        sminus(j) += l*v(i);
      end;
    end;
  end;
end;  

# <GZ| S+S- |GZ> :
transpose(sminus) * sminus # = 0, wie erwartet

# Für 1AZ:
sminus = zeros(16,1);
v = v1az;

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
        sminus(j) += l*v(i);
      end;
      if N == Zspinanr(:,j);
        sminus(j) -= l*v(i);
      end;
      if O == Zspinanr(:,j);
        sminus(j) -= l*v(i);
      end;
      if P == Zspinanr(:,j);
        sminus(j) -= l*v(i);
      end;
      if Q == Zspinanr(:,j);
        sminus(j) += l*v(i);
      end;
      if R == Zspinanr(:,j);
        sminus(j) += l*v(i);
      end;
    end;
  end;
end;  
 
# <GZ| S+S- |GZ> :
transpose(sminus) * sminus # = 2, wie erwartet