# Zwei-Elektronen-Vier-Niveau-System, Berechnung der Eigenwerte/ -vektoren
# Zustandsmatrix:
Z_1 = [1,2,3,4]; # Vier entartete Niveaus
Z = zeros(2,6);

a = 1;
for i=1:3;
  for j=(i+1):4;
    Z(1,a) = Z_1(i);
    Z(2,a) = Z_1(j);
    a += 1;
  end;
end;

Z 
save('Zustand_2e4N.txt', 'Z');

# Hüpfoperator (für das Elektron k, Hüpfrichtung r):
function f = J(Z,k,r);
  s = (r*2)-3; # -1 oder +1
  f = Z;
  if s==1;
    if Z(k)==4; #Randbedingung 
      f(k) -= 3;
    else;
      f(k) += s;
    end;  
  else;
    if Z(k)==1; #Randbedingung 
      f(k) += 3;
    else;
      f(k) += s;
    end;
  end;
end;

# Hamiltonian:
H_J = zeros(6,6);

for i=2:6;
  for j=1:(i-1);
    for k=1:2; # für das Elektron k
      for r=1:2; # Hüpfrichtung r
        M = J(Z(:,i),k,r);
        if M == Z(:,j);
          H_J(i,j) += 1;
        elseif flip(M) == Z(:,j);
          H_J(i,j) -= 1; # Vertauschung zweier Fermionen
        end;
      end;
    end;
  end;
end;  

# H ist symmetrisch:
for j = 2:6;
  for i = 1:(j-1);
    H_J(i,j) = H_J(j,i);
  end;
end;
# Vorfaktor J ist negativ

# Diagonalelemente:
a = 1*10^(-15)
H_a = zeros(6,6);
for i = 1:6;
  if Z(1,i) == 1 || Z(1,i) == 3;
      H_a(i,i) += a;
  else;
      H_a(i,i) -=a;
  end;
  if Z(2,i) == 1 || Z(2,i) == 3;
      H_a(i,i) += a;
  else;
      H_a(i,i) -=a;
  end;
end;

H = -H_J + H_a

# Eig:
[v, lambda]=eig(H)