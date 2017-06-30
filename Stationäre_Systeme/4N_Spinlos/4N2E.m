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
H = zeros(6,6);

for i=2:6;
  for j=1:(i-1);
    for k=1:2; # für das Elektron k
      for r=1:2; # Hüpfrichtung r
        M = J(Z(:,i),k,r);
        if M == Z(:,j);
          H(i,j) += 1;
        elseif flip(M) == Z(:,j);
          H(i,j) -= 1; # Vertauschung zweier Fermionen
        end;
      end;
    end;
  end;
end;  

# H ist symmetrisch:
for j = 2:6;
  for i = 1:(j-1);
    H(i,j) = H(j,i);
  end;
end;
H
# Eig:
[v, lambda]=eig(H)