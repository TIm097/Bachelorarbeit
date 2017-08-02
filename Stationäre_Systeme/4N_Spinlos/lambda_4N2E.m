# Zwei-Elektronen-Vier-Niveau-System, Berechnung der Eigenwerte/ -vektoren mit Diagonalhüpfen
# Zustandsmatrix:

lambda = 0.2
Z = cell2mat(struct2cell(load('Zustand_2e4N.txt')))

# Hüpfoperator (für das Elektron k, Hüpfrichtung r):
function f = J(Z,k,r);
  r = r-1;
  f = Z;
  if r > 0;
    s = (r*2)-3; # -1 oder +1
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
  else;
    if Z(k) == 1;
      f(k) = 3;
    end;
    if Z(k) == 3;
      f(k) = 1;
    end;
    if Z(k) == 4;
      f(k) = 2;
    end;
    if Z(k) == 2;
      f(k) = 4;
    end;
  end;    
end;

# Hamiltonian:
H_J = zeros(6,6);

for i=2:6;
  for j=1:(i-1);
    for k=1:2; # für das Elektron k
      for r=1:3; # Hüpfrichtung r
        M = J(Z(:,i),k,r);
        l = 1; #Faktor
        if r==1;
          l = lambda;
        end; 
        if M == Z(:,j);
          H_J(i,j) += l;
        elseif flip(M) == Z(:,j);
          H_J(i,j) -= l; # Vertauschung zweier Fermionen
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

H = -H_J

# Eig:
[v, lambda]=eig(H)