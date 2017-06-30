# Heisenberg-4-Elektronen-System:
# Zustände:
Z = [1,1,1,-1,-1,-1; 1,-1,-1,-1,1,1; -1,1,-1,1,-1,1; -1,-1,1,1,1,-1];

function f = period(i); # für die periodische Randbedingung
  if i==4; # periodische Randbedingung
    f=1;
  else;
    f=i+1;
  end;
end;

# Spin-Z-Operator:

function f = Sz(Z,i); 
  i_p = period(i); # = i+1
  if Z(i)+Z(i_p) == 0; #Dann Spins antiparallel
    f = -0.25;
  else;
    f = 0.25;
  end;
end;

H = zeros(6,6);

for k=1:6; #Zustand
  H(k,k) -= 1; # Konstante!!
  for i=1:4;
    H(k,k) += Sz(Z(:,k),i);
  end;
end;
    
# Spin +- Operatoren:

function f = Spm(Z,i);
  i_p = period(i);  
  f = Z;
  f(i) = Z(i_p);
  f(i_p) = Z(i);
end;
  
for k=2:6;
  for l=1:(k-1);
    for i=1:4;
      if Spm(Z(:,k),i) == Z(:,l);
        H(k,l) += 0.5;
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
[v,lambda] = eig(H)