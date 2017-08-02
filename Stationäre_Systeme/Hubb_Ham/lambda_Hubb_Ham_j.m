# Hubbard-4-Elektronen-System: Nichtdiagonelelemente des Hamiltonian mit Diagonalhüpfen

lambda = 0.2
# Zustandsvektoren:
Z = cell2mat(struct2cell(load('Hubb_Zust.txt')));
z = linspace(1,36,36);

[z;Z]

# Hüpfoperator (für das Elektron k, Hüpfrichtung r, System sy):
function f = J(Z,k,r,sy);
  r = r-1;
  f = Z;
  ksy = k+sy;
  if r > 0;
    s = (r*2)-3; # -1 oder +1 
    if s==1;
      if Z(ksy)==4; #Randbedingung 
        f(ksy) -= 3;
      else;
        f(ksy) += s;
      end;  
    else;
      if Z(ksy)==1; #Randbedingung 
        f(ksy) += 3;
      else;
        f(ksy) += s;
      end;
    end;
  else;
    if Z(ksy) == 1;
      f(ksy) = 3;
    end;
    if Z(ksy) == 3;
      f(ksy) = 1;
    end;
    if Z(ksy) == 4;
      f(ksy) = 2;
    end;
    if Z(ksy) == 2;
      f(ksy) = 4;
    end;
  end;    
end;

# Hamiltonian(Hüpfterme):
H_j = zeros(36,36);

for i=2:36;
  for j=1:(i-1);
    for sys=1:2; # für das System sy
      for k=1:2; # für das Elektron k
        for r=1:3; # Hüpfrichtung r
          l = 1; #Faktor
          if r==1;
            l = lambda;
          end; 
          sy = sys*2 -1; # 1 oder 3
          M = J(Z(:,i),k,r,sy-1);
          N = M;
          N(sy:sy+1) = flip(M(sy:sy+1));
          if M == Z(:,j);
            H_j(i,j) += l;
          elseif N == Z(:,j);
            H_j(i,j) -= l; # Vertauschung zweier Fermionen
          end;
        end;
      end;
    end;
  end;
end;

# H ist symmetrisch:
for j = 2:36;
  for i = 1:(j-1);
    H_j(i,j) = H_j(j,i);
  end;
end;
H_j;

save('lambda_Hubb_Ham_j.txt', 'H_j');

H = -H_j
[lambda,v] = eig(H) 