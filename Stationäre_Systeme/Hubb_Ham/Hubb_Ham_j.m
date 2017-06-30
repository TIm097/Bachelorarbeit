# Hubbard-4-Elektronen-System: Nichtdiagonelelemente des Hamiltonian

# Zustandsvektoren:
Z_a = load('Hubb_Zust.txt');
Z = cell2mat(struct2cell(Z_a));

# Hüpfoperator (für das Elektron k, Hüpfrichtung r, System sy):
function f = J(Z,k,r,sy);
  s = (r*2)-3; # -1 oder +1 
  ksy = k+sy;
  f = Z;
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
end;

# Hamiltonian(Hüpfterme):
H_j = zeros(36,36);

for i=2:36;
  for j=1:(i-1);
    for sys=1:2; # für das System sy
      for k=1:2; # für das Elektron k
        for r=1:2; # Hüpfrichtung r
          sy = sys*2 -1; # 1 oder 3
          M = J(Z(:,i),k,r,sy-1);
          N = M;
          N(sy:sy+1) = flip(M(sy:sy+1));
          if M == Z(:,j);
            H_j(i,j) += 1;
          elseif N == Z(:,j);
            H_j(i,j) -= 1; # Vertauschung zweier Fermionen
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
H_j

save('Hubb_Ham_j.txt', 'H_j');