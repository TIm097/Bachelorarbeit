# Spinanregung des Hubbard-4-Elektronen Systems
# Hier: 1 Spin Down, 3 Spin Up
# -> enfacher mit einem Loch, einem Elektron

Z_down = [1,2,3,4];
Z_up_loch = [1,2,3,4];

# Kombination:
Z = zeros(2, 16); # (Z(1): Spin down, Z(2): Spin up Loch)
for i=1:4;
  for j=1:4;
    Z(1,(i-1)*4+j) = Z_down(i);
    Z(2,(i-1)*4+j) = Z_up_loch(j);
  end;
end;
Z
save('Spinanr_Zust.txt', 'Z');

# Hamiltonian:
# Diagonalelemente:
H_diag = zeros(16,16);
for i=1:16;
  if Z(1,i) != Z(2,i); #Bedingung für Loch und e- nicht aufeinander -> Potential
    H_diag(i,i) = 1;
  end;
end;
H_diag;

# Hüpfterme:
function f = J(Z,k,direction); #Hüpfoperator
  if direction == 1;
    if Z(k) == 4;
      Z(k) = 1;
    else;
      Z(k) += 1;
    end;
  else; #direction == -1
    if Z(k) == 1;
      Z(k) = 4;
    else;
      Z(k) -= 1;
    end;
  end;
  f = Z;
end;

H_hupf = zeros(16,16);
for i=2:16; # Zeilen
  for j=1:(i-1); # Spalten
    for k=1:2; # beide Elektronen
      for direction=1:2; # beide Hüpfrichtungen
        direction = direction*2 -3; #-1,1
        M = J(Z(:,j),k,direction);
        if M == Z(:,i); # Vergleich
          H_hupf(i,j) += 1;
        end;
      end;
    end;
  end;
end;
      
# H ist symmetrisch:
for j = 2:16;
  for i = 1:(j-1);
    H_hupf(i,j) = H_hupf(j,i);
  end;
end;

H = -H_hupf + 4*H_diag

[v,lambda] = eig(H);
gzv = v(:,1);
save('Spinanr_gzv.txt','gzv'); # Für S^2

# Verschiedene WW-Werte:
a = 2000; # Anzahl Werte  (U/J=a/100)
gzeig = zeros(a,2);

gzeig(:,1) = linspace(1,a/100,a);

for o=1:a;
  H = H_hupf + gzeig(o,1)*H_diag;
  [v,lambda] = eig(H);
  gzeig(o,2) = lambda(1,1);
end;

save('Spinanr_gze.txt','gzeig')