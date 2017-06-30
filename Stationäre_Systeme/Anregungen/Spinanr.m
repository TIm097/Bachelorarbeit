# Spinanregung des Hubbard-4-Elektronen Systems
# Hier: 1 Spin Down, 3 Spin Up

# Zustände:
# System Spin down:
Z_down = [1,2,3,4];

# System Spin up:
Z_up = [1,1,1,2; 2,2,3,3; 3,4,4,4];

# Kombination:
Z = zeros(4, 16); # (erste Komponente: Spin down, Z(2-4): Spin up)
for i=1:4;
  for j=1:4;
    Z(1,(i-1)*4+j) = Z_down(i);
    Z(2:4,(i-1)*4+j) = Z_up(:,j);
  end;
end;
Z
save('Spinanr_Zust.txt', 'Z');


# Hamiltonian:
# Diagonalelemente:
H_diag = zeros(16,16);
for i=1:16;
  if 10-sum(Z(2:4,i)) != Z(1,i); #Bedingung für zwei e-(up,down) auf einem Atom
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
      Z(k) += direction;
    end;
  else; #direction = -1
    if Z(k) == 1;
      Z(k) = 4;
    else;
      Z(k) += direction;
    end;
  end;
  f = Z;
end;

H_hupf = zeros(16,16);
for i=2:16;
  for j=1:(i-1);
    s=0; # Abbruchparameter, falls H schon erhöht wurde
    #Das Spin down e- kann immer hüpfen
    D = Z(:,i);
    if Z(1,i) == 1; # nach unten hüpfen
      D(1) = 4;
      if D == Z(:,j);
        H_hupf(i,j) += 1;
        s = 1;
      end;
    else;
      D(1) -= 1; 
      if D == Z(:,j);
        H_hupf(i,j) += 1;
        s = 1;
      end;
    end;
    if s==0; # Abbruchbedingung
      D = Z(:,i);
      if D(1) == 4; # nach oben hüpfen
        D(1) -= 3;
        if D == Z(:,j);
          H_hupf(i,j) += 1;
          s = 1;
        end;
      else;
        D(1) += 1; 
        if D == Z(:,j);
          H_hupf(i,j) += 1;
          s = 1;
        end;
      end;
      if s==0; # Abbruchbedingung
        for k=2:4;
          for direction=1:2;
            direction = direction*2 -3; #-1,1            
            M = J(Z(:,i),k,direction);
            if M(2)!=M(3);
              if M(2)!=M(4);
                if M(3)!=M(4);
                  N = [M(1);flip(M(2:3));M(4)]; #flip 1
                  O = [M(1);M(2);flip(M(3:4))]; #flip 2 
                  P = [M(1);M(4); M(3); M(2)]; #flip 3
                  Q = [M(1);N(2); flip(N(3:4))]; #doppelflip1
                  R = [M(1);N(4); N(3); N(2)]; #doppelflip2
                  if M == Z(:,j);
                    H_hupf(i,j) += 1;
                  end;
                  if N == Z(:,j);
                    H_hupf(i,j) -= 1;
                  end;
                  if O == Z(:,j);
                    H_hupf(i,j) -= 1;
                  end;
                  if P == Z(:,j);
                    H_hupf(i,j) -= 1;
                  end;
                  if Q == Z(:,j);
                    H_hupf(i,j) += 1;
                  end;
                  if R == Z(:,j);
                    H_hupf(i,j) += 1;
                  end; 
                end;
              end;
            end;
          end;
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

H = H_hupf + 10*H_diag;
[v,lambda] = eig(H);
gzv = v(:,1);
save('Spinanr_gzv.txt','gzv'); # Für S^2

# Verschiedene WW-Werte:
a = 1000; # Anzahl Werte  (U/J=a)
gzeig = zeros(a,2);

gzeig(:,1) = linspace(1,a,a);

for o=1:a;
  H = H_hupf + gzeig(o,1)*H_diag;
  [v,lambda] = eig(H);
  gzeig(o,2) = lambda(1,1);
end;

save('Spinanr_gze.txt','gzeig')