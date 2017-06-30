# Berechnung von <GZ| S^2 |GZ> in der Spinanregung, erwartet S(S+1) = 2

# <GZ| (Sz^2 + Sz) |GZ>  = 2
# es fehlt:  <GZ| S-S+ |GZ>

# Spinanregung-Hilbertraum:
Z = cell2mat(struct2cell(load('Spinanr_Zust.txt')));
vgz = cell2mat(struct2cell(load('Spinanr_gzv.txt'))); #Grundzustand

# S+ Hilbertraum:
Zplus = [1;2;3;4];
splus = 0;

function perm = f(Z);
  l = length(Z)-1;
  if l == 1;
    perm = Z;
  else M = zeros(l+1, factorial(l)); # +1 für das VZ
    M(l+1, 1:factorial(l-1)) = Z(l+1); # pos VZ
    M(l+1, (factorial(l-1)+1):factorial(l)) = -Z(l+1); # neg VZ
    for i = 1:l;
      N = Z;
      N(i) = Z(1); # max eine Permutation
      M(1,(factorial(l-1)*(i-1)+1):(factorial(l-1)*i)) = Z(i);
      K = f(N(2:(l+1))); # Rekursion
      M(2:l,(factorial(l-1)*(i-1)+1):(factorial(l-1)*i)) = K(1:(l-1),:);
      for u = 1:factorial(l-1);
        M(l+1,factorial(l-1)*(i-1)+u) *= K(l,u);
      end;
    end;
    perm = M;
  end;
end;

for i = 1:16;
  M = f([Z(:,i);1]); #Alle Permutationen, insgesamt 24
  for k = 1:24; #Anzahl Permutation
    if M(1:4,k) == Zplus;
      splus += M(5,k)*vgz(i);
      break;      
    end;
  end;
end;

#  <GZ| S-S+ |GZ> = 
splus*splus # = 0, wie erwartet
