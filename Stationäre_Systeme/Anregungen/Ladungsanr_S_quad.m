# Berechnung von <GZ| S^2 |GZ> in der Ladungsanregung, erwartet S(S+1) = 0.75

# <GZ| (Sz^2 + Sz) |GZ>  = 0.75
# es fehlt:  <GZ| S-S+ |GZ>

# Ladungsanregung-Hilbertraum:
Z = cell2mat(struct2cell(load('Ladungsanr_Zust.txt')));
vgz = cell2mat(struct2cell(load('Ladungsanr_gzv.txt'))); #Grundzustand

# S+ Hilbertraum:
Zplus = [1,1,1,2;2,2,3,3;3,4,4,4];
splus = zeros(4,1);

function perm = f(Z); # Permutationsfunktion
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

for i = 1:24;
  M = f([Z(:,i);1]); #Alle Permutationen, insgesamt 6
  for k = 1:6; #Anzahl Permutationen
    for u = 1:4;
      if M(1:3,k) == Zplus(:,u);
        splus(u) += M(4,k)*vgz(i);
        break;
      end;
    end;
  end;
end;

#  <GZ| S-S+ |GZ> = 
transpose(splus) * splus # = 0, wie erwartet