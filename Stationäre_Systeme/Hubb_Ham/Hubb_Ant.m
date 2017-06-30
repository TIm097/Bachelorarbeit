# Übergangsmöglichkeiten der einzelnen Zustände
# Matrizen mit den Zustandsvektoren für das jeweilige Spin-UP/DOWN System:
SYS = [1,1,1, 0,0,0; 1,0,0, 1,1,0; 0,1,0, 1,0,1; 0,0,1,0,1,1];

# Zustandsvektoren Gesamtsystem:
Z = zeros(8,36);
for i=0:5;
  for j=1:6;
    Z(1:8,i*6+j) = [SYS(1:4,i+1); SYS(1:4,j)];
  end;
end;
Z;

# Funktion für den Hüpfoperator:
function f = J(x);
  f = eye(8,8);
  A = eye(8,8);
  if x==4;
    f(x,1:8) = A(x-3,1:8);
    f(x-3,1:8) = A(x,1:8);
  elseif x==8;
    f(x,1:8) = A(x-3,1:8);
    f(x-3,1:8) = A(x,1:8);
  else;
    f(x,1:8) = A(x+1,1:8);
    f(x+1,1:8) = A(x,1:8);
  end;
end;

C = zeros(36,2);
C(:,1) = linspace(1,36,36); 

for i=1:36;
  for j=1:36;
    if i != j;
      for k=1:4;
        if J(k)*Z(:,i) == Z(:,j); # Hüpfterm im Spin Up System
          C(i,2) +=1;
          break;
        end;
        p=k+4;
        if J(p)*Z(:,i) == Z(:,j); # Hüpfterm im Spin Down System
          C(i,2) +=1;
          break;
        end;
      end;
    end;
  end;
end;

C;
d = sum(C(:,2));
C(:) = 0.25*C(:);
for i = 1:36;
  C(i) = 1/sqrt(C(i));
end;
Cd = diag(C(:));
#save('Hubb_Ant.txt', 'Cd');