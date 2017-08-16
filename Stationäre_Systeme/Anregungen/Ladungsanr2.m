# Ladungsanregung des Hubbard-4-Elektronen Systems
# Hier: 2 Spin-up, 1 Spin-down

# Zustände:
Z_down = [1,2,3,4];
Z_up_a = load('Zustand_2e4N.txt');
Z_up = cell2mat(struct2cell(Z_up_a));

Z = zeros(3,24); # 24-dim Hilbertraum
for i=1:4;
  for j=1:6;
    Z(:,(i-1)*6 + j)=[Z_down(i);Z_up(:,j)];
  end;
end;

Z
save('Ladungsanr_Zust.txt', 'Z');

# Hamiltonian:
H_diag = 2*eye(24,24)
# Diagonalelemente:
for i=1:24;
  if Z(1,i) == Z(2,i); # Wenn Loch und ein Elektron aufeinander sitzen
    H_diag(i,i) -= 1;
  elseif Z(1,i) == Z(3,i); # Wenn Loch und ein Elektron aufeinander sitzen
    H_diag(i,i) -= 1;
  end;
end;

# Nichtdiagonalelemente:
H_hupf = zeros(24,24);
function f = period_plus(i); # für die periodische Randbedingung
  if i==4; # periodische Randbedingung nach oben
    f=1;
  else;
    f=i+1;
  end;
end;

function f = period_minus(i); # für die periodische Randbedingung
  if i==1; # periodische Randbedingung nach unten
    f=4;
  else;
    f=i-1;
  end;
end;

function f = J_hoch(Z,i); # Ein Atom hochhüpfen
  f = Z;
  Z_p = period_plus(Z(i));
  f(i) = Z_p;
end;

function f = J_runter(Z,i); # Ein Atom runterhüpfen
  f = Z;
  Z_m = period_minus(Z(i));
  f(i) = Z_m;
end;

for k=2:24;
  for l=1:(k-1);
    for i=1:3;
      M = J_hoch(Z(:,k),i);
      M_f = [M(1); flip(M(2:3))];
      N = J_runter(Z(:,k),i);
      N_f = [N(1); flip(N(2:3))];
      if M == Z(:,l);
        H_hupf(k,l) += 1;
      end;
      if M_f == Z(:,l);
        H_hupf(k,l) -= 1;
      end;
      if N == Z(:,l);
        H_hupf(k,l) += 1;
      end;
      if N_f == Z(:,l);
        H_hupf(k,l) -= 1;
      end;
    end;
  end;
end;
    
# H ist symmetrisch:
for j = 2:24;
  for i = 1:(j-1);
    H_hupf(i,j) = H_hupf(j,i);
  end;
end;
      
[v,lambda] = eig(-H_hupf + 10000*H_diag);

lambda
