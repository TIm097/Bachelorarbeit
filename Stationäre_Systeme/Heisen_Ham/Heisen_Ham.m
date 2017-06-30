# Heisenberg-4-Elektronen-System:
# Zustände:
Z = [1,1,1,3,2,3; 2,3,4,4,4,2; 3,2,3,1,1,1; 4,4,2,2,3,4]

# Hamiltonian:

function f = period(i); # für die periodische Randbedingung
  if i==4; # periodische Randbedingung
    f=1;
  else;
    f=i+1;
  end;
end;

function f = Sz(Z,i); 
  i_p = period(i); 
  for d=1:4; # Womit ist Platz i/i+1 besetzt?
    if i==Z(d);
      u=d;
    elseif i_p==Z(d);
      v=d;
    end;
  end;
  if u<3; # up
    if v<3; # up/up
      f = 0.25;
    else; # up/down
      f = -0.25;
    end;
  else; # down
    if v<3; # down/up 
      f = -0.25;
    else;
      f = 0.25; # down/down
    end;
  end;
end;

H = zeros(6,6);

for k=1:6; #Zustand
  H(k,k) -= 1; # Konstante!!
  for i=1:4;
    H(k,k) += Sz(Z(:,k),i);
  end;
end;

function f = Spm(Z,i);
  i_p = period(i);  
  for d=1:4; # Womit ist Platz i/i+1 besetzt?
    if i==Z(d);
      u=d;
    elseif i_p==Z(d);
      v=d;
    end;
  end;
  Z(v)=i; # Tausch der Spins
  Z(u)=i_p;
  f = Z;
end;

for k=2:6;
  M = zeros(16,4); # (4 Spintauschmöglichkeiten * 4 Flipmöglichkeiten, 4 Niveaus)
  for i=1:4;
    M(1:4,i) = Spm(Z(:,k),i);
    M(5:8,i) = [flip(M(1:2,i));M(3:4,i)];
    M(9:12,i) = [M(1:2,i);flip(M(3:4,i))];
    M(13:16,i) = [flip(M(1:2,i));flip(M(3:4,i))];
  end;
  for l=1:(k-1);
    for i=1:4;
      if M(1:4,i) == Z(:,l);
        H(k,l) += 1/2;
      end;
      if M(5:8,i) == Z(:,l);
        H(k,l) -= 1/2; # einfacher Fermiontausch
      end;
      if M(9:12,i) == Z(:,l);
        H(k,l) -= 1/2; # einfacher Fermiontausch
      end;
      if M(13:16,i) == Z(:,l);
        H(k,l) += 1/2; # doppelter Fermiontausch
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
[v, lambda] = eig(H)


                                                                    