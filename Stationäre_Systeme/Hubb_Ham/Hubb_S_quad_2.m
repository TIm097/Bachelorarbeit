# Hubbard-4-Elektronen-System: Berechnung von <psi| S^2 |psi> für alle eigs 

# S_z wird 0 
# nur S_+ |Psi> zu berechnen

# Zustände:
Z_alt = cell2mat(struct2cell(load('Hubb_Zust.txt'))); # 36D 
z_alt = linspace(1,36,36);
[z_alt;Z_alt]

Z_neu = cell2mat(struct2cell(load('Spinanr_Zust.txt'))); # 16D
z_neu = linspace(1,16,16);
[z_neu;Z_neu]

# Hamiltonian:
H_diag = cell2mat(struct2cell(load('Hubb_Ham_d.txt')));
H_j = cell2mat(struct2cell(load('Hubb_Ham_j.txt')));

# Eigenwerte + Eigenvektoren:
H = -H_j + 4*H_diag;
[ev,lambda] = eig(H);
vgz = ev(:,1); 

# S+ Matrix:
Splus = zeros(16,36);

# Vorzeichenfunktion:
function f = V(Up); # Up hat Länge drei, Funktion liefert Vorzeichen von Permutierung für Sortierung
  U = Up;
  f = 1;
  for n = 1:2;
    for i = 1:(3-n);
      if U(i) > U(i+1); # dann sortieren
        # Tauschen
        o = U(i);
        U(i) = U(i+1);
        U(i+1) = o;
        # Vorzeichen
        f*= -1;
      end;
    end;
  end;
end;      
      
for i = 1:36; # Linkomb des Eigenzustands besteht aus 36 Basis-Zuständen (36D)
  for k = 3:4; # Beide Down-e- können geflippt werden
    if Z_alt(k,i) != Z_alt(1,i); # nur wenn auf dem Gitterplatz kein Up-e- sitzt
      if Z_alt(k,i) != Z_alt(2,i); # nur wenn auf dem Gitterplatz kein Up-e- sitzt
        x_h = 10 - Z_alt(1,i) - Z_alt(2,i) - Z_alt(k,i); # Ort des Spin-Up-Lochs 
        x_e = sum(Z_alt(3:4,i))-Z_alt(k,i); # Ort des Spin-Down-Elektrons
        M = [x_e;x_h];
        vzk = k*2 -7; # +- Vorzeichen durch k
        Up = [Z_alt(1:2,i);Z_alt(k,i)];
        vzp = V(Up);
        for j = 1:16; # Linkomb des Endzustands besteht aus 16 Basis-Zuständen (16D)
          if M == Z_neu(:,j);
            Splus(j,i) += vzp*vzk;
          end;
        end;
      end;
    end;
  end;
end;

save('Hubb_S_plus.txt', 'Splus');

for k = 1:36;
  transpose(Splus*ev(:,k))*(Splus*ev(:,k));
end;

Splus;

Splus(1:16,19:36)

transpose(Splus)*Splus;

# Spin Für die ersten m Eigenvektoren in Abhängigkeit von U, für Tabelle:
a = 10; #(U/J=a/2)
m = 6;
l = linspace(0.5,a/2,a);
spin = zeros(a,m+1);
spin(:,1) = l; 

for o = 1:a;
  H = -H_j + spin(o,1)*H_diag;
  [v,lambda] = eig(H);
  S_quad = zeros(m,1);
  for k = 1:m;
    S_quad(k) = transpose(Splus*v(:,k))*(Splus*v(:,k));
    # Lösen der quadratischen Gleichung a = S(S+1):
    spin(o,k+1) = -0.5 + sqrt(0.25 + S_quad(k));
  end;
end;

save('Hubb_S_tab.txt','spin');

# Spin Für die Bestimmung des ersten Knicks:
a = 100; #(U/J= 2.4 + a/1000/10)
l = linspace(2.44945,2.44955,a);
spin = zeros(a,2);
spin(:,1) = l; 

for o = 1:a;
  H = -H_j + spin(o,1)*H_diag;
  [v,lambda] = eig(H);
  S_quad = transpose(Splus*v(:,4))*(Splus*v(:,4));
  # Lösen der quadratischen Gleichung a = S(S+1):
  spin(o,2) = -0.5 + sqrt(0.25 + S_quad);
end;


