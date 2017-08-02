# Eigs des Hamiltonian im Hubbardmodell:

# Diagonalelemente:
H_a = load('Hubb_Ham/Hubb_Ham_d.txt');
H_d = cell2mat(struct2cell(H_a));

# Nebendiagonalelemente:
H_a = load('Hubb_Ham/Hubb_Ham_j.txt');
H_j = cell2mat(struct2cell(H_a));

# Zustände:
Z = cell2mat(struct2cell(load('Hubb_Ham/Hubb_Zust.txt')));
z = linspace(1,36,36);
[z;Z]

#Prüfen der Hamiltonians:

# Ohne WW-Term:
#H_hupf = H_j;
#[v,lambda] = eig(H_hupf)

# Ohne Hüpf-Term:
#H_ww = H_d;
#[v,lambda] = eig(H_ww)


# Eigs für verschiedene WW-Werte:
a = 2000; # Anzahl Werte  (U/J=a/100)
gzv = zeros(a,4);
gzeig = zeros(a,9);

gzv(:,1) = linspace(0.01,a/100,a);
gzeig(:,1) = linspace(0.01,a/100,a);


# Für den Plot von Eigenvektorwerten, Eigenwerten:
for o = 1:a;
  H = -H_j + gzv(o,1)*H_d;
  [v, lambda] = eig(H);
  gzv(o,2) = v(11,1)^2 * 2; #Beitrag aller diag-symmetrischer min-Energie Beiträge
  gzv(o,3) = v(6,1)^2 * 4; #Beitrag aller para-symmetrischer min-Energie Beiträge
  gzv(o,4) = 1 -gzv(o,2) -gzv(o,3);
  
  gzeig(o,2:9) = diag(lambda)(1:8);
end;

H = -H_j + 4*H_d; # realistische Werte
[v,lambda] = eig(H);
[transpose(z),diag(lambda)]
v0 = v(:,1);
[transpose(z),v0];
save('Hubb_Eig_Ergebn/Hubb_Ham_psi0.txt', 'v0')

save('Hubb_Eig_Ergebn/Hubb_eplot.txt', 'gzeig');
save('Hubb_Eig_Ergebn/Hubb_gzv.txt', 'gzv');

# Für den Plot der Differenz zwischen Grundzustand und erstem angeregten Zustand:

a = 1000; # Anzahl Werte  (U=a/10)
gzdiff = zeros(a,3);
gzdiff(:,1) = linspace(0.1,a/10,a);

for o = 1:a;
  H = -H_j + gzdiff(o,1)*H_d;
  [v, lambda] = eig(H);
  gzdiff(o,2) = lambda(2,2) - lambda(1,1);
  gzdiff(o,3) = lambda(2,2); 
end;

save('Hubb_Eig_Ergebn/Hubb_anr_diff.txt','gzdiff');

# Für den Plot der Energiedifferenz zwischen den ersten beiden Zuständen mit S=0 (1 und 3):

a = 100; # Anzahl Werte  (U=a/10)
ediff = zeros(a,2);
ediff(:,1) = linspace(0.1,a/10,a);

for o = 1:a;
  H = -H_j + ediff(o,1)*H_d;
  [v, lambda] = eig(H);
  ediff(o,2) = lambda(3,3) - lambda(1,1);
end;

save('Hubb_Eig_Ergebn/Hubb_ediff.txt','ediff');

# Für die Tabelle der ersten 6 Eigenenergien in Abh. von U:

a = 20; # Anzahl Werte (U=a)
eigU = zeros(a,7);
eigU(:,1) = linspace(1,a,a);

for o = 1:a;
  H = -H_j + eigU(o,1)*H_d;
  [v, lambda] = eig(H);
  eigU(o,2:7) = diag(lambda)(1:6);
end;

save('Hubb_Eig_Ergebn/Hubb_eigU.txt','eigU')

# Für die Lokalisierung des Knicks (mit differenzenquotient):

a = 200; # Anzahl Werte (U=a)
diff = zeros(a,2);
lin = linspace(2.4493,2.4496,a);
diff(:,1) = lin;

for o = 1:a;
  H = -H_j + diff(o,1)*H_d;
  [v, lambda] = eig(H);
  diff(o,2) = lambda(4,4);
end;

for d = 1:3; # zwei 'Generationen' des Differenzenquotienten
  diff_neu = zeros(a-d*2,2); # zwei Werte am Rand verschwinden
  I_1 = diff(1,1); # untere Grenze
  I_2 = diff(a-(d-1)*2,1); # obere Grenze
  I = (I_2-I_1)/(a-(d-1)*2)/2; # halbe Intervalllänge
  diff_neu(:,1) = linspace(I_1+I,I_2-I,a-d*2);
  for u = 1:(a-d*2);
    diff_neu(u,2) = (diff(u,2)-diff(u+1,2))/(diff(u+1,1)-diff(u,1)); # Differenzenquotient
  end;
  diff = diff_neu;
end;

diff