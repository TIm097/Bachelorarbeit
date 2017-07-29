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

# Für die Tabelle der ersten 8 Eigenenergien in Abh. von U:

a = 10; # Anzahl Werte (U=a)
eigU = zeros(a,5);
eigU(:,1) = linspace(1,a,a);

for o = 1:a;
  H = -H_j + eigU(o,1)*H_d;
  [v, lambda] = eig(H);
  eigU(o,2:5) = diag(lambda)(1:4);
end;

save('Hubb_Eig_Ergebn/Hubb_eigU.txt','eigU')
