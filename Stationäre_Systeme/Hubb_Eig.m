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
a = 250; # Anzahl Werte  (U/J=a/10)
gzv = zeros(a,4);
gzeig = zeros(a,2);

gzv(:,1) = linspace(0.1,a/10,a);
gzeig(:,1) = linspace(0.1,a/10,a);


# Für den Plot von Eigenvektorwerten, Eigenwerten:
for o = 1:a;
  H = H_j + gzv(o,1)*H_d;
  [v, lambda] = eig(H);
  gzv(o,2) = v(11,1)^2 * 2; #Beitrag aller diag-symmetrischer min-Energie Beiträge
  gzv(o,3) = v(6,1)^2 * 4; #Beitrag aller para-symmetrischer min-Energie Beiträge
  gzv(o,4) = 1 -gzv(o,2) -gzv(o,3);
  
  gzeig(o,2) = lambda(1,1);
end;

H = 0.3*H_j + 3*H_d; # realistische Werte
[v,lambda] = eig(H);
diag(lambda);
[transpose(z),v(:,1)];

save('Hubb_Eig_Ergebn/Hubb_gze.txt', 'gzeig');
save('Hubb_Eig_Ergebn/Hubb_gzv.txt', 'gzv');

# Für den Plot der Differenz zwischen Grundzustand und erstem angeregten Zustand:

a = 1000; # Anzahl Werte  (U=a/10)
gzdiff = zeros(a,3);
gzdiff(:,1) = linspace(0.1,a/10,a);

for o = 1:a;
  H = H_j + gzdiff(o,1)*H_d;
  [v, lambda] = eig(H);
  gzdiff(o,2) = lambda(2,2) - lambda(1,1);
  gzdiff(o,3) = lambda(2,2); 
end;

save('Hubb_Eig_Ergebn/Hubb_anr_diff.txt','gzdiff');




