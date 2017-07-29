# Grenzfall des Hubbard-4-Elektronen-Systems

H_j = cell2mat(struct2cell(load('Hubb_Ham_j.txt')));
H_d = cell2mat(struct2cell(load('Hubb_Ham_d.txt')));

# Verschiedene WW-Werte:
a = 1000; # Anzahl Werte  (U=a)
EW = zeros(a,2);
EW(:,1) = linspace(1,a,a);

for o = 1:a;
  H = -H_j + EW(o,1)*H_d;
  [v, lambda] = eig(H);
  EW(o,2) = lambda(1,1);
end;

save('Hubb_Grenz.txt', 'EW');