# Grenzfall des Hubbard-4-Elektronen-Systems: Approximationsgerade

H_j = cell2mat(struct2cell(load('Hubb_Ham_j.txt')));
H_d = cell2mat(struct2cell(load('Hubb_Ham_d.txt')));

# Verschiedene WW-Werte:
a = 40000; # Anzahl Werte  (U/J=a)
EW = zeros(a,2);
EW(:,1) = linspace(a,a+10000,a);

for o = 1:a;
  H = -H_j + EW(o,1)*H_d;
  [v, lambda] = eig(H);
  EW(o,2) = lambda(1,1);
end;

save('Hubb_Grenz_lin.txt', 'EW');