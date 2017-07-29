# Wie oft ist jeder Gitterplatz bei dem jeweiligen Zustand besetzt?

Z = cell2mat(struct2cell(load('Hubb_Zust.txt')));

M = zeros(36,4);
for i = 1:36;
  for k = 1:4;
    M(i,Z(k,i))+=1;
  end;
end;

z = linspace(1,36,36);

[z;Z]
[transpose(z),M]

save('Hubb_Besetzung.txt','M')