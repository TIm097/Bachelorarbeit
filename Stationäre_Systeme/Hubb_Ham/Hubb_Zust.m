# Alle mögichen Zustände im Hubbard-4-Elektronen-System

# Zwei-Elektron-Vier-Niveau-System:
SYS_a = load('Zustand_2e4N.txt');
SYS = cell2mat(struct2cell(SYS_a));

# Zustandsvektoren Gesamtsystem:
Z = zeros(4,36);
for i=0:5;
  for j=1:6;
    Z(1:4,i*6+j) = [SYS(1:2,i+1); SYS(1:2,j)];
  end;
end;
z = linspace(1,36,36);
[z;Z]


save('Hubb_Zust.txt', 'Z');