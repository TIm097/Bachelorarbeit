# Hubbard-4-Elektronen-System: Diagonelelemente des Hamiltonian

# Zustandsvektoren:
Z = cell2mat(struct2cell(load('Hubb_Zust.txt')));
z = linspace(1,36,36);

[z;Z]

# Hamiltonian (Elemente der lokalen WW):
H_d = zeros(36,36);
for i=1:36;
  for k=1:2;
    if Z(k,i)==Z(3,i); # Vgl zwischen beiden Systemen
      H_d(i,i) += 1; 
    elseif Z(k,i)==Z(4,i);
      H_d(i,i) += 1; 
    end;
  end;
end;
H_d

save('Hubb_Ham_d.txt', 'H_d');