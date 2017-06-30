# Hubbard-4-Elektronen-System: Diagonelelemente des Hamiltonian

# Zustandsvektoren:
Z_a = load('Hubb_Zust.txt');
Z = cell2mat(struct2cell(Z_a));

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