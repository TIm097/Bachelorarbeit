# Hubbard Hamiltonian:
H_j = cell2mat(struct2cell(load('Hubb_Ham/Hubb_Ham_j.txt')));
H_d = cell2mat(struct2cell(load('Hubb_Ham/Hubb_Ham_d.txt')));
Z = cell2mat(struct2cell(load('Hubb_Ham/Hubb_Zust.txt')));
H_0 = 0.3*H_j + 3*H_d;


function a = f(x,t);
  a = [0;x(1)]; # real;imag
end;

x = lsode('f',[0;1],t=linspace(0,4*pi,10));
x