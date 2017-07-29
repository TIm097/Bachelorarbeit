# 2 Niveau system mit Spin up/down:
H_d = zeros(4,4);
U = 1000;
H_d(1,1) = U;
H_d(4,4) = U;

H_j = [0,1,1,0; 1,0,0,1; 1,0,0,1; 0,1,1,0];

H = H_d - H_j

[v,lambda] = eig(H)