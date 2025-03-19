% create JW transformations
Z = @(i) Fermion(i,0,1)+Fermion(i,1,1);
Y = @(i) Fermion(i,0,-1i)+Fermion(i,1,1i);
X = @(i) Fermion([],[],1)-2*Fermion([i,i],[1,0],1);

% creating the ground state of TFIM

N = 4;

kk = pi/N:2*pi/N:(N-1)*pi/N;




