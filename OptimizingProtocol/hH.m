function [H,Hloc] = hH(J1,J2,N,flg)
% create 2^N dimensional matrix representation of the global Hamiltonian of
% Next-nearest-neighbour Heisenberg model

if nargin == 3
    flg = 0;
end

H = zeros(2^N);
Hloc = zeros(2^N,2^N,N);

for j = 1:N
    hj = HlocHeis(J1,J2,N,j,flg);
    H = H + hj;
end


end