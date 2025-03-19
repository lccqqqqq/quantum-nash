function [Hloc] = HlocHeis(J1,J2,N,i,flagLocalHamiltonian)
% create ith local term for Heisenberg model
% priorotizing either summation or NE

X = [0,1;1,0];
Y = 1i * [0,-1;1,0];
Z = [1,0;0,-1];
I = eye(2);

switch flagLocalHamiltonian
    case 0 % prioritizing decomposition
        a = 1/2;
    case 1 % prioritizing state being NE
        a = 1;
end

Hloc = zeros(2^N,2^N);
J = [J1,J2];
% n^th NN interaction
for iNN = 1:2
    jl = mod(i-1-iNN,N)+1;
    jr = mod(i-1+iNN,N)+1;
    
    tmpl = cell(1,3);
    for l = 1:3
        tmpl{l} = 1;
    end
    for l = 1:N
        switch l
            case {jl,i}
                tmpl{1} = kron(tmpl{1},X);
                tmpl{2} = kron(tmpl{2},Y);
                tmpl{3} = kron(tmpl{3},Z);
            otherwise
                tmpl{1} = kron(tmpl{1},I);
                tmpl{2} = kron(tmpl{2},I);
                tmpl{3} = kron(tmpl{3},I);
        end
    end
    Hloc = Hloc - a*J(iNN)*(tmpl{1}+tmpl{2}+tmpl{3});
    
    tmpr = cell(1,3);
    for l = 1:3
        tmpr{l} = 1;
    end
    for l = 1:N
        switch l
            case {jr,i}
                tmpr{1} = kron(tmpr{1},X);
                tmpr{2} = kron(tmpr{2},Y);
                tmpr{3} = kron(tmpr{3},Z);
            otherwise
                tmpr{1} = kron(tmpr{1},I);
                tmpr{2} = kron(tmpr{2},I);
                tmpr{3} = kron(tmpr{3},I);
        end
    end
    Hloc = Hloc - a*J(iNN)*(tmpr{1}+tmpr{2}+tmpr{3});
end

end