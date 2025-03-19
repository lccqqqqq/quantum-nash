function [Xi,Yi,Zi] = PauliM(N,i)
% create Pauli matrix string with ith operator nontrivial
X = [0,1;1,0];
Y = 1i * [0,-1;1,0];
Z = [1,0;0,-1];
I = eye(2);

Xi = 1;
Yi = 1;
Zi = 1;
for j = 1:N
    if j~=i
        Xi = kron(Xi,I);
        Yi = kron(Yi,I);
        Zi = kron(Zi,I);
    else
        Xi = kron(Xi,X);
        Yi = kron(Yi,Y);
        Zi = kron(Zi,Z);
    end
end

end
