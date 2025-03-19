function [H,Hloc] = hIM(J,Gm,h,N,bc,flg)
% create 2^N dimensional matrix representation of the global Hamiltonian of
% Ising model
if nargin == 5
    flg = 0;
end

H = zeros(2^N);
Hloc = zeros(2^N,2^N,N);


I = [1,0;0,1];
X = [0,1;1,0];
% Y = 1i * [0,-1;1,0];
Z = [1,0;0,-1];

for i = 1:N
    switch bc % different boundary conditions
        case 'PBC'
            jr = mod(i,N)+1;
            jl = mod(i-2,N)+1;
        case 'OBC'
            if i == N
                jr = 0;
                jl = N-1;
            elseif i == 1
                jl = 0;
                jr = 2;
            else
                jr = i+1;
                jl = i-1;
            end
    end
    
    % ZZ term for sites i,i+1
    if jr~=0
        tmpr = 1;
        for k = 1:N
            if and(k~=i,k~=jr)
                tmpr = kron(tmpr,I);
            else
                tmpr = kron(tmpr,Z);
            end
        end
        H = H - 1/2*J*tmpr;
    else
        tmpr = zeros(2^N);
    end
    if jl~=0
        tmpl = 1;
        for k = 1:N
            if and(k~=i,k~=jl)
                tmpl = kron(tmpl,I);
            else
                tmpl = kron(tmpl,Z);
            end
        end
        H = H - 1/2*J*tmpl;
    else
        tmpl = zeros(2^N);
    end

    
    % X term for site i
    tmpx = 1;
    for k = 1:N
        if k~=i
            tmpx = kron(tmpx,I);
        else
            tmpx = kron(tmpx,X);
        end
    end
    H = H - Gm*tmpx;
    
    % Z term for site i
    tmpz = 1;
    for k = 1:N
        if k~=i
            tmpz = kron(tmpz,I);
        else
            tmpz = kron(tmpz,Z);
        end
    end
    H = H - h*tmpz;
    
    switch flg
        case 0 % decomposition of global Hamiltonian
            Hloc(:,:,i) = -1/2*J*(tmpl+tmpr) - Gm*tmpx - h*tmpz;
        case 1 % feature minimization of local energy
            Hloc(:,:,i) = -J*(tmpl+tmpr) - Gm*tmpx - h*tmpz;
    end
end


end

