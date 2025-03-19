%%

% Variational method to finding Nash Equilibrium

% sample data
J = 1;
Gm = 1.7;
h = 0.1;
N = 5;
bc = 'PBC';

% initial guess
[H,Hloc] = hIM(J,Gm,h,N,bc);

[spec,lam] = eig(H);
vac = spec(:,1);
gs = vac;
iter = 1;
eps = 1;

% vac = (1-2*rand(2^N,1));
% vac = vac/norm(vac);

% vac = zeros(2^N,1);
% vac(1) = 1;
% iter = 1;
% eps = 1;
while and(or(abs(lam(1))>1e-5,eps>1e-2),iter<20)
% construct K
    K = zeros(2^N);

    for j = 1:N
        for a = 1:3
            g = G(j,a,Hloc);
            assert(isequal(g,g'),'g is not Hermitian')
            m = vac'*g*vac;
            K = K + m*g;
        end
    end
    try
        assert(isempty(null(K)),'K has non-empty kernel')
    catch
        break
    end
    
    [spec,lam] = eig(K,'vector');
    [~,C] = sort(abs(lam));
    lam = lam(C);
    
    % K typically has non-empty kernel, if ground state is used
    % ker = null(K);
    
    % find the vector having maximal overlap with original vacuum state
    % i.e. compute the projection of vac on ker.
    
    vac_new = spec(:,C(1));
    
    
    
    
    

    iter = iter+1;    
    eps = min([abs(norm(vac_new-vac)),abs(norm(vac_new+vac))]); 
    vac = vac_new;
    
    

    
%     
    data = zeros(N,3);
    for j = 1:N
        for a = 1:3
            g = G(j,a,Hloc);
            data(j,a) = real(vac'*g*vac);
        end
    end
    l = sum(sum(data.^2));

    disp(['i=',num2str(iter),' lam=',num2str(lam(1)),' eps=',num2str(eps), ...
    ' l=',num2str(l)]);
    
    % assert(isequal(lam(1),vac'*K*vac),'eq2.7 violation')
end


% check that the state is indeed a Nash eq
%data = zeros(N,3);
for j = 1:N
    for a = 1:3
        g = G(j,a,Hloc);
        
        %data(j,a) = gs'*g*gs;
        try
            assert(norm(vac'*g*vac)<1e-10,'stationary condition not satisfied');
        catch
            disp([j,a])
        end
    end
end

K = zeros(2^N);

for j = 1:N
    for a = 1:3
        g = G(j,a,Hloc);
        assert(isequal(g,g'),'g is not Hermitian')
        m = vac'*g*vac;
        K = K + m*g;
    end
end



% compare the statements..
k1 = 0; k2 = 0; K11 = zeros(2^N,2^N);
for j = 1:N
    for a = 1:3
        K11 = K11 + (vac'*G(j,a,Hloc)*vac)*G(j,a,Hloc);
        
        k2 = k2 + (vac'*G(j,a,Hloc)*vac)*(vac'*G(j,a,Hloc)*vac);
    end
end
k1 = vac'*K11*vac;



l1 = vac'*K*vac;

    K0 = zeros(2^N);

    for j = 1:N
        for a = 1:3
            g = G(j,a,Hloc);
            assert(isequal(g,g'),'g is not Hermitian')
            K0 = K0 + (gs'*g*gs)*g;
        end
    end
    
%%


for j = 1:N
    for a = 1:3
        g = G(j,a,Hloc);
        
    end
end



g12 = G(1,2,Hloc);
g11 = G(1,1,Hloc);
[spec,lam] = eig(g,'vector');

%% tests

% vac'*G(1,2,Hloc)*vac
% vac'*G(2,2,Hloc)*vac
% vac'*G(3,2,Hloc)*vac
% vac'*G(4,2,Hloc)*vac

% randomized state psi
N = 5;
psi = 1-2*rand(2^N,1);
psi = psi/norm(psi);

psi = zeros(2^N,1);
psi(1) = 1;

K = fK(Hloc,psi);
kerK = null(K);
try
    assert(isempty(kerK),'kernel of K is not empty')
catch
    disp(kerK)
end
        

%% functions

function g = G(j,a,Hloc)
% find the matrix g
c = @(A,B) A*B-B*A; % commutator

[Xj,Yj,Zj] = PauliM(size(Hloc,3),j);
switch a
    case 1
        g = 1i*c(Xj,Hloc(:,:,j));
    case 2
        g = 1i*c(Yj,Hloc(:,:,j));
    case 3
        g = 1i*c(Zj,Hloc(:,:,j));
end

end

function K = fK(Hloc,v)
% construct the matrix K
N = size(Hloc,3);
K = zeros(2^N);

for j = 1:N
    for a = 1:3
        g = G(j,a,Hloc);
        assert(isequal(g,g'),'g is not Hermitian')
        m = v'*g*v;
        K = K + m*g;
    end
end

end

function [Xi,Yi,Zi] = PauliM(N,i)
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
