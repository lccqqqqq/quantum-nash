nkx = 100;
kx = linspace(-pi,pi,nkx);
t = 1;
D = 0.3;
mu = -2;
gamma = -0.5;

N = 12; % number of chains in y direction


X = [0,1;1,0]; Y = [0,-1i;1i,0]; Z = [1,0;0,-1]; I = eye(2);
p = @(m,n) proj(m,n,N,'pbc');

Spec = zeros(2*N,nkx);
for nk = 1:nkx
    % construct Hamiltonian
    H = zeros(2*N);
    k = kx(nk);
    for i = 1:N
        H = H - gamma * kron((Z + 1i*X), p(i,i+1));
        H = H - gamma * kron((Z - 1i*X), p(i+1,i));
        H = H + kron(-(2*t*cos(k)+mu)*Z + D*sin(k)*Y, p(i,i));
    end
    [~,lam] = eig(H,'vector');
    lam = sort(lam);
    Spec(:,nk) = lam;
end

plot(kx,Spec)
ylim([-2.1,2.1])
xlim([-pi,pi])
pbaspect([1,1,1])

% 3D case - Weyl semimetal

% add another dimension:
kz = linspace(-pi,pi,nkx);
k0 = pi/3;

% hopping parameters
tx = 1;
ty = 1;
tz = 1;
m = 4;
%%

% same N = 6;
% the mid gap ribbon of edge states
nkx = 301;
kz = linspace(-pi/3,pi/3,nkx);
kx = linspace(-0.2,0.2,nkx);
Spec = zeros(nkx,nkx,2*N);
for ikx = 1:nkx
    for ikz = 1:nkx
        H = zeros(2*N);
        k1 = kx(ikx);
        k2 = kz(ikz);
        for i = 1:N
            H = H + kron((2*tz*(cos(k0)-cos(k2))+m*(2-cos(k1)))*Z + 2*tx*sin(k1)*X, p(i,i));
            H = H + kron(-m/2*Z+1i*ty*Y, p(i+1,i));
            H = H + kron(-m/2*Z-1i*ty*Y, p(i,i+1));
        end
        [~,lam] = eig(H,'vector');
        lam = sort(lam);
        Spec(ikx,ikz,:) = lam;
    end
end
[Xmesh,Ymesh] = meshgrid(kx,kz);
figure()
hold on
for i = N:N+1
    surf(Xmesh,Ymesh,Spec(:,:,i),EdgeColor="none",FaceAlpha=1);
end
hold off
grid off


% Energy function as
% E = zeros(size(Xmesh));
% for i = 1:size(Xmesh,1)
%     for j = 1:size(Xmesh,2)
%         E(i,j) = sqrt(4)


% slice at some kz
%%
nkx = 151;
kz_slice = 0;
kx = linspace(-1,1,nkx);
kz = linspace(-1,1,nkx);
[Xmesh,Ymesh] = meshgrid(kx,kz);
E = @(kx,kz,ky) sqrt(4*sin(kx)^2 + 4*sin(ky)^2 + (2*(cos(k0)-cos(kz)) + m*(2-cos(kx)-cos(ky))));
kyset = 0:pi/(N+1):pi*(N-1)/N;

Bands = zeros(nkx,nkx,N);
for i = 1:nkx
    for j = 1:nkx
        for l = 1:N
            Bands(i,j,l) = real(E(kx(i),kz(j),0));
        end
    end
end
figure()
hold on
for i = 1
    surf(Xmesh,Ymesh,Bands(:,:,i));
    %surf(Xmesh,Ymesh,-Bands(:,:,i));
end
hold off

%%
kz_slice = 0;
nk_slice = 501;
kx = linspace(-pi,pi,nk_slice);
spec_slice = zeros(nk_slice,2*N);

p = @(m,n) proj(m,n,N,'obc');
for j = 1:nk_slice
    H = zeros(2*N);
    k1 = kx(j);
    for i = 1:N
        H = H + kron((2*tz*(cos(k0)-cos(kz_slice))+m*(2-cos(k1)))*Z + 2*tx*sin(k1)*X, p(i,i));
        H = H + kron(-m/2*Z+1i*ty*Y, p(i+1,i));
        H = H + kron(-m/2*Z-1i*ty*Y, p(i,i+1));
    end
    [~,lam] = eig(H,'vector');
    lam = sort(lam);
    spec_slice(j,:) = lam;
end

figure()
plot(kx,spec_slice);
axis tight

%% finding mid-gap ribbon




function p = proj(m,n,N,bc)
if (m>0) && (m<=N) && (n>0) && (n<=N)
    p = zeros(N);
    p(m,n) = 1;
else % for open boundaries
    switch bc
        case 'obc'
            p = zeros(N);
        case 'pbc'
            while m>N
                m = m-N;
            end
            while m<0
                m = m+N;
            end
            while n>N
                n = n-N;
            end
            while n<=0
                n = n+N;
            end
            p = zeros(N);
            p(m,n) = 1;
    end

end
end

