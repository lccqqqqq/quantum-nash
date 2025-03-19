function g = fG(j,a,Hloc)
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