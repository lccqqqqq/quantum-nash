function [H] = locInt2(dij,op,N)
% Creating 2-local terms for N-site spin-1/2 Hamiltonians
% specifically for N-th neighbour interactions
% without numerical factor
% positive sign
I = eye(2);
H = zeros(2^N);
for i = 1:N
    if i > N-dij
        j = i-N+dij;
    else
        j = i+dij;
    end
    h = 1;
    for l = 1:N
        if and(l~=i,l~=j)
            h = kron(h,I);
        else
            h = kron(h,op);
        end
    end
    H = H + h;
end

end