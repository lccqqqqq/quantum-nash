using LinearAlgebra
using Symbolics

function rand_herm(n::Int64)
    """
    Create a random hermitian matrix of complex entries of size n-by-n
    """
    A = ones(n, n) - 2 * rand(Complex{Float64}, n, n)
    # A = ones(n, n) - 2 * rand(Float64, n, n)
    A = 1/2 * (A + A')
    return A
end

function rand_realsym(n::Int)
    """
    Create a real symmetric matrix of size n-by-n
    """
    A = ones(n, n) - 2 * rand(Float64, n, n)
    # A = ones(n, n) - 2 * rand(Float64, n, n)
    A = 1/2 * (A + transpose(A))
    return A
end


function pauli(N::Int64, i::Int64, a::Int64)
    """
    Create a Pauli string I...ISI...I where
    # N is the total number of sites
    # i is the site at which the Pauli matrix is acted
    # a = x, y, z specify which Pauli is acted on this site
    """
    S = 1.0
    for j in 1:N
        if j != i
            S = kron(S, I(2))
        else
            if a == 1
                S = kron(S, [0 1; 1 0])  # Pauli X
            elseif a == 2
                S = kron(S, 1im * [0 -1; 1 0])  # Pauli Y
            elseif a == 3
                S = kron(S, [1 0; 0 -1])  # Pauli Z
            end
        end
    end
    return S
end

function com(A, B)
    """
    The commutator of matrices A and B
    """
    return A * B - B * A
end


function ent_ent_v1(rho::Matrix{ComplexF64})
    rhoA = zeros(ComplexF64, 2, 2)
    for ia in 1:2
        for ja in 1:2
            rhoA[ia, ja] = rho[2*(ia-1)+1, 2*(ja-1)+1] + rho[2*ia, 2*ja]
        end
    end
    lam = eigvals(rhoA)
    # for zero eigenvalues...
    e = 0
    for i = 1:size(lam)[1]
        if lam[i] <= 0 && abs(lam[i]) < 1e-8
            e += 0
        else
            e += -lam[i] * log(lam[i]);
        end
    end

    # for precision 1e-8
    if abs(e) < 1e-8
        e = 0;
    end


    # return -dot(lam, log.(lam))
    return e
end



function dec2bin(n::Int64, min_digit::Int64)
    n = n-1 # start from 0
    T = []
    while n > 0
        if n%2 == 1
            n = (n-1)/2
            push!(T,1)
        else 
            n = n/2
            push!(T,0)
        end
    end
    foreach(print, reverse(T))
end

function compute_hessian(H::Matrix{ComplexF64}, N::Int64, i::Int64, psi::Vector{ComplexF64})
    M = zeros(ComplexF64, 3, 3)
    for a = 1:3
        for b = 1:3
            M[a, b] = psi' * (pauli(N, i, a) * H * pauli(N, i, b) - H * isequal(a, b)) * psi
        end
    end
    # check whether all eigenvalues are positive definite
    lam = real(eigvals(M))
    IsAllPositive = all(x -> x>0, lam);
    return M, lam, IsAllPositive
end
# 
function disp(a)
    return show(IOContext(stdout, :limit => false), "text/plain", a)
end