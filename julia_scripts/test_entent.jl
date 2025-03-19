using LinearAlgebra
using HomotopyContinuation
include("list_of_functions.jl")

psi = 1/4 * kron(Matrix{ComplexF64}(I, 2, 2), Matrix{ComplexF64}(I, 2, 2));

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

rho = M_C4[15] * M_C4[15]'
ent_ent_v1(rho)