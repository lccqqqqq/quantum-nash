using HomotopyContinuation
using LinearAlgebra
using Symbolics
include("list_of_functions.jl")
H1 = H1; H2 = - H2;
# define two player games with non-commutative payoff matrices
N = 2;
#H1 = pauli(2, 1, 1) * pauli(2, 2, 2); # H1 = X1 * X2
#H2 = pauli(2, 1, 3) + pauli(2, 2, 3); # H2 = Z1 + Z2

# keeping variables H1, H2; udisncomment to renew the payoff matrices


H1 = diagm([1,1,1,1]) + 0.001*rand_herm(4)
H2 = diagm([1,1,1,1]) + 0.001*rand_herm(4)

H1 = kron(rand_herm(2), diagm([1, 1])) + 1e-10*rand_herm(4)
H2 = kron(diagm([1, 1]), rand_herm(2)) + 1e-10*rand_herm(4)

H1 = rand_herm(4);
H2 = rand_herm(4);

@var x[1:2^N] y[1:2^N]
f_vec = [];

# should be 6 real quadratic form constraints
for a in 1:3
    G1 = 1im * com(H1, pauli(N, 1, a));
    #push!(f_vec, sum(G1[i, j] * (x[i] - 1im * y[i]) * (x[j] + 1im * y[j]) for i in 1:4, j in 1:4));
    push!(f_vec, (x'-im*y') * G1 * (x+im*y))
    
    G2 = 1im * com(H2, pauli(N, 2, a));
    #push!(f_vec, sum(G2[i, j] * (x[i] - 1im * y[i]) * (x[j] + 1im * y[j]) for i in 1:4, j in 1:4));
    push!(f_vec, (x'-im*y') * G2 * (x+im*y))
end

# 8 variables. For finite number (0 dimensional) solutions we need 2 more constraints
# normalization
push!(f_vec, x'*x + y'*y - 1);
# phase fixing, let's fix the 1st variable to be real
push!(f_vec, y[1]);

F = System(f_vec);
Ivs = GroupActions(x -> -x)
result = solve(F; start_system = :total_degree)


# compute entanglement entropy of the states 
M = solutions(result; only_nonsingular=true, only_real=true);
M = unique_points(M, atol=1e-5, group_actions=Ivs);
# convert M[i] from R8 to C4
ns = Vector{Vector{ComplexF64}}() # initialize Nash states
for i in 1:size(M)[1]
    m = Vector{ComplexF64}()
    for j in 1:4
        push!(m, M[i][j] + 1im * M[i][j+4])
    end
    push!(ns, m)
end

# expected utility functions for the two players, for each Nash state
eu = zeros(Float64, size(M)[1], 2)
isNE = zeros(Bool, size(M)[1], 2)
hes_eigvals = Array{Vector{Float64}, 2}(undef, size(M)[1], 2)
hes_eigvals = [zeros(Float64, 3) for _ in 1:size(M)[1], _ in 1:2]
hes = [zeros(ComplexF64, 3, 3) for _ in 1:size(M)[1], _ in 1:2 ]
# Local Hamiltonian, summarizing all the payoff functions
Hloc = [H1, H2]
for player = 1:2
    for i = 1:size(M)[1]
        eu[i, player] = real(ns[i]' * Hloc[player] * ns[i]);
        # verify that these states are indeed Nash stationary

        # check Hessian, whether a Nash stationary point is a Nash equilibrium
        M0, lam, isPosDef = compute_hessian(Hloc[player], N, player, ns[i]);
        isNE[i, player] = isPosDef;
        hes_eigvals[i, player] = lam;
        hes[i, player] = M0;
    end
end





ents = Vector{Float64}()
for i in 1:size(M)[1]
    push!(ents, real(ent_ent_v1(ns[i] * ns[i]')))
end

# find equal elements in the code and identify the solutions
# proved that Nash stationary points exists for some fine-tuned initial states

# Let's consider the entanglement of the eigenstates of Hamiltonian

best_v1 = eigvecs(H1)[:, 4];
best_l1 = eigvals(H1)[4];
best_v2 = eigvecs(H2)[:, 4];
best_l2 = eigvals(H2)[4]

ent_best_v1 = ent_ent_v1(best_v1 * best_v1')
ent_best_v2 = ent_ent_v1(best_v2 * best_v2')


