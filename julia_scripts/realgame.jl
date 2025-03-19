using HomotopyContinuation
using LinearAlgebra, IterTools
using CSV, DelimitedFiles

include("list_of_functions.jl")

N = 2;
d = 2^N;
H = [rand_realsym(d) for _ in 1:N]
Z = [1 0; 0 -1]
X = [0 1; 1 0];
# H = [(-1)^N*kron(Z, Z)+1e-10*rand_realsym(d) for _ in 1:N]
# ZZ or XX game, finding Nash states? and the Nash minimal states among them?
# H = [kron(Z, Z) for _ in 1:N]
# H = [kron(X, X) for _ in 1:N]
# h1 = [3 0 0 0; 0 0 0 0; 0 0 5 0; 0 0 0 1]+1e-10*rand_realsym(d);
# h2 = [3 0 0 0; 0 5 0 0; 0 0 0 0; 0 0 0 1]+1e-10*rand_realsym(d);
# H = [h1, h2];

# for real and symmetric payoff functions, we would have only 2 nontrivial conditions
G = real([map(i -> [1im * com(H[i], pauli(N, i, 2))], 1:N)][1])

@var x[1:d]
cstr = [sum(x.*x)-1,] # the normalization condition

for i = 1:N
    global cstr = vcat(cstr, transpose(x) * G[i][1] * x)
end

###########################
# The following: solve for the maximally entangled orbit in S3 stereographical projection

# Now attempt to create a dense sample of Nash variety
# compute bottlenecks

D = length(cstr) # dimension of constaints
@var λ[1:D] μ[1:D]; # lagrange multipliers
∇ = differentiate(cstr, x)
@var y[1:d] # create replicated system
replica_sys = subs(cstr, x => y)
∇y = subs(∇, x => y)
system = [cstr; replica_sys; map(j -> x[j]-y[j]-dot(λ, ∇[:, j]), 1:d); map(j -> x[j]-y[j]-dot(μ, ∇y[:, j]), 1:d)]

result = solve(system)


# pick the smallest bottleneck
bottlenecks = map(s -> (s[1:d], s[d+1:2*d]), real_solutions(nonsingular(result)))
bn_lengths = sort!(map(b -> (norm(b[1]-b[2]), b), bottlenecks), by = a -> a[1])
δ = max(bn_lengths[1][1]/(2*sqrt(d)), 0.1) # grid-size

q = bn_lengths[end][2][1] + (bn_lengths[end][2][2]-bn_lengths[end][2][1])/2
system = [cstr; map(j -> x[j]-q[j]-dot(λ, ∇[:, j]), 1:d)]
result1 = solve(system)

critical_points = sort!(map(c -> (norm(c[1:d]-q), c[1:d]), real_solutions(nonsingular(result1))), by = a -> a[1])
b = critical_points[end][1] # length of bounding box
indices = [i for i in -b:δ/10:b]

# samples obtained by intersecting the variety with linear subspaces in the 4 dimensional hypercube
samples = []
k = d - D; 
@var p[1:d]
for s in IterTools.subsets(1:d, k)
    Ft = [cstr; map(i -> x[s[i]]-p[i]-q[s[i]], 1:k)]
    p₀ = randn(ComplexF64, k)
    S_p₀ = solutions(solve(subs(Ft, p[1:k] => p₀)))
    params = [[indices[jj] for jj in p1] for p1 in Iterators.product(map(j-> 1:length(indices), s)...)][:]

    result = solve(
           Ft,
           S_p₀;
           parameters = p[1:k],
           start_parameters = p₀,
           target_parameters = params,
           transform_result = (_r, _p) -> real_solutions(_r),
           flatten = true
       )
    samples = vcat(samples, result)
end



