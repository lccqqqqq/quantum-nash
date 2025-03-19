using HomotopyContinuation
using LinearAlgebra, IterTools
using CSV, DelimitedFiles

include("list_of_functions.jl")


# The payoff functions 
h1 = [3 0 0 0; 0 0 0 0; 0 0 5 0; 0 0 0 1];
h2 = [3 0 0 0; 0 5 0 0; 0 0 0 0; 0 0 0 1];
d = 6;
@var x[1:6]

ψ = [x[1], x[2]+1im*x[3], x[4], x[5]+1im*x[6]];
ψT = transpose([x[1], x[2]-1im*x[3], x[4], x[5]-1im*x[6]]);

cstr = [sum(x.*x)-1,
    ψT * 1im*com(h1, pauli(2, 1, 2)) * ψ,
    ψT * 1im*com(h2, pauli(2, 2, 2)) * ψ,
    ψT * com(h1, pauli(2, 1, 1)) * ψ,
    ψT * com(h2, pauli(2, 2, 1)) * ψ];

D = 5; # Length of Constraints
@var λ[1:D] μ[1:D]; 
∇ = differentiate(cstr, x)
@var y[1:6];
χ = [y[1], y[2]+1im*y[3], y[4], y[5]+1im*y[6]];
χT = transpose([y[1], y[2]-1im*y[3], y[4], y[5]-1im*y[6]]);

replica_sys = subs(cstr, x => y)
∇y = subs(∇, x => y)
system = [cstr; replica_sys; map(j -> x[j]-y[j]-dot(λ, ∇[:, j]), 1:6); map(j -> x[j]-y[j]-dot(μ, ∇y[:, j]), 1:6)]

result = solve(system)

writedlm("QPDbottleneckdata1.csv", real_solutions(nonsingular(result)), ',')


bottlenecks = map(s -> (s[1:d], s[d+1:2*d]), real_solutions(nonsingular(result)))
bn_lengths = sort!(map(b -> (norm(b[1]-b[2]), b), bottlenecks), by = a -> a[1])
δ = max(bn_lengths[1][1]/(2*sqrt(d)), 0.1)

reach = 1;
indices = [i for i in -reach:δ/10:reach]

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

writedlm("QPDNVexplicit1.csv", samples, ',')