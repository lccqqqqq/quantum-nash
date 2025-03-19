using HomotopyContinuation
using LinearAlgebra, IterTools
using CSV, DelimitedFiles

include("list_of_functions.jl")


# The payoff functions 
h1 = [3 0 0 0; 0 0 0 0; 0 0 5 0; 0 0 0 1];
h2 = [3 0 0 0; 0 5 0 0; 0 0 0 0; 0 0 0 1];
d = 5;
@var x[1:5]

ψ = [x[1], x[2], x[3], x[4]+1im*x[5]];
ψT = transpose([x[1], x[2], x[3], x[4]-1im*x[5]]);
cstr = [sum(x.*x)-1,
    ψT * 1im*com(h1, pauli(2, 1, 2)) * ψ,
    ψT * 1im*com(h2, pauli(2, 2, 2)) * ψ,
    ψT * com(h1, pauli(2, 1, 1)) * ψ,
    ψT * com(h2, pauli(2, 2, 1)) * ψ];

D = 5; # Length of Constraints
@var λ[1:D] μ[1:D]; 
∇ = differentiate(cstr, x)
@var y[1:d];
χ = [y[1], y[2], y[3], y[4]+1im*y[5]];
χT = transpose([y[1], y[2], y[3], y[4]-1im*y[5]]);

replica_sys = subs(cstr, x => y)
∇y = subs(∇, x => y)
system = [cstr; replica_sys; map(j -> x[j]-y[j]-dot(λ, ∇[:, j]), 1:d); map(j -> x[j]-y[j]-dot(μ, ∇y[:, j]), 1:d)]

result = solve(system)

solve(cstr)

reach = 1; δ = 0.1;
indices = [i for i in -reach:δ/10:reach];
