using HomotopyContinuation
using LinearAlgebra, IterTools
using CSV, DelimitedFiles

include("list_of_functions.jl");

h1 = readdlm("nash_payoff_h1XXlinked.csv", ',')
h2 = readdlm("nash_payoff_h2XXlinked.csv", ',')

d = 7;
@var x[1:d]

ψ = [x[1], x[2]+1im*x[3], x[4]+1im*x[5], x[6]+1im*x[7]];
ψT = transpose([x[1], x[2]-1im*x[3], x[4]-1im*x[5], x[6]-1im*x[7]]);

h1 = h1 + 1e-10*rand_herm(4);
h2 = h2 + 1e-10*rand_herm(4);



cstr = [sum(x.*x)-1,
    ψT * 1im*com(h1, pauli(2, 1, 2)) * ψ,
    ψT * 1im*com(h2, pauli(2, 2, 2)) * ψ,
    ψT * 1im*com(h1, pauli(2, 1, 1)) * ψ,
    ψT * 1im*com(h2, pauli(2, 2, 1)) * ψ,
    ψT * 1im*com(h1, pauli(2, 1, 3)) * ψ,
    ψT * 1im*com(h2, pauli(2, 2, 3)) * ψ];

result = solve(cstr)
realsols = real_solutions(nonsingular(result))

writedlm("realsols3.csv", realsols, ',')
disp(subs(cstr, x => realsols[end]))