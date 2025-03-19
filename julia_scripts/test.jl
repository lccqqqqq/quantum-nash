using HomotopyContinuation, DynamicPolynomials, LinearAlgebra, IterTools

include("list_of_functions.jl")


# Using witness set to determine the dimension of Nash manifold for 3 person games
N = 3;
n = 2^N;
H = [rand_herm(n) for _ in 1:N]
G = [map(i -> [1im * com(H[i], pauli(N, i, a)) for a in 1:3], 1:N)][1]
@var x[1:2*n]
cstr = [x' * x - 1; x[n+1]]



# new model: Ising
for i = 1:N
    for a = 1:3
        global cstr = vcat(cstr, (x[1:n]'-1im*x[n+1:2*n]') * (G[i][a] * (x[1:n]+1im*x[n+1:2*n])))
    end
end

L = rand_subspace(2*n; dim = length(cstr), affine=true, real=true)
W = witness_set(System(cstr), L)

# we can get some samples for out manifold via a predictor-corrector procedure by linearizing around the obtained soliution.
W_real = real_solutions(results(W))

# linearization
@var χ[1:2*n]



# Can only get a coordinate patch and analyze the patch when slightly perturbing the game payoffs
Dim = 5
sol = W_real[1][1:n] + 1im * W_real[1][n+1:2*n]
iter_range = 1:2 # so that there will be 3^5=243 sample points
δx = 0.1
iters = IterTools.product(ntuple(_ -> iter_range, Dim)...)

sample_points = []

struct CoordinatePatchData
    index::Tuple{Vararg{Int}}
    x::Vector{ComplexF64}
    tangent_space::Matrix{ComplexF64}
end

# result from dimension counting, insert for general formula

hamming_distance(x::Tuple, y::Tuple) = sum([abs(x[i]-y[i]) for i in 1:length(x)])
hamming(x::Tuple, y::Tuple) = [abs(x[i]-y[i]) for i in 1:length(x)]
for iter in iters
    if sum(iter) == Dim # meaning we are at the starting points
        sol_data = CoordinatePatchData(iter, sol, A)
        push!(sample_points, sol_data)
    else
        # choose the previous point with Hamming distance 1
        i = 1
        rev_sample_points = reverse(sample_points)
        while i <= length(sample_points)
            index_prev = rev_sample_points[i].index
            if hamming_distance(index_prev, iter) == 1
            # when we find a nearest neighbour whose information is known
            # we would like to find the move direction
            previous_point = rev_sample_points[i]
                for (move_direction, val) in enumerate(hamming(index_prev, iter))
                    if val == 1
                        break
                    end
                end
            break
        end
    i += 1
    end

    # use Newton's method to calculate the new starting point
    sol = previous_point.x
    A = previous_point.tangent_space

    newton_start_point = normalize!(sol + δx * (A[1:n, move_direction] + 1im*A[n+1:2*n, move_direction]))
    # generate a new point via newton() in the Homotopy continuation package
    newton_result = newton(System(cstr), [real(newton_start_point);imag(newton_start_point)], atol=1e-10)

    new_sample = newton_result.x[1:n] + 1im * newton_result.x[n+1:2*n]

    # need information on the new sample:
    # index is current iter
    # sample point stored in newton_result
    # consider the tangent space: need linearization around the new point

    ### The linearized conditions:

    # Compute all linearized constraints stem from the Nash constraints
    linearized_cstr = [sol'*(G[i][a]*(χ[1:n]+1im*χ[n+1:2*n])) + (χ[1:n]'-1im*χ[n+1:2*n]')*(G[i][a]*sol) for i in 1:length(G) for a in 1:3]
    # Implement normalization
    linearized_cstr = vcat(linearized_cstr, [sol'*(χ[1:n]+1im*χ[n+1:2*n]) + (χ[1:n]'-1im*χ[n+1:2*n]')*sol])
    # extract all coefficients
    linearized_cstr_coefs = [coefficients(linearized_cstr[i], χ; expanded=false) for i = 1:length(linearized_cstr)]
    # The coefficients for the normalization constraints has one empty entry
    if size(linearized_cstr_coefs[end], 1) < size(linearized_cstr_coefs[1], 1)
        nmlz = copy(linearized_cstr_coefs[end])
        insert!(nmlz, n, 0)
        linearized_cstr_coefs[end] = nmlz
    end
    linearized_cstr_coefs = mapreduce(permutedims, vcat, linearized_cstr_coefs) # the matrix of coefficients
    linearized_cstr_coefs = vcat(linearized_cstr_coefs, [(x==n+1) ? 1 : 0 for x in 1:2*n]') 
    # the obtained set of basis vectors on the tangent space of the manifold
    A_new = nullspace(linearized_cstr_coefs)

    # projection of previous point's basis vectors onto the new embedded subspace.
    # Ideally the basis vectors should be made consistent via parallel transport, but for now we are moving only a small 
    # distance around, so projection should do the job with reasonable accuracy.
    new_tangent_space = zeros(2*n, Dim)
    new_basis_vec_coefs = transpose(previous_point.tangent_space) * A_new # a Dim-by-Dim matrix of coefficients
    for j = 1:Dim
        new_tangent_space[:, j] = sum(map(k -> new_basis_vec_coefs[j, k] * A_new[:, k], 1:Dim))
    end

    sol_data = CoordinatePatchData(iter, new_sample, new_tangent_space)
    end
    println(iter)
end