using ITensors
using LinearAlgebra

function rand_SU(n::Int64)
    A = randn(ComplexF64, n, n) + im * randn(ComplexF64, n, n)
    Q, R = qr(A)

    # Make the matrix special unitary (determinant = 1)
    λ = det(Q)
    U = Q / (λ^(1/n))  # Adjust each column to make det(U) = 1

    return U
end

# Define the number of sites
N = 2

# Create a site set for spin-1/2 sites
sites = siteinds("S=1/2", N)

# Initialize a random MPS
psi = randomMPS(sites)

# Define the Ising model Hamiltonian parameters
J = 1.0  # Coupling constant
h = 1.5  # Magnetic field strength
local op
# Construct the Hamiltonian
op = OpSum()
for j = 1:N-1
    op += J, "Sz", j, "Sz", j+1
end
for j = 1:N
    op += h, "Sx", j
end
H = MPO(op, sites)

# Use DMRG to find the ground state
energy, psi_ground = dmrg(H, psi; nsweeps=5, maxdim=10)

# Output the ground state energy
println("Ground state energy: $energy")

anc_site = siteind("S=1/2","Aux")
anc_state = ITensor(anc_site)
anc_state[anc_site => "Up"] = 1.0

# include an ancilla qubit
site_amp = unioninds(sites, anc_site)
psi_amp = MPS(site_amp)
for i = 1:N
    psi_amp[i] = psi_ground[i]
end
psi_amp[N+1] = anc_state

# consrtuct a unitary acting as an ITensor
A = rand_SU(4)
U = ITensor(A, sites[2]', anc_site', sites[2], anc_site)
apply_su!(psi_amp, U, (2, 3))

# pass the total MPS to the 


# Example: Print the MPS tensors
for i = 1:length(psi_ground)
    println(psi_ground[i])
end


