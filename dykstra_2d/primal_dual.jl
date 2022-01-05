# Primal-Dual algorithm

include("operators.jl")
include("projections.jl")

# Steps sigma_0 tau_0 such that sigma_0 tau_0 * L_0^2 < 1
# where L_0 is the operator norm of the discretized gradient.
const sigma_0 = 1/L_0
const tau_0 = 1/L_0

function primal_dual()
    v = zeros(N,N,M)
    projection_C!(v)

    v_next = zeros(N,N,M)
    v_bar = copy(v)

    sigma = zeros(N,N,M,3) 

    k = 0
    error = 1
    energy = E_operator(v,sigma)

    while error > E_algo && k < K_algo
        A_sigma = A_operator(v_bar)
        sigma = projection_K!(sigma + sigma_0 * A_sigma)

        B_v = B_operator(sigma)
        v_next = projection_C!(v - tau_0 * B_v)

        v_bar = 2 * v_next - v

        energy_next = E_operator(v_next,sigma)
        error = abs(energy_next - energy)
        k = k + 1

        energy = energy_next
        v = copy(v_next)
    end 
    return v, sigma
end
