# Primal-Dual algorithm

include("operators.jl")
include("projections.jl")

# Steps sigma_0 tau_0 such that sigma_0 tau_0 * L_0^2 < 1,
# where L_0 is the operator norm of the discretized gradient.
const sigma_0 = 1/L_0
const tau_0 = 1/L_0

function primal_dual()
    v = zeros(N,M)
    xi = zeros(N,M,M+1)
    eta = zeros(N,M,M+1)
    mu = zeros(N,M)
    projection_C!(v,xi,eta)

    v_next = zeros(N,M)
    xi_next = zeros(N,M,M+1)
    eta_next = zeros(N,M,M+1)
    mu_next = zeros(N,M)

    v_bar = copy(v)
    xi_bar = copy(xi)
    eta_bar = copy(eta)
    mu_bar = copy(mu)

    sigma = zeros(N,M,2) 
    m = zeros(N,M,M+1) 
    p = zeros(N,M+1)
    for i = 1:N
        for k = 1:M
            for l = k+2:M+1
                m[i,k,l] = psi[i,k,l-1] * M
            end
        end
    end

    k = 0
    error = 1
    energy = E_operator(v,xi,eta,mu,sigma,m,p)

    while error > E_algo && k < K_algo
        A_sigma, A_m, A_p = A_operator(v_bar,xi_bar,eta_bar,mu_bar)
        sigma = projection_K!(sigma + sigma_0 * A_sigma)
        m = m + sigma_0 * A_m
        for i = 1:N
            for k = 1:M
                for l = k+2:M+1
                    if m[i,k,l] > psi[i,k,l-1] * M
                        m[i,k,l] = psi[i,k,l-1] * M
                    end
                end
            end
        end
        p = p + sigma_0 * A_p

        B_v, B_xi, B_eta, B_mu = B_operator(sigma,m,p)
        v_next, xi_next, eta_next = projection_C!(v - tau_0 * B_v, xi - tau_0 * B_xi, eta - tau_0 * B_eta)
        mu_next =  mu - tau_0 * B_mu

        v_bar = 2 * v_next - v
        xi_bar = 2 * xi_next - xi
        eta_bar = 2 * eta_next - eta
        mu_bar = 2 * mu_next - mu

        energy_next = E_operator(v_next,xi_next,eta_next,mu_next,sigma,m,p)
        error = abs(energy_next - energy)
        k = k+1

        energy = energy_next
        v = copy(v_next)
        xi = copy(xi_next)
        eta = copy(eta_next)
        mu = copy(mu_next)
    end 
    return v, sigma
end
