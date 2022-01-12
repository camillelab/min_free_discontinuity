using JLD
include("primal_dual.jl")

# We minimize
# E(u) = \alpha \int \abs{\nabla u}^2 + \int_{J_u} psi(x,u^-,u^+) + \int \rho(x,u)
# where $u:[0,1]^2 \to \R$ is a piecewise smooth function whose jump set is $J_u$.

# The algorithm is the primal-dual algorithm designed by
# Bischof, Chambolle, Cremers and Pock in the article
# "An Algorithm for Minimizing the Mumford-Shah Functional".
# However, the method also extends to other free-discontinuity functionals,
# provided one adapts the definition of $\alpha$, $\psi$ and $\rho$.

## Discretization of the domain and the color channel
const N, M = 16, 16
#const N, M = 32, 32

## Precision constants in the Dykstra projection algorithm
const E_dykstra = 0.001
const K_dykstra = 25
## Precision constants in the primal-dual algorithm
const E_algo = 0.0001
const K_algo = 5000

## Dirichlet conditions
# We define an indicator function S:[0,1]^2 \to {0,1}$ and we define
# a function U0:[0,1]^2 \to [0,1] such that for all coordinates (i,j),
# the condition S(i,j) = 1 imposes u(i,j) = U0(i,j).
const S = Array{Int8}(undef,N,N)
const U0 = Array{Float64}(undef,N,N)
for i = 1:N
    for j = 1:N
        if i == 1 || i == N || j == 1 || j == N
            S[i,j] = 1
            U0[i,j] = 0
        else
            S[i,j] = 0
        end
    end
end

## Definition of the energy (gamma = 0 for the homogeneous functional)
# Set gamma = 0 for the homogeneous functional
const alpha, beta, gamma = 1, 1/4, 0
# Set the surface energy
const psi = Array{Float64}(undef,N,N,M,M)
for i = 1:N
    for j = 1:N
        for k = 1:M
            for l = 1:M
                # Mumford-Shah functional
                psi[i,j,k,l] = beta
                # Thermal Insulation functional
                # psi[i,j,k,l] = beta * ((k/M)^2 + (l/M)^2)
            end
        end
    end
end
# Set the error term
# For the Mumford-Shah functional, we need to define an image f first.
function f(x,y)
    if sqrt((x-1/2)^2 + (y-1/2)^2) < 1/8
        return 1
    else
        return 0
    end
end
const rho = Array{Float64}(undef,N,N,M)
for i = 1:N
    for j = 1:N
        for k = 1:M
            # Mumford-Shah functional
            rho[i,j,k] = gamma * (k/M - f(i/N))^2 
            # Thermal Insulation functional
            #if k == 1
            #    rho[i,k] = 0
            #else
            #    rho[i,k] = gamma^2
            #end
        end
    end
end

const v, sigma = primal_dual()
save("data.jld", "v", v, "sigma", sigma, "N", N, "M", M)
