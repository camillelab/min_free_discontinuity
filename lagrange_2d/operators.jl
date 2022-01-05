# The bilinear operator E and two linear operators A, B such that
# E(v,xi,mu,sigma,m,p) = <A(v,xi,mu),(sigma,m,p)> = <(v,xi,mu),B(sigma,m,p)>.

## Operator norm of A
const L_0 = 5
#const L_0 = sqrt(12)
#const L_0 = sqrt(8)

# v = Array(N,N,M), xi = Array(N,N,M,M+1,2), eta = Array(N,N,M,M+1), mu = Array(N,N,M,2),
# sigma = Array(N,N,M,3), m = Array(N,N,M,M+1), p = Array(N,N,M+1,2).

function E_operator(v,xi,eta,mu,sigma,m,p)
    E = 0
    for i = 1:N
        for j = 1:N
            for k = 1:M
                if i < N
                    E = E + (v[i+1,j,k] - v[i,j,k]) * sigma[i,j,k,1]
                end
                if j < N
                    E = E + (v[i,j+1,k] - v[i,j,k]) * sigma[i,j,k,2]
                end
                if k < M
                    E = E + (v[i,j,k+1] - v[i,j,k]) * sigma[i,j,k,3]
                end
                E = E + mu[i,j,k,1] * (p[i,j,k+1,1] - p[i,j,k,1] - sigma[i,j,k,1])
                E = E + mu[i,j,k,2] * (p[i,j,k+1,2] - p[i,j,k,2] - sigma[i,j,k,2])
            end
            for k = 1:M
                for l = k+2:M+1
                    E = E + xi[i,j,k,l,1] * (p[i,j,l,1] - p[i,j,k,1])
                    E = E + xi[i,j,k,l,2] * (p[i,j,l,2] - p[i,j,k,2])
                    E = E + eta[i,j,k,l] * m[i,j,k,l]
                end
            end
        end
    end
    return E
end

function A_operator(v,xi,eta,mu)
    sigma = zeros(N,N,M,3)
    m = zeros(N,N,M,M+1)
    p = zeros(N,N,M+1,2)
    for i = 1:N
        for j = 1:N
            for k = 1:M
                if i < N
                    sigma[i,j,k,1] = v[i+1,j,k] - v[i,j,k]
                end
                if j < N
                    sigma[i,j,k,2] = v[i,j+1,k] - v[i,j,k]
                end
                if k < M
                    sigma[i,j,k,3] = v[i,j,k+1] - v[i,j,k]
                end
                    sigma[i,j,k,1] = sigma[i,j,k,1] - mu[i,j,k,1]
                    sigma[i,j,k,2] = sigma[i,j,k,2] - mu[i,j,k,2]
            end
            for k = 1:M+1
                for l = 1:M+1
                    if l <= k-2
                        p[i,j,k,1] = p[i,j,k,1] + xi[i,j,l,k,1]
                        p[i,j,k,2] = p[i,j,k,2] + xi[i,j,l,k,2]
                    elseif l >= k+2
                        p[i,j,k,1] = p[i,j,k,1] - xi[i,j,k,l,1]
                        p[i,j,k,2] = p[i,j,k,2] - xi[i,j,k,l,2]
                        m[i,j,k,l] = eta[i,j,k,l]
                    end
                end
                if k == 1
                    p[i,j,k,1] = p[i,j,k,1] - mu[i,j,k,1]
                    p[i,j,k,2] = p[i,j,k,2] - mu[i,j,k,2]
                elseif k == M+1
                    p[i,j,k,1] = p[i,j,k,1] + mu[i,j,k-1,1]
                    p[i,j,k,2] = p[i,j,k,2] + mu[i,j,k-1,2]
                else
                    p[i,j,k,1] = p[i,j,k,1] + mu[i,j,k-1,1] - mu[i,j,k,1]
                    p[i,j,k,2] = p[i,j,k,2] + mu[i,j,k-1,2] - mu[i,j,k,2]
                end
            end
        end
    end
    return sigma, m, p
end

function B_operator(sigma,m,p)
    v = zeros(N,N,M)
    xi = zeros(N,N,M,M+1,2)
    eta = zeros(N,N,M,M+1)
    mu = zeros(N,N,M,2)
    for i = 1:N
        for j = 1:N
            for k = 1:M
                if i == 1
                    v[i,j,k] = -sigma[i,j,k,1]
                elseif i == N
                    v[i,j,k] = sigma[i-1,j,k,1]
                else
                    v[i,j,k] = sigma[i-1,j,k,1] - sigma[i,j,k,1]
                end
                if j == 1
                    v[i,j,k] = v[i,j,k] - sigma[i,j,k,2]
                elseif j == N
                    v[i,j,k] = v[i,j,k] + sigma[i,j-1,k,2]
                else
                    v[i,j,k] = v[i,j,k] + sigma[i,j-1,k,2] - sigma[i,j,k,2]
                end
                if k == 1
                    v[i,j,k] = v[i,j,k] - sigma[i,j,k,3]
                elseif k == M
                    v[i,j,k] = v[i,j,k] + sigma[i,j,k-1,3]
                else
                    v[i,j,k] = v[i,j,k] + sigma[i,j,k-1,3] - sigma[i,j,k,3]
                end
                for l = k+2:M+1
                    xi[i,j,k,l,1] = xi[i,j,k,l,1] + p[i,j,l,1] - p[i,j,k,1]
                    xi[i,j,k,l,2] = xi[i,j,k,l,2] + p[i,j,l,2] - p[i,j,k,2]
                    eta[i,j,k,l] = m[i,j,k,l]
                end
                mu[i,j,k,1] = p[i,j,k+1,1] - p[i,j,k,1] - sigma[i,j,k,1]
                mu[i,j,k,2] = p[i,j,k+1,2] - p[i,j,k,2] - sigma[i,j,k,2]
            end
        end
    end
    return v, xi, eta, mu
end
