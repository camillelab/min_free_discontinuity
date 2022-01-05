# The bilinear operator E and two linear operators A, B such that
# E(v,xi,mu,sigma,m,p) = <A(v,xi,mu),(sigma,m,p)> = <(v,xi,mu),B(sigma,m,p)>.

# Operator norm of A
const L_0 = 5
#const L_0 = sqrt(12)
#const L_0 = sqrt(8)

# v = Array(N,M), xi = Array(N,M,M+1), eta = Array(N,M,M+1), mu = Array(N,M),
# sigma = Array(N,M,2), m = Array(N,M,M+1), p = Array(N,M+1).

function E_operator(v,xi,eta,mu,sigma,m,p)
    E = 0
    for i = 1:N
        for k = 1:M
            if i < N
                E = E + (v[i+1,k] - v[i,k]) * sigma[i,k,1]
            end
            if k < M
                E = E + (v[i,k+1] - v[i,k]) * sigma[i,k,2]
            end
            E = E + mu[i,k] * (p[i,k+1] - p[i,k] - sigma[i,k,1])
        end
        for k = 1:M
            for l = k+2:M+1
                E = E + xi[i,k,l] * (p[i,l] - p[i,k])
                E = E + eta[i,k,l] * m[i,k,l]
            end
        end
    end
    return E
end

function A_operator(v,xi,eta,mu)
    sigma = zeros(N,M,2)
    m = zeros(N,M,M+1)
    p = zeros(N,M+1)
    for i = 1:N
        for k = 1:M
            if i < N
                sigma[i,k,1] = v[i+1,k] - v[i,k]
            end
            if k < M
                sigma[i,k,2] = v[i,k+1] - v[i,k]
            end
            sigma[i,k,1] = sigma[i,k,1] - mu[i,k]
        end
        for k = 1:M+1
            for l = 1:M+1
                if l <= k-2
                    p[i,k] = p[i,k] + xi[i,l,k]
                elseif l >= k+2
                    p[i,k] = p[i,k] - xi[i,k,l]
                    m[i,k,l] = eta[i,k,l]
                end
            end
            if k == 1
                p[i,k] = p[i,k] - mu[i,k]
            elseif k == M+1
                p[i,k] = p[i,k] + mu[i,k-1]
            else
                p[i,k] = p[i,k] + mu[i,k-1] - mu[i,k]
            end
        end
    end
    return sigma, m, p
end

function B_operator(sigma,m,p)
    v = zeros(N,M)
    xi = zeros(N,M,M+1)
    eta = zeros(N,M,M+1)
    mu = zeros(N,M)
    for i = 1:N
        for k = 1:M
            if i == 1
                v[i,k] = -sigma[i,k,1]
            elseif i == N
                v[i,k] = sigma[i-1,k,1]
            else
                v[i,k] = sigma[i-1,k,1] - sigma[i,k,1]
            end
            if k == 1
                v[i,k] = v[i,k] - sigma[i,k,2]
            elseif k == M
                v[i,k] = v[i,k] + sigma[i,k-1,2]
            else
                v[i,k] = v[i,k] + sigma[i,k-1,2] - sigma[i,k,2]
            end
            for l = k+2:M+1
                xi[i,k,l] = xi[i,k,l] + (p[i,l] - p[i,k])
                eta[i,k,l] = m[i,k,l]
            end
            mu[i,k] = p[i,k+1] - p[i,k] - sigma[i,k,1]
        end
    end
    return v, xi, eta, mu
end
