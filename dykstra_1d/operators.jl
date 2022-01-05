# The bilinear operator E and two linear operators A, B such that
# E(v,sigma) = <A(v),sigma> = <v,B(sigma)>.

# Pour les conditions au bord, on fait comme si v était constante en dehors du domaine.

## Operator norm of A

const L_0 = sqrt(12)
#const L_0 = sqrt(8)

# v = Array(N,M), sigma = Array(N,M,2).

function E_operator(v,sigma)
    F = 0
    for i = 1:N
        for k= 1:M
            if i < N
                F = F + (v[i+1,k] - v[i,k]) * sigma[i,k,1]
            end
            if k < M
                F = F + (v[i,k+1] - v[i,k]) * sigma[i,k,2]
            end
        end
    end
    return F
end

# Operator A such that F = <A(v),sigma)>.
# v = Array(N,M), sigma = Array(N,M,3).
# Using |a+b|^2 <= 2(|a|^2+|b|^2), we get |A(v)|^2 <= 12 |v|^2.
function A_operator(v)
    sigma = zeros(N,M,2)
    for i = 1:N
        for k = 1:M
            if i < N
                sigma[i,k,1] = v[i+1,k] - v[i,k]
            end
            if k < M
                sigma[i,k,2] = v[i,k+1] - v[i,k]
            end
        end
    end
    return sigma
end

# Operator B such that F = <v,B(sigma)>.
# v = Array(N,M), sigma = Array(N,M,3).
function B_operator(sigma)
    v = zeros(N,M)
    for i = 1:N
            for k = 1:M
                if i == 1
                    v[i,k] = - sigma[i,k,1]
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
            end
    end
    return v
end
