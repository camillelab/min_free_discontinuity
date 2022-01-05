# The bilinear operator E and two linear operators A, B such that
# E(v,sigma) = <A(v),sigma> = <v,B(sigma)>.

# Pour les conditions au bord, on fait comme si v était constante en dehors du domaine.

## Operator norm of A
const L_0 = sqrt(12)
#const L_0 = sqrt(8)

# v = Array(N,N,M), sigma = Array(N,N,M,3).
function E_operator(v,sigma)
    E = 0
    for i = 1:N
        for j = 1:N
            for k= 1:M
                if i < N
                    E = E + (v[i+1,j,k] - v[i,j,k]) * sigma[i,j,k,1]
                end
                if j < N
                    E = E + (v[i,j+1,k] - v[i,j,k]) * sigma[i,j,k,2]
                end
                if k < M
                    E = E + (v[i,j,k+1] - v[i,j,k]) * sigma[i,j,k,3]
                end
            end
        end
    end
    return E
end

function A_operator(v)
    sigma = zeros(N,N,M,3)
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
            end
        end
    end
    return sigma
end

function B_operator(sigma)
    v = zeros(N,N,M)
    for i = 1:N
        for j = 1:N
            for k = 1:M
                if i == 1
                    v[i,j,k] = - sigma[i,j,k,1]
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
            end
        end
    end
    return v
end
