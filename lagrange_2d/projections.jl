# Orthogonal projection onto C
# Lagrange multipliers !!

## Two-dimensional case

# v = Array(N,N,M), xi = Array(N,N,M,M+1,2), eta = Array(N,N,M,M+1)

function projection_C!(v,xi,eta)
    for i = 1:N
        for j = 1:N
            v[i,j,1] = 1
            v[i,j,M] = 0
            if S[i,j] == 1 # Dirichlet condition (u = U0 on S)
                for k = 2:M-1
                    if k/M < U0[i,j]
                        v[i,j,k] = 1
                    else
                        v[i,j,k] = 0
                    end
                end
            else
                for k = 2:M-1
                    v[i,j,k] = min(1,max(0,v[i,j,k]))
                end
            end
        end
    end

    # Orthogonal projection onto the set { norm(xi) <= eta }
    for i = 1:N
        for j = 1:N
            for k = 1:M
                for l = 1:M+1
                    xi_norm = sqrt(xi[i,j,k,l,1]^2 + xi[i,j,k,l,2]^2)
                    if eta[i,j,k,l] < xi_norm
                        if eta[i,j,k,l] <= - xi_norm
                            eta[i,j,k,l] = 0
                            xi[i,j,k,l,1] = 0
                            xi[i,j,k,l,2] = 0
                        else
                            eta[i,j,k,l] = 1/2 * (eta[i,j,k,l] + xi_norm)
                            xi[i,j,k,l,1] = eta[i,j,k,l] * xi[i,j,k,l,1] / xi_norm
                            xi[i,j,k,l,2] = eta[i,j,k,l] * xi[i,j,k,l,2] / xi_norm
                        end
                    end
                end
            end
        end
    end
    return v, xi, eta
end

# Orthogonal projection onto
#     {(y1,y2,y3) | y3 >= a1 y1^2 + a2 y2^2 + c },
# where (y1,y2,y3) are orthonormal coordinates.
# The idea is that if y is the orthogonal projection of x,
# then there exists a scalar t such that
#     (y1 - x1, y2 - x2, y3 - x3) = t (2 a1 y1, 2a2 y2, -1)
# and since y3 = a1 y1^2 + a2 y2^2 + c, we deduce an equation f(t) = 0.
const E_parabola = 0.000001
const K_parabola = 15
function projection_parabola(x1,x2,x3,a1,a2,c)
    if x3 < a1 * x1^2 + a2 * x2^2 + c - E_parabola
        k = 0
        t = 0
        ft = a1 * (x1)^2 + a2 * x2^2 + c - x3 
        while abs(ft) > E_parabola && k < K_parabola
            k = k + 1
            ft = t + a1 * (x1 / (1 - 2 * a1 * t))^2 + a2 * (x2 / (1 - 2 * a2 * t))^2 + c - x3 
            dft = 1 + (4 * a1^2 * x1^2) / (1 - 2 * a1 * t)^3 + (4 * a2^2 * x2^2) / (1 - 2 * a2 * t)^3
            t = t - ft/dft
        end
        x1 = x1 / (1 - 2 * a1 * t)
        x2 = x2 / (1 - 2 * a2 * t)
        x3 = x3 - t
    end
    return x1, x2, x3
end

function projection_K!(sigma)
    a = 1 / (4 * alpha)

    for i = 1:N
        for j = 1:N
            for k = 1:M
                sigma[i,j,k,1], sigma[i,j,k,2], sigma[i,j,k,3] = projection_parabola(sigma[i,j,k,1],sigma[i,j,k,2],sigma[i,j,k,3],a,a,-rho[i,j,k])
            end
        end
    end
    return sigma
end
