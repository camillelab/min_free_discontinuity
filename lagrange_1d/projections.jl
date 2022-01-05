# Orthogonal projections

# Orthogonal projection onto C
# v = Array(N,M), xi = Array(N,M,M+1), eta = Array(N,M,M+1)
function projection_C!(v,xi,eta)
    for i = 1:N
        v[i,1] = 1
        v[i,M] = 0
        if i == 1 # Dirichlet condition u = m1 on i = 1
            for k = 2:M-1
                if k/M < m1
                    v[i,k] = 1
                else
                    v[i,k] = 0
                end
            end
        elseif i == N # Dirichlet condition u = m2 on i = N
            for k = 2:M-1
                if k/M < m2
                    v[i,k] = 1
                else
                    v[i,k] = 0
                end
            end
        else
            for k = 2:M-1
                v[i,k] = min(1,max(0,v[i,k]))
            end
        end
    end

    # Orthogonal projection onto the set {norm(xi) <= eta}
    for i = 1:N
        for k = 1:M
            for l = 1:M+1
                xi_norm = abs(xi[i,k,l])
                if eta[i,k,l] < xi_norm
                    if eta[i,k,l] <= - xi_norm
                        eta[i,k,l] = 0
                        xi[i,k,l] = 0
                    else
                        eta[i,k,l] = 1/2 * (eta[i,k,l] + xi_norm)
                        xi[i,k,l] = eta[i,k,l] * xi[i,k,l] / xi_norm
                    end
                end
            end
        end
    end
    return v, xi, eta
end

# Orthogonal projection onto 
#     {(y1,y2) | y2 >= a1 y1^2 + c },
# where (y1,y2) are orthonormal coordinates.
# The idea is that if y is the orthogonal projection of x,
# then there exists a scalar t such that
#     (y1 - x1, y2 - x2) = t (2 a1 y1, -1)
# and since y2 = a1 y1^2 + c, we deduce an equation f(t) = 0.
const E_parabola = 0.000001
const K_parabola = 15
function projection_parabola(x1,x2,a1,c)
    if x2 < a1 * x1^2 + c - E_parabola
        k = 0
        t = 0
        ft = a1 * x1^2 + c - x2
        while abs(ft) > E_parabola && k < K_parabola
            k = k + 1
            ft = t + a1 * (x1 / (1 - 2 * a1 * t))^2 + c - x2 
            dft = 1 + (4 * a1^2 * x1^2) / (1 - 2 * a1 * t)^3
            t = t - ft/dft
        end
        x1 = x1 / (1 - 2 * a1 * t)
        x2 = x2 - t
    end
    return x1, x2
end

# Orthogonal projection onto K
# sigma = Array(N,M,2)
function projection_K!(sigma)
    a = 1 / (4 * alpha)

    for i = 1:N
        for k = 1:M
            sigma[i,k,1], sigma[i,k,2] = projection_parabola(sigma[i,k,1],sigma[i,k,2],a,-rho[i,k])
        end
    end
    return sigma
end
