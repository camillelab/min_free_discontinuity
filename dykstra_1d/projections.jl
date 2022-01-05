# Orthogonal projections

# Orthogonal projection onto C
# v = Array(N,M)
function projection_C!(v)
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
    return v
end

# Orthogonal projection onto the cylinder
#     abs(y) <= R.
function projection_cylinder(x,R)
    norm = abs(x)
    if norm > R
        x = R * x / norm
    end
    return x
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
# sigma = Array(N,M,M)
function projection_K!(sigma)
    a = 1 / (4 * alpha)

    I = zeros(N,M,2)
    J = zeros(N,M,M,M)

    error = 1
    k = 0

    # Dykstra algorithm
    while error > E_dykstra && k < K_dykstra
    k = k + 1
    # Orthogonal projections onto the parabolas
    for i = 1:N
        for k = 1:M
            x1, x2 = sigma[i,k,1], sigma[i,k,2]
            z1, z2 = sigma[i,k,1] - I[i,k,1], sigma[i,k,2] - I[i,k,2]
            y1, y2 = projection_parabola(z1,z2,a,-rho[i,k])
            I[i,k,1], I[i,k,2] = I[i,k,1] + y1 - x1, I[i,k,2] + y2 - x2
            sigma[i,k,1], sigma[i,k,2] = y1, y2
        end
    end

    # Orthogonal projections onto the cylinders
    for i = 1:N
        for k = 1:M
            for l = k+1:M
                e_norm = sqrt(l-k+1)
                R = M * psi[i,k,l] / e_norm
                x = 0
                z = 0
                for t = k:l
                    x = x + sigma[i,t,1]
                    z = z + sigma[i,t,1] - J[i,k,l,t]
                end
                x = x / e_norm
                z = z / e_norm
                y = projection_cylinder(z,R)

                for t = k:l
                    J[i,k,l,t] = J[i,k,l,t] + (y - x) / e_norm
                    sigma[i,t,1] = sigma[i,t,1] + (y - x) / e_norm
                end
            end
        end
    end

    # Computation of the error
    error = 0
    for i = 1:N
        for k = 1:M
            x1, x2 = sigma[i,k,1], sigma[i,k,2]
            error = max(error, - min(x2 - a * x1^2 + rho[i,k],0))
        end
    end

    for i = 1:N
        for k = 1:M
            for l = k+1:M
                e_norm = sqrt(l-k+1)
                R = M * psi[i,k,l] / e_norm
                x = 0
                for t = k:l
                    x = x + sigma[i,t,1]
                end
                x = x / e_norm
                error = max(error, - min(R - abs(x),0))
            end
        end
    end
    end
    return sigma
end
