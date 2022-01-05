# Orthogonal projections

# Orthogonal projection onto C
# v = Array(N,N,M)
function projection_C!(v)
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
    return v
end

# Orthogonal projection onto the cylinder
#     sqrt(y1^2 + y2^2) <= R,
# where (y1,y2) are orthonormal coordinates.
function projection_cylinder(x1,x2,R)
    norm = sqrt(x1^2 + x2^2)
    if norm > R
        x1 = R * x1 / norm
        x2 = R * x2 / norm
    end
    return x1, x2
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

# Orthogonal projection onto K
# sigma = Array(N,N,M,M)
function projection_K!(sigma)
    a = 1 / (4 * alpha)

    I = zeros(N,N,M,3)
    J = zeros(N,N,M,M,M,2)

    error = 1
    k = 0
    ## Dykstra algorithm
    while error > E_dykstra && k < K_dykstra
    k = k + 1
    # parabolic constraints
    for i = 1:N
        for j = 1:N
            for k = 1:M
                x1, x2, x3 = sigma[i,j,k,1], sigma[i,j,k,2], sigma[i,j,k,3]
                z1, z2, z3 = sigma[i,j,k,1] - I[i,j,k,1], sigma[i,j,k,2] - I[i,j,k,2], sigma[i,j,k,3] - I[i,j,k,3]
                y1, y2, y3 = projection_parabola(z1,z2,z3,a,a,-rho[i,j,k])
                I[i,j,k,1], I[i,j,k,2], I[i,j,k,3] = I[i,j,k,1] + y1 - x1, I[i,j,k,2] + y2 - x2, I[i,j,k,3] + y3 - x3
                sigma[i,j,k,1], sigma[i,j,k,2], sigma[i,j,k,3] = y1, y2, y3
            end
        end
    end

    # cylindrical constraints
    for i = 1:N
        for j = 1:N
            for k = 1:M
                for l = k+1:M
                    e_norm = sqrt(l-k+1)
                    R = M * psi[i,j,k,l] / e_norm
                    x1, x2 = 0, 0
                    z1, z2 = 0, 0
                    for t = k:l
                        x1, x2 = x1 + sigma[i,j,t,1], x2 + sigma[i,j,t,2]
                        z1, z2 = z1 + sigma[i,j,t,1] - J[i,j,k,l,t,1], z2 + sigma[i,j,t,2] - J[i,j,k,l,t,2]
                    end
                    x1, x2 = x1 / e_norm, x2 / e_norm
                    z1, z2 = z1 / e_norm, z2 / e_norm
                    y1, y2 = projection_cylinder(z1,z2,R)

                    for t = k:l
                        J[i,j,k,l,t,1], J[i,j,k,l,t,2] = J[i,j,k,l,t,1] + (y1 - x1) / e_norm, J[i,j,k,l,t,2] + (y2 - x2) / e_norm
                        sigma[i,j,t,1] = sigma[i,j,t,1] + (y1 - x1) / e_norm
                        sigma[i,j,t,2] = sigma[i,j,t,2] + (y2 - x2) / e_norm
                    end
                end
            end
        end
    end

    error = 0
    for i = 1:N
        for j = 1:N
            for k = 1:M
                x1, x2, x3 = sigma[i,j,k,1], sigma[i,j,k,2], sigma[i,j,k,3]
                error = max(error, - min(x3 - a * x1^2 - a * x2^2 + rho[i,j,k],0))
            end
        end
    end

    for i = 1:N
        for j = 1:N
            for k = 1:M
                for l = k+1:M
                    e_norm = sqrt(l-k+1)
                    R = M * psi[i,j,k,l] / e_norm
                    x1, x2 = 0, 0
                    for t = k:l
                        x1, x2 = x1 + sigma[i,j,t,1], x2 + sigma[i,j,t,2]
                    end
                    x1, x2 = x1 / e_norm, x2 / e_norm
                    error = max(error, - min(R - sqrt(x1^2 + x2^2),0))
                end
            end
        end
    end
    end
    return sigma
end
