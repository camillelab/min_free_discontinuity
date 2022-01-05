# README

## Introduction
We minimize
```
    E(u) = \alpha \int \abs{\nabla u}^2 + \int_{J_u} psi(x,u^-,u^+) + \int \rho(x,u)
```
where `u:[0,1]^d \to [0,1]` is a piecewise smooth function whose jump set is `J_u` and `d = 1, 2`.

In the folders `dykstra/`, we implement primal-dual algorithm designed by Bischof, Chambolle, Cremers and Pock in the article _An Algorithm for Minimizing the Mumford-Shah Functional_.
In the folders `lagrange/`, we implement the primal-dual algorithm designed by Chambolle, Cremers and Strekalovskiy in the article _A Convex Representation for the Vectorial Mumford-Shah Functional_.

We underline that the method extends easily to other free-discontinuity functionals than Mumford-Shah by adapting the definition of `\alpha`, `\psi` and `\rho`.
In particular, we have written the definition of `\psi` and `\rho` for both the Mumford-Shah functional and the Thermal Insulation Functional.

It is not clear which method produces the best result in the least time.
The `lagrange` method has fastest iteratations but requires more iterations to produce a relevant output.
In 2D, the computations are very demanding so it is is better to use an HPC computation system. There are also better alternative to Julia.

## How to run the program

Launch `julia main.jl`.
The program outputs a file `data.jld` containing a Julia dictionary with 4 keys.
The keys "v" and "sigma" are an optimal pair (v,sigma) of the dual problem.
The keys "N" and "M" are the numbers used to discretize the domain and the color channel.
More precisely, the functions `v: [0,1]^d \times [0,1] \to \R` and `\sigma: [0,1]^d \times [0,1] \to \R^3` are discretized as functions `v:{1,...,N}^d \times {1,...,M} \to [0,1]` and `\sigma: {1,...,N}^d \times {1,...,M} \to \R^3`.

The corresponding minimizer `u` can be defined as the 0.5-isosurface of `v`.
Hence, one can plot `u` in a Jupyter session with the following code.

```
    using Plots
    using JLD

    data = load("./data.jld")
    v = data["v"]
    sigma = data["sigma"]
    N = data["N"]
    M = data["M"]

    I = Array{Float64}(undef,N)
    J = Array{Float64}(undef,N)
    for i = 1:N
        I[i] = i/N
    end
    for k = 1:M
        J[k] = k/M
    end

    ## 1D
    u = zeros(N)
    for i = 1:N
        k = 1
        while k < M && v[i,k] > 0.5
            k = k + 1
        end
        u[i] = k/M
    end
    p = plot(I,u,xlim=(1/M,1),ylim=(1/M,1),zlim=(1/M,1))
    display(p)

    ## 2D
    #u = zeros(N,N);
    #for i = 1:N
    #    for j = 1:N
    #        k = 1
    #        while k < M && v[i,j,k] > 0.5
    #            k = k + 1
    #        end
    #        u[i,j] = k/M
    #    end
    #end
    #p = plot(I,I,u,seriestype=:surface,xlim=(1/M,1),ylim=(1/M,1),zlim=(1/M,1))
    #display(p)
```

## Examples of applications

### 1D Homogeneous Mumford-Shah functional

We minimize
```
    E(u) = \int (u')^2 + \beta \HH^0(J_u)
```
where `u:[0,1] \to [0,1]` is a piecewise smooth function whose jump set is `J_u`.
Say that `u(0) = m` and `u(1) = M`.
There are two possible minimizers.
If
```
(M - m)^2 \leq \beta
```
then the competitor which joins `m` and `M` linearly is a minimizer.
Otherwise, the piewise constant function with one jump is a minimizer.


### 1D Homogenous Thermal insulation functional

We minimize
```
    E(u) = \int (u')^2 + \beta \int_{J_u} (u^-)^2 + (u^+)^2
```
where `u:[0,1] \to [0,1]` is a piecewise smooth function whose jump set is `J_u`.
Say that `u(0) = m` and `u(1) = M`.
There are two possible minimizers.
If
```
    E(u) = \alpha \int \abs{\nabla u}^2 + \int_{J_u} psi(x,u^-,u^+) + \int \rho(x,u)
```
where `\delta = \frac{M}{1 + \beta}`, then the affine competitor is a mimimizer.
Otherwise, the function which jumps from `m` to `\delta` at `x = 0` and which
joins `\delta` and `M` linearly is a minimizer.

### Crack-tip

For the 2D Mumford-Shah functional, one can set the crack-tip as a Dirichlet condition. We recall its expression in polar coordinates
```
u(r,\theta) = \sqrt{\frac{2r}{\pi}} \sin(\tfrac{\theta}{2})
```
where `r > 0`, `\theta \in ]-\pi,\pi[` or in cartesian coordinates
```
u(x,y) = \mathrm{sgn}(y) \sqrt{frac{r-x}{\pi}}.
```
