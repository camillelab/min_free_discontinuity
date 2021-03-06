# README

## Introduction
We minimize
```
    E(u) = \alpha \int \abs{\nabla u}^2 + \int_{J_u} psi(x,u^-,u^+) + \int \rho(x,u)
```
where `u:[0,1]^d \to [0,1]` is a piecewise smooth function whose jump set is `J_u` and `d = 1, 2`.

The folders `dykstra_*/` contain an implementation of the primal-dual algorithm designed by Pock, Cremers, Bischof and Chambolle in [_An Algorithm for Minimizing the Mumford-Shah Functional_](https://doi.org/10.1109/iccv.2009.5459348).
The folders `lagrange_*/` contain an alternative method presented by Strekalovskiy, Chambolle and Cremers in [_A Convex Representation for the Vectorial Mumford-Shah Functional_](https://doi.org/10.1109/cvpr.2012.6247866).

We underline that the algorithm extends easily to other free-discontinuity functionals than Mumford-Shah by adapting the definition of `\alpha`, `\psi` and `\rho`.
In particular, the files `main.jl` provides the definition of `\psi` and `\rho` for the Thermal Insulation Functional.

It is not clear which method produces the best result in the least time.
The `lagrange` method has fastest iterations but requires more iterations to produce a relevant output.
In 2D, the computations are very demanding so it is is better to use an HPC computation system. There are also better alternative to Julia.

## How to run the program

Launch `julia main.jl`.
The program outputs a file `data.jld` containing a Julia dictionary with 4 keys.
The keys `"v"` and `"sigma"` are an optimal pair `(v,sigma)` of the dual problem.
The keys `"N"` and `"M"` are the numbers used to discretize the domain and the color channel.
More precisely, the functions
```
v: [0,1]^d \times [0,1] \to \R
\sigma: [0,1]^d \times [0,1] \to \R^3
```
are discretized as functions
```
v:{1,...,N}^d \times {1,...,M} \to [0,1]
\sigma: {1,...,N}^d \times {1,...,M} \to \R^3.
```

The corresponding minimizer `u` can be defined as the 0.5-isosurface of `v`.
For example, one can plot `u` in a Jupyter session with the following code.

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
(M - m)^2 \leq (M - \delta)^2 + \beta(m^2 + \delta^2)
```
where `\delta = \frac{M}{1 + \beta}`, then the affine competitor is a mimimizer.
Otherwise, the function which jumps from `m` to `\delta` at `x = 0` and which
joins `\delta` and `M` linearly is a minimizer.

### Crack-tip

For the 2D homogeneous Mumford-Shah functional, one can set the crack-tip as a Dirichlet condition. We recall its expression in polar coordinates
```
u(r,\theta) = \sqrt{\frac{2r}{\pi}} \sin(\tfrac{\theta}{2})
```
where `(r,\theta) \in ]0,\infty[ \times ]-\pi,\pi[` or in cartesian coordinates
```
u(x,y) = \mathrm{sgn}(y) \sqrt{frac{r-x}{\pi}}.
```
where `(x,y) \in \R^2 \setminus {y = 0, x \leq 0}`.
