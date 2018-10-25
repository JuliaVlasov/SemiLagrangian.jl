# Semi-Lagrangian method

Let us consider an abstract scalar advection equation of the form

$ \frac{∂f}{∂t}+ a(x, t) ⋅ ∇f = 0. $

The characteristic curves associated to this equation are the solutions of 
the ordinary differential equations

$ \frac{dX}{dt} = a(X(t), t) $

We shall denote by ``X(t, x, s)`` the unique solution of this equation 
associated to the initial condition ``X(s) = x``.

The classical semi-Lagrangian method is based on a backtracking of 
characteristics. Two steps are needed to update the distribution function 
``f^{n+1}`` at ``t^{n+1}`` from its value ``f^n`` at time ``t^n`` :

1. For each grid point ``x_i`` compute ``X(t^n; x_i, t^{n+1})`` the value 
   of the characteristic at ``t^n`` which takes the value ``x_i`` at 
   ``t^{n+1}``.
2. As the distribution solution of first equation verifies
   ``f^{n+1}(x_i) = f^n(X(t^n; x_i, t^{n+1})),``
   we obtain the desired value of ``f^{n+1}(x_i)`` by computing 
   ``f^n(X(t^n;x_i,t^{n+1})`` by interpolation as ``X(t^n; x_i, t^{n+1})`` 
   is in general not a grid point.

*[Eric Sonnendrücker - Numerical methods for the Vlasov equations](http://www-m16.ma.tum.de/foswiki/pub/M16/Allgemeines/NumMethVlasov/Num-Meth-Vlasov-Notes.pdf)*

```@meta
CurrentModule = Splittings
```

## Advection functions

```@docs
advection!(::Array{Float64,2}, ::UniformMesh, v, ::CubicSpline, ::Float64, ::Int64)
```

```@docs
advection!(::Array{Float64,2}, ::UniformMesh, v, ::BSpline, ::Float64, ::Int64)
```
