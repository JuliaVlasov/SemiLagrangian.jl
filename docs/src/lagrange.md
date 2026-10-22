# Lagrange interpolation


To compute a Lagrange interpolation you have to instantiate
a `Lagrange` type variable:

`Lagrange` and `Hermite` accept a `bc` keyword. The default,
`PeriodicBC()`, wraps interpolation stencils around the ends of the data line.
Use `ExtrapolationBC()` to keep stencils within the line and extrapolate near
its ends:

```julia
using SemiLagrangian

periodic_interp = Lagrange(5; bc = PeriodicBC())
extrapolating_interp = Lagrange(5; bc = ExtrapolationBC())
```

```@autodocs
Modules = [SemiLagrangian]
Private = false
Order   = [:type, :function]
Pages   = ["lagrange.jl", "lagrange_interp.jl"]
```
