## Meshes

`UniformMesh` represents a one-dimensional, evenly spaced grid over a
half-open interval `[a, b)`. Construct one with the interval boundaries and
the number of cells:

```julia
using SemiLagrangian

mesh = UniformMesh(-1.0, 1.0, 100)
```

It contains 100 points from `-1.0` up to, but not including, `1.0`. The
spacing is `(b - a) / 100`; equivalently, the mesh width is the number of
points times the spacing. This convention is useful for periodic grids,
where storing the right endpoint would duplicate the left endpoint. The
constructor uses the same point convention for non-periodic interpolation;
the interpolation boundary condition determines how data near the ends is
handled.

Use `points(mesh)` to access coordinates, `step(mesh)` for grid spacing,
`start(mesh)` and `stop(mesh)` for the interval boundaries, and `length(mesh)`
for the point count. `stdtomesh` and `meshtostd` convert between physical
coordinates and grid units. `vec_k_fft` returns Fourier wavenumbers scaled to
the mesh width.

The `traitmodbegin!` and `traitmodend!` helpers adjust periodic-valued arrays
in place: the former unwraps neighboring values to avoid jumps across a
periodic boundary, while the latter wraps values into the mesh interval.

```@autodocs
Modules = [SemiLagrangian]
Private = false
Order   = [:type, :function]
Pages   = ["mesh.jl"]
```
