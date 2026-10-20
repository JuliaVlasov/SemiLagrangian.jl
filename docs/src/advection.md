# Advection

`advection!` performs one in-place, directional semi-Lagrangian step. It
interpolates along one axis of an array; callers control the time loop and
compose directional steps when they need operator splitting.

```julia
using SemiLagrangian

mesh_x = UniformMesh(0.0, 1.0, 128)
mesh_v = UniformMesh(-6.0, 6.0, 96)
interp_x = Lagrange(5)
interp_v = Lagrange(5)
f = [exp(-v^2 / 2) * (1 + 0.01cos(2pi * x)) for x in mesh_x.points, v in mesh_v.points]

# The x velocity varies by velocity line; advect along array dimension 1.
vx = mesh_v.points
advection!(f, mesh_x, vx, 0.01, interp_x; dim = 1)

# The velocity is constant along each velocity line and varies with x.
electric_field = sin.(2pi .* mesh_x.points)
advection!(f, mesh_v, electric_field, 0.01, interp_v; dim = 2)
```

The departure-point shift in grid units is `-dt * velocity / step(mesh)`.
Velocity may be a scalar, an array with one value per transverse line, or an
array matching the data shape for pointwise velocities (the latter requires
an `AbstractInterpolation`). The mesh length must match the selected data
dimension. Interpolation edge behavior is determined by the interpolant.

Compose a Strang-split step explicitly by applying the x direction for half a
step, the v direction for a full step, then x for half a step. More advanced
splitting coefficients are available from the `strangsplit`, `triplejumpsplit`,
and related helpers.

```@autodocs
Modules = [SemiLagrangian]
Order   = [:type, :function]
Pages   = ["splitting.jl", "advection.jl"]
```
