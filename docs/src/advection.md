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
dimension. Interpolation boundary behavior is determined by its boundary
condition.

For repeated two-operator compositions, describe the directional operators
once and let `split_advection!` alternate them according to a splitting
coefficient vector. The first operator is applied at odd stages; the second at
even stages. This works with `strangsplit`, `triplejumpsplit`, and the other
coefficient helpers.

```julia
op_x = AdvectionOperator(mesh_x, vx, interp_x; dim = 1)
op_v = AdvectionOperator(mesh_v, electric_field, interp_v; dim = 2)
split_advection!(f, op_x, op_v, strangsplit(0.01))
```

```@autodocs
Modules = [SemiLagrangian]
Private = false
Order   = [:type, :function]
Pages   = ["splitting.jl", "advection.jl"]
```
