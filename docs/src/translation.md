# Translation

Use `advection!` directly for each directional translation. Here `f` is stored
with position on its first axis and velocity on its second axis:

```julia
using SemiLagrangian

mesh_x = UniformMesh(-5.0, 5.0, 128)
mesh_v = UniformMesh(-5.0, 5.0, 96)
interp_x = Lagrange(5)
interp_v = Lagrange(5)
f = rand(length(mesh_x), length(mesh_v))

# Translate in x using the velocity coordinate for each column.
advection!(f, mesh_x, mesh_v.points, 0.01, interp_x; dim = 1)

# Translate in v using one acceleration for each x row.
electric_field = zeros(length(mesh_x))
advection!(f, mesh_v, electric_field, 0.01, interp_v; dim = 2)
```

Users choose and compose the directional steps themselves, so coupled models
can update fields between advections without an advection-specific callback
object.
