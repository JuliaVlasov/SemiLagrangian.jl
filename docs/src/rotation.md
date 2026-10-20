# Rotation

A two-dimensional rotation can be composed from directional advections. With
the velocity convention used by `advection!`, the x velocity is `y` and the y
velocity is `-x`:

```julia
using SemiLagrangian

mesh_x = UniformMesh(-5.0, 5.0, 256)
mesh_y = UniformMesh(-5.0, 5.0, 256)
interp_x = Lagrange(11)
interp_y = Lagrange(11)
x, y = mesh_x.points, mesh_y.points
f = [exp(-2 * (x0^2 + (y0 + 1.2)^2)) for x0 in x, y0 in y]

dt = 0.01
for _ in 1:100
    advection!(f, mesh_x, y, dt / 2, interp_x; dim = 1)
    advection!(f, mesh_y, -x, dt, interp_y; dim = 2)
    advection!(f, mesh_x, y, dt / 2, interp_x; dim = 1)
end
```

The directional API lets a simulation update coupled quantities between
steps. For example, a Vlasov-Poisson solver can recompute the electric field
after the x advection and before the v advection.
