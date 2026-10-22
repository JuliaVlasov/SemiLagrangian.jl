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
initial = [exp(-2 * (x0^2 + (y0 + 1.2)^2)) for x0 in x, y0 in y]

op_x = AdvectionOperator(mesh_x, y, interp_x; dim = 1)
op_y = AdvectionOperator(mesh_y, -x, interp_y; dim = 2)

function rotation_error(initial, op_x, op_y, nsteps; rotations = 10, splitting = magicsplit)
    dt = 2pi / nsteps
    coefficients = splitting(dt)
    f = copy(initial)
    for _ in 1:(rotations * nsteps)
        split_advection!(f, op_x, op_y, coefficients)
    end
    return maximum(abs, f .- initial)
end
```

### Strang splitting

Strang splitting is a second-order approximation of the rotation flow. The
error is measured after ten full rotations, when the exact solution returns
to its initial state:

```julia
strang_error = rotation_error(initial, op_x, op_y, 100; splitting = strangsplit)
```

### Magic splitting

For this rotation, `magicsplit(dt)` gives the shear coefficients
`[tan(dt / 2), sin(dt), tan(dt / 2)]`. Their composition is exactly the
rotation by `dt` at the continuous operator level:

```julia
magic_error = rotation_error(initial, op_x, op_y, 100; splitting = magicsplit)
```

The magic coefficients introduce no splitting error for this rotation; the
measured numerical error still includes interpolation and finite-domain errors.

The directional API lets a simulation update coupled quantities between
steps. For example, a Vlasov-Poisson solver can recompute the electric field
after the x advection and before the v advection.
