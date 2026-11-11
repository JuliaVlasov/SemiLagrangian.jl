"""
Landau damping with a reusable PoissonAdvectionOperator.

Run with `julia --project=docs examples/landau-damping-poisson.jl`.
"""

using SemiLagrangian

function landau_poisson(nx, nv, dt, nt)
    meshx = UniformMesh(0.0, 4pi, nx)
    meshv = UniformMesh(-6.0, 6.0, nv)
    f = [
        (1 + 0.001cos(0.5x)) * exp(-v^2 / 2) / sqrt(2pi)
            for x in meshx.points, v in meshv.points
    ]
    position = AdvectionOperator(meshx, meshv.points, Spectral(nx); dim = 1)
    velocity = PoissonAdvectionOperator(meshx, meshv)
    coefficients = strangsplit(dt)
    time = zeros(nt + 1)
    energy = zeros(nt + 1)
    update_field!(velocity, f)
    energy[1] = sum(abs2, velocity.field) * step(meshx)

    for it in 1:nt
        # The velocity stage refreshes the field after the position half-step.
        split_advection!(f, position, velocity, coefficients)
        time[it + 1] = (it - 0.5) * dt
        energy[it + 1] = sum(abs2, velocity.field) * step(meshx)
    end
    return time, energy
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Plots
    time, energy = landau_poisson(64, 128, 0.05, 2000)
    display(plot(
        time, energy;
        xlabel = "Time",
        ylabel = "Electric field energy",
        yaxis = :log,
        label = "Spectral",
        title = "Landau damping: cached Poisson operator",
    ))
end
