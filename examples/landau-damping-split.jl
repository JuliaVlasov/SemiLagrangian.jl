"""
Landau damping with AdvectionOperator and split_advection!.

Run from the repository root with
`julia --project=docs examples/landau-damping-split.jl`.
"""

using SemiLagrangian

include(joinpath(@__DIR__, "compute_rho.jl"))
include(joinpath(@__DIR__, "compute_e.jl"))

function landau_split(nx, nv, dt, nt)
    meshx = UniformMesh(0.0, 4π, nx)
    meshv = UniformMesh(-6.0, 6.0, nv)
    x, v = meshx.points, meshv.points
    f = (1.0 .+ 0.001 .* cos.(0.5 .* x)) ./ sqrt(2π) .*
        transpose(exp.(-0.5 .* v .^ 2))
    f_transposed = similar(f, nv, nx)

    interpolant_x = Spectral(nx)
    interpolant_v = Spectral(nv)
    position = AdvectionOperator(meshx, v, interpolant_x; dim = 2)
    coefficients = strangsplit(dt)

    e = compute_e(meshx, compute_rho(meshv, f))
    energy = [sum(abs2, e) * step(meshx)]
    time = [0.0]

    for it in 1:nt
        # The self-consistent field must be updated after the first half-step.
        advection!(f, meshx, v, coefficients[1], interpolant_x; dim = 1)
        e = compute_e(meshx, compute_rho(meshv, f))
        push!(energy, sum(abs2, e) * step(meshx))
        push!(time, (it - 0.5) * dt)

        permutedims!(f_transposed, f, (2, 1))
        velocity = AdvectionOperator(meshv, e, interpolant_v; dim = 1)
        split_advection!(f_transposed, velocity, position, @view(coefficients[2:3]))
        permutedims!(f, f_transposed, (2, 1))
    end

    return time, energy
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Plots
    time, energy = landau_split(64, 128, 0.05, 2000)
    display(plot(
        time,
        energy;
        xlabel = "Time",
        ylabel = "Electric field energy",
        yaxis = :log,
        label = "Spectral",
        title = "Landau damping",
    ))
end
