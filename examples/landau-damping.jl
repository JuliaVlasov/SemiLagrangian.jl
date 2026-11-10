"""
Compare interpolation methods in a 1D-1V Landau-damping simulation.

The Vlasov-Poisson model evolves a distribution function `f(x, v, t)` under
transport in position and velocity, coupled through the electric field. The
initial condition is a small sinusoidal perturbation of a Maxwellian:
`f(x, v, 0) = (1 + 0.001 cos(0.5x)) exp(-v^2 / 2) / sqrt(2pi)`.

Each time step uses Strang splitting: advect half a step in position, update
the charge density and electric field, advect a full step in velocity, then
finish with a half-step in position. Each advection is semi-Lagrangian: trace
the characteristics backward and interpolate at their departure points.

The plot compares electric-field energy for the same grid, initial condition,
and time step, changing only the interpolation method. It includes spectral,
FFT-based Lagrange, local Lagrange, and cubic spline interpolation through the
`CubicSpline` wrapper. Reported timings include Julia compilation and depend
on the machine.

Run from the repository root with `julia --project examples/landau-damping.jl`.
"""

using FFTW
using LinearAlgebra
using Plots
using SemiLagrangian

include(joinpath(@__DIR__, "compute_rho.jl"))
include(joinpath(@__DIR__, "compute_e.jl"))

function landau(nx, nv, dt, nt, interpolant_x, interpolant_v)
    meshx = UniformMesh(0.0, 4π, nx)
    meshv = UniformMesh(-6.0, 6.0, nv)

    x = meshx.points
    v = meshv.points
    dx = step(meshx)
    # Store f with position along rows and velocity along columns.
    f = (1.0 .+ 0.001 .* cos.(0.5 .* x)) ./ sqrt(2π) .*
        transpose(exp.(-0.5 .* v .^ 2))
    fᵗ = similar(f, nv, nx)

    # Store the initial energy and then sample the field at each step midpoint.
    rho = compute_rho(meshv, f)
    e = compute_e(meshx, rho)
    energy = [sum(abs2, e) * dx]
    time = [0.0]

    for it in 1:nt
        # Strang splitting: half-step in x, full-step in v, then half-step in x.
        # Recompute the electric field between the two advection directions.
        advection!(f, meshx, v, 0.5dt, interpolant_x; dim = 1)
        rho = compute_rho(meshv, f)
        e = compute_e(meshx, rho)
        push!(energy, sum(abs2, e) * dx)
        # This field is evaluated at the midpoint and drives the velocity step.
        push!(time, (it - 0.5) * dt)
        # Transpose so that velocity is the interpolated (row) direction.
        transpose!(fᵗ, f)
        advection!(fᵗ, meshv, e, dt, interpolant_v; dim = 1)
        transpose!(f, fᵗ)
        advection!(f, meshx, v, 0.5dt, interpolant_x; dim = 1)
    end

    return time, energy
end

nx, nv = 64, 128
# Keep the final time fixed while choosing dt so local Lagrange's characteristic
# shift in x stays within one grid cell: dt * maximum(abs, v) / (2dx) < 1.
dt, nt = 0.05, 2000
# Keep the simulation parameters fixed and vary only the interpolation method.
methods = [
    ("Spectral", Spectral(nx), Spectral(nv)),
    ("Lagrange FFT", LagrangeFFT(nx, 7), LagrangeFFT(nv, 7)),
    ("Local Lagrange", LocalLagrange(7), LocalLagrange(7)),
    ("Cubic spline", CubicSpline(nx), CubicSpline(nv)),
]

results = map(methods) do (name, interpolant_x, interpolant_v)
    @info "Landau damping with $name interpolation"
    @time time, energy = landau(nx, nv, dt, nt, interpolant_x, interpolant_v)
    (name, time, energy)
end

comparison_plot = plot(;
    xlabel = "Time",
    ylabel = "Electric field energy",
    yaxis = :log,
    legend = :best,
    title = "Landau damping: interpolation method comparison",
)
for (name, time, energy) in results
    plot!(comparison_plot, time, energy; label = name)
end
comparison_plot
