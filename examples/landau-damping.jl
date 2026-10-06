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
periodic Lagrange, periodic B-spline, fast Lagrange, and FastInterpolations
cubic interpolation. Reported timings include Julia compilation and depend on
the machine.

Run from the repository root with `julia --project examples/landau-damping.jl`.
"""

import FastInterpolations
using FFTW
using LinearAlgebra
using Plots
using SemiLagrangian

include(joinpath(@__DIR__, "uniform_mesh.jl"))
include(joinpath(@__DIR__, "compute_rho.jl"))
include(joinpath(@__DIR__, "compute_e.jl"))

# FastInterpolations uses its own API, so this marker dispatches to it from
# the common advection routine without changing the other interpolation calls.
struct FastInterpolationsCubic end

function interpolate_column!(out, interpolant, values, alpha, xi)
    interpolate!(out, interpolant, values, alpha)
end

function interpolate_column!(out, ::FastInterpolationsCubic, values, alpha, xi)
    # FastInterpolations expects a uniform grid; the exclusive periodic
    # endpoint means every sample represents a distinct grid point.
    FastInterpolations.cubic_interp!(
        out,
        xi,
        values,
        xi .+ alpha;
        bc=FastInterpolations.PeriodicBC(endpoint=:exclusive),
    )
end

function advection!(f, mesh::UniformMesh, velocity, dt, interpolant)
    xi = 1.0:mesh.nx
    fp = similar(view(f, :, 1))

    for j in eachindex(velocity)
        fi = view(f, :, j)
        # Trace characteristics backward to obtain the shift in grid-cell units.
        alpha = -dt * velocity[j] / mesh.dx
        interpolate_column!(fp, interpolant, fi, alpha, xi)
        fi .= fp
    end
end

function landau(nx, nv, dt, nt, interpolant_x, interpolant_v)
    meshx = UniformMesh(0.0, 4π, nx)
    meshv = UniformMesh(-6.0, 6.0, nv)

    x = meshx.x
    v = meshv.x
    dx = meshx.dx
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
        advection!(f, meshx, v, 0.5dt, interpolant_x)
        rho = compute_rho(meshv, f)
        e = compute_e(meshx, rho)
        push!(energy, sum(abs2, e) * dx)
        # This field is evaluated at the midpoint and drives the velocity step.
        push!(time, (it - 0.5) * dt)
        # Transpose so that velocity is the interpolated (row) direction.
        transpose!(fᵗ, f)
        advection!(fᵗ, meshv, e, dt, interpolant_v)
        transpose!(f, fᵗ)
        advection!(f, meshx, v, 0.5dt, interpolant_x)
    end

    return time, energy
end

nx, nv = 64, 128
# Keep the final time fixed while choosing dt so FastLagrange's characteristic
# shift in x stays within one grid cell: dt * maximum(abs, v) / (2dx) < 1.
dt, nt = 0.05, 2000
# Keep the simulation parameters fixed and vary only the interpolation method.
methods = [
    ("Spectral", Spectral(nx), Spectral(nv)),
    ("Periodic Lagrange", PeriodicLagrange(nx, 7), PeriodicLagrange(nv, 7)),
    ("Periodic B-spline", PeriodicBSpline(nx, 6), PeriodicBSpline(nv, 6)),
    ("Fast Lagrange", FastLagrange(7), FastLagrange(7)),
    ("FastInterpolations cubic", FastInterpolationsCubic(), FastInterpolationsCubic()),
]

results = map(methods) do (name, interpolant_x, interpolant_v)
    @info "Landau damping with $name interpolation"
    @time time, energy = landau(nx, nv, dt, nt, interpolant_x, interpolant_v)
    (name, time, energy)
end

comparison_plot = plot(;
    xlabel="Time",
    ylabel="Electric field energy",
    yaxis=:log,
    legend=:best,
    title="Landau damping: interpolation method comparison",
)
for (name, time, energy) in results
    plot!(comparison_plot, time, energy; label=name)
end
comparison_plot
