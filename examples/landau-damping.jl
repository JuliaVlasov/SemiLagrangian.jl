"""
    Vlasov-Poisson 1D1V Simulation Example

This example demonstrates the simulation of the **Vlasov-Poisson system** in 1D phase space
(1 spatial dimension + 1 velocity dimension) and compares the performance and accuracy of
different interpolation methods provided by PeriodicInterpolation1D.

## Problem Description

The Vlasov-Poisson system models the evolution of a collisionless plasma:

    ∂f/∂t + v·∇ₓf + E·∇ᵥf = 0    (Vlasov equation)
    ∇²ϕ = ρ = ∫f dv                 (Poisson equation)
    E = -∇ₓϕ                        (Electric field)

where:
- `f(x,v,t)`: distribution function in phase space
- `x`: spatial coordinate
- `v`: velocity coordinate  
- `E`: electric field
- `ρ`: charge density

## Numerical Method

**Operator Splitting**: The time evolution is split into alternating advection steps:

1. **X-advection**: Move particles in position space with velocity `v`
   `∂f/∂t + v·∇ₓf = 0`

2. **V-advection**: Move particles in velocity space with acceleration from electric field
   `∂f/∂t + E·∇ᵥf = 0`

The solution is computed by alternating these advections (Strang splitting):
- Advect in x by half time-step
- Advect in v by full time-step (updating E based on current charge density)
- Advect in x by half time-step

## Interpolation Role

Each advection step is implemented via **semi-Lagrangian interpolation**:
- Trace characteristics backward along velocity vectors
- Interpolate distribution function at characteristic endpoints
- This approach avoids numerical dissipation from upwinding schemes

## Test Case: Landau Damping

This example implements **Landau damping**, a classic weakly nonlinear plasma phenomenon:

    f₀(x,v) = [1 + ε·cos(kₓ·x)] / √(2π) · exp(-v²/2)

where:
- `ε = 0.001`: Small perturbation amplitude
- `kₓ = 0.5`: Wavenumber

**Expected behavior**: The electric field energy decays exponentially as
    E(t) ∝ exp(-γ_L · t)
where `γ_L ≈ 0.1533` is the Landau damping rate.

## Comparison of Methods

The script compares four interpolation methods:

1. **FastLagrange**: Locally computed 7-point Lagrange polynomials
2. **Lagrange**: Global FFT-based Lagrange (7-point stencil)
3. **Spectral**: Pure Fourier interpolation
4. **BSpline**: B-spline interpolation (order 6)

Each method shows:
- Execution time for the simulation
- Energy evolution over time
- Accuracy compared to theoretical decay rate (linear fit: `-0.1533*t - 5.5`)

## References

- Vlasov-Poisson system: Krall & Trivelpiece (1973)
- Landau damping: Landau (1946)
- Semi-Lagrangian schemes: Cheng & Knorr (1976), Sonnendrucker et al. (1999)
- SeLaLib (Fortran library): https://selalib.github.io/

## Expected Output

- Performance comparison of different interpolation methods
- Energy decay plots comparing theoretical and numerical solutions
- Demonstration of semi-Lagrangian advection with periodic boundary conditions
"""

using Plots
using PeriodicInterpolation1D
using LinearAlgebra

include("uniform_mesh.jl")
include("compute_rho.jl")
include("compute_e.jl")

function advection!(
        f::Array{Float64, 2},
        p::Int64, mesh::UniformMesh, v::Vector{Float64},
        nv::Int64, dt::Float64, interpolant
    )

    nx = mesh.nx
    dx = mesh.dx

    fp = zeros(nx)
    for j in 1:nv
        fi = view(f, :, j)
        alpha = - dt * v[j] / dx
        interpolate!(fp, interpolant, fi, alpha)
        f[:, j] .= fp
    end

    return
end

function landau(nx, nv, dt, nt::Int64, interpolant_x, interpolant_v)

    # Set grid
    xmin, xmax = 0.0, 4pi
    vmin, vmax = -6.0, 6.0
    meshx = UniformMesh(xmin, xmax, nx)
    meshv = UniformMesh(vmin, vmax, nv)
    x = meshx.x
    v = meshv.x
    dx = meshx.dx

    eps, kx = 0.001, 0.5
    f = zeros(Float64, (nx, nv))
    f .= (1.0 .+ eps * cos.(kx * x)) / sqrt(2π) .* transpose(exp.(-0.5 * v .^ 2))
    fᵗ = zeros(Float64, (nv, nx))

    rho = compute_rho(meshv, f)
    e = compute_e(meshx, rho)
    ℰ = Float64[0.5 * log(sum(e .* e) * dx)]
    t = Float64[0.0]

    for it in 1:nt
        advection!(f, p, meshx, v, nv, 0.5dt, interpolant_x)
        rho = compute_rho(meshv, f)
        e = compute_e(meshx, rho)
        push!(ℰ, 0.5 * log(sum(e .* e) * dx))
        push!(t, (it - 0.5) * dt)
        transpose!(fᵗ, f)
        advection!(fᵗ, p, meshv, e, nx, dt, interpolant_v)
        transpose!(f, fᵗ)
        advection!(f, p, meshx, v, nv, 0.5dt, interpolant_x)
    end

    return t, ℰ

end

@info "Fast Lagrange"
p = 7
nx, nv = 64, 128
dt, nt = 0.01, 10000
landau(nx, nv, dt, 1, FastLagrange(p), FastLagrange(p))
@time t, nrj = landau(nx, nv, dt, nt, FastLagrange(p), FastLagrange(p))
plot(t, -0.1533 * t .- 5.5; label = "-0.1533t.-5.5")
plot!(t, nrj; label = "Fast Lagrange")
@info "Lagrange"
dt, nt = 0.1, 1000
landau(nx, nv, dt, 1, Lagrange(nx, p), Lagrange(nv, p))
@time t, nrj = landau(nx, nv, dt, nt, Lagrange(nx, p), Lagrange(nv, p))
plot(t, -0.1533 * t .- 5.5; label = "-0.1533t.-5.5")
plot!(t, nrj; label = "Lagrange")
@info "Spectral"
dt, nt = 0.1, 1000
landau(nx, nv, dt, 1, Spectral(nx), Spectral(nv))
@time t, nrj = landau(nx, nv, dt, nt, Spectral(nx), Spectral(nv))
plot!(t, nrj; label = "Spectral")
@info "B-splines"
dt, nt = 0.1, 1000
landau(nx, nv, dt, 1, BSpline(nx, 6), BSpline(nv, 6))
@time t, nrj = landau(nx, nv, dt, nt, BSpline(nx, 6), BSpline(nv, 6))
plot!(t, nrj; label = "B-splines")
@info "Cubic-splines"
dt, nt = 0.1, 1000
landau(nx, nv, dt, 1, CubicSpline(nx), CubicSpline(nv))
@time t, nrj = landau(nx, nv, dt, nt, CubicSpline(nx), CubicSpline(nv))
plot!(t, nrj; label = "Cubic-splines")
