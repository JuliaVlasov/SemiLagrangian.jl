# Vlasov-Poisson 1D1V Simulation Example

This example demonstrates the simulation of the **Vlasov-Poisson system** in 1D phase space
(1 spatial dimension + 1 velocity dimension) and compares the performance and accuracy of
different periodic interpolation methods available in `SemiLagrangian`.

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
2. **PeriodicLagrange**: Global FFT-based Lagrange (7-point stencil)
3. **Spectral**: Pure Fourier interpolation
4. **PeriodicBSpline**: B-spline interpolation (order 6)

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

```@example vp1d1v
using Plots
using FFTW
using SemiLagrangian
using LinearAlgebra
import Statistics: mean

```


```@example vp1d1v
function compute_rho(
        meshv::UniformMesh,
        f::Array{Float64, 2}
    )

    dv = step(meshv)
    rho = dv * sum(f, dims = 2)
    return vec(rho .- mean(rho))
end

" compute Ex using that -ik*Ex = rho "
function compute_e(meshx::UniformMesh, rho::Vector{Float64})
    nx = length(meshx)
    k = 2pi / meshx.width
    modes = zeros(Float64, nx)
    modes .= k * fftfreq(nx, nx)
    modes[1] = 1.0
    rhok = fft(rho) ./ modes
    rhok .*= -1im
    ifft!(rhok)
    return real(rhok)
end
```

# Landau Damping

[Landau damping - Wikipedia](https://en.wikipedia.org/wiki/Landau_damping)

```@example vp1d1v
function landau(dt, nt::Int64)

    # Set grid
    p = 7
    nx, nv = 128, 256
    xmin, xmax = 0.0, 4pi
    vmin, vmax = -6.0, 6.0
    meshx = UniformMesh(xmin, xmax, nx)
    meshv = UniformMesh(vmin, vmax, nv)
    x = meshx.points
    v = meshv.points
    dx = step(meshx)

    eps, kx = 0.001, 0.5
    f = zeros(Float64, (nx, nv))
    f .= (1.0 .+ eps * cos.(kx * x)) / sqrt(2π) .* transpose(exp.(-0.5 * v .^ 2))
    fᵗ = zeros(Float64, (nv, nx))
    interpolant = FastLagrange(p)

    rho = compute_rho(meshv, f)
    e = compute_e(meshx, rho)
    ℰ = Float64[0.5 * log(sum(e .* e) * dx)]
    t = Float64[0.0]
    for it in 1:nt
        advection!(f, meshx, v, 0.5dt, interpolant; dim = 1)
        rho = compute_rho(meshv, f)
        e = compute_e(meshx, rho)
        push!(ℰ, 0.5 * log(sum(e .* e) * dx))
        push!(t, (it - 0.5) * dt)
        transpose!(fᵗ, f)
        advection!(fᵗ, meshv, e, dt, interpolant; dim = 1)
        transpose!(f, fᵗ)
        advection!(f, meshx, v, 0.5dt, interpolant; dim = 1)
    end

    return t, ℰ

end

```

```@example vp1d1v
nt = 5000
dt = 0.01
@time t, nrj = landau(dt, nt)
plot(t, nrj; label = "E")
plot!(t, -0.1533 * t .- 5.5; label = "-0.1533t.-5.5")
```
