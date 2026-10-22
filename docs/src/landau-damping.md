# Landau Damping

This chapter presents the one-dimensional, one-velocity-dimensional
(1D-1V) Vlasov-Poisson problem, its semi-Lagrangian solution with Strang
splitting, and a comparison of interpolation methods for the Landau-damping
test case.

## Problem description

Landau damping is the collisionless decay of a small electric-field
perturbation in a plasma. The one-dimensional, one-velocity-dimensional
(1D-1V) Vlasov-Poisson model evolves the particle distribution
``f(x,v,t)``, where ``x`` is position, ``v`` is velocity, and ``E(x,t)`` is
the self-consistent electric field:


```math
\frac{\partial f}{\partial t}
 + v\frac{\partial f}{\partial x}
 + E\frac{\partial f}{\partial v} = 0,
\qquad
\frac{\partial E}{\partial x} = \rho,
\qquad
\rho(x,t) = \int f(x,v,t)\,dv
- \left\langle \int f(x,v,t)\,dv \right\rangle_x.
```

Subtracting the spatial mean of the density enforces the neutrality
condition for the periodic Poisson problem. In the Vlasov equation, the
term ``v\,\partial f/\partial x`` transports the distribution in position,
while ``E\,\partial f/\partial v`` transports it in velocity under the
self-consistent electric field.

The Landau-damping initial condition is a weakly perturbed Maxwellian:

```math
f(x,v,0) =
\frac{1 + \epsilon\cos(k_x x)}{\sqrt{2\pi}}
\exp\left(-\frac{v^2}{2}\right),
\qquad
\epsilon = 0.001,\quad k_x = 0.5.
```

For this perturbation, linear theory predicts exponential decay of the
electric-field amplitude, ``E \propto \exp(-\gamma_L t)``, with
``\gamma_L \approx 0.1533``. The plotted electric-field energy is the
spatial integral of ``E^2``, so its linear-theory decay rate is twice the
field-amplitude rate.

## Numerical method: Strang splitting

The Vlasov equation is split into its position- and velocity-advection
operators. A Strang step composes them symmetrically:

1. Advect in position for half a time step with velocity ``v``.
2. Integrate ``f`` over velocity to obtain ``\rho``, then solve the periodic
   Poisson equation to update ``E``.
3. Advect in velocity for a full time step with acceleration ``E``.
4. Complete the time step with another half-step in position.

Each directional transport step uses the semi-Lagrangian method: trace
characteristics backward from each arrival grid point, then interpolate
``f`` at the resulting departure point. For a direction with velocity ``u``,
time step ``\Delta t``, and grid spacing ``\Delta x``, the shift is
``\alpha = -\Delta t\,u/\Delta x``, measured in grid cells.

## Interpolation methods compared

The comparison keeps the problem, initial condition, and simulation
parameters fixed, changing only the interpolation method used for both
advection directions:

- **Spectral** — Fourier-based interpolation.
- **LagrangeFFT** — periodic FFT-based Lagrange interpolation with a 7-point
  stencil.
- **FastLagrange** — local 7-point Lagrange interpolation.
- **FastInterpolations cubic** — periodic cubic interpolation from
  FastInterpolations.jl.

The run uses 64 position points, 128 velocity points, a time step of 0.05,
and 2,000 steps. The smaller step keeps the largest spatial characteristic
shift below one grid cell for the local FastLagrange interpolator:
``\Delta t\,v_{\max}/(2\Delta x) < 1``. Its interpolation weights are
intended for fractional shifts of about one cell or less; using a larger
shift directly can make the high-order polynomial weights unstable. The
step count is doubled to keep the final simulation time unchanged.

The logarithmic plot compares electric-field energy over time, and the
script reports elapsed time for each method. Timings are indicative: they
depend on the machine and include compilation overhead.

## Running the example

The full, commented source is shared with the standalone repository example
so that the documentation and runnable script stay in sync:

```@eval
import Markdown
Markdown.parse(
    "```julia\n" *
    read(joinpath(@__DIR__, "..", "..", "examples", "landau-damping.jl"), String) *
    "\n```",
)
```

Run the comparison from the repository root with:

```sh
julia --project examples/landau-damping.jl
```

```@example landau
include(joinpath(@__DIR__, "..", "..", "examples", "landau-damping.jl"))
```

## References

- Vlasov-Poisson system: Krall & Trivelpiece (1973)
- Landau damping: Landau (1946)
- Semi-Lagrangian schemes: Cheng & Knorr (1976), Sonnendrucker et al. (1999)
- SeLaLib (Fortran library): https://selalib.github.io/
