# Landau damping: interpolation comparison

## Problem description

Landau damping is the collisionless decay of a small electric-field
perturbation in a plasma. The one-dimensional, one-velocity-dimensional
(1D-1V) Vlasov-Poisson model evolves the particle distribution
\(f(x,v,t)\), where \(x\) is position, \(v\) is velocity, and \(E(x,t)\) is
the self-consistent electric field:

\[
\frac{\partial f}{\partial t}
 + v\frac{\partial f}{\partial x}
 + E\frac{\partial f}{\partial v} = 0,
\qquad
\frac{\partial E}{\partial x} = \rho,
\qquad
\rho(x,t) = \int f(x,v,t)\,dv
- \left\langle \int f(x,v,t)\,dv \right\rangle_x.
\]

Subtracting the spatial mean of the density enforces the periodic Poisson
problem's neutrality condition. The initial distribution is a weakly
perturbed Maxwellian,

\[
f(x,v,0) =
\frac{1 + \epsilon\cos(k_x x)}{\sqrt{2\pi}}
\exp\left(-\frac{v^2}{2}\right),
\qquad
\epsilon = 0.001,\quad k_x = 0.5.
\]

For this perturbation, linear theory predicts exponential decay of the
electric-field amplitude, \(E \propto \exp(-\gamma_L t)\), with
\(\gamma_L \approx 0.1533\). The plotted electric-field energy is the
spatial integral of \(E^2\), so its linear-theory decay rate is twice the
field-amplitude rate.

## Numerical method

The simulation uses Strang operator splitting to separate transport in
position and velocity:

1. Advect in position for half a time step with velocity \(v\).
2. Integrate \(f\) over velocity to obtain \(\rho\), then solve the periodic
   Poisson equation to update \(E\).
3. Advect in velocity for a full time step with acceleration \(E\).
4. Complete the time step with another half-step in position.

Each transport step is semi-Lagrangian. For each grid line, the method traces
characteristics backward from the arrival grid points, then interpolates
\(f\) at the resulting departure points. In this implementation, the
departure-point shift in grid-cell units is
\(\alpha = -\Delta t\,u/\Delta x\), where \(u\) is the velocity or electric
field for the corresponding advection.

## Interpolation methods compared

The comparison keeps the initial condition and simulation parameters fixed,
changing only the interpolation method used for both advection directions:

- **Spectral** — Fourier-based interpolation.
- **Periodic Lagrange** — periodic Lagrange interpolation with a 7-point
  stencil.
- **Periodic B-spline** — periodic B-spline interpolation of order 6.
- **Fast Lagrange** — local 7-point Lagrange interpolation.
- **FastInterpolations cubic** — periodic cubic interpolation from
  FastInterpolations.jl.

The run uses 64 position points, 128 velocity points, a time step of 0.1,
and 1,000 steps. The logarithmic plot compares electric-field energy over
time, and the script reports elapsed time for each method. Timings are
indicative: they depend on the machine and include compilation overhead.

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
