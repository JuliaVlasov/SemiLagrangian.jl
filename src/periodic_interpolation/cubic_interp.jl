using FastInterpolations

export CubicSpline, interpolate!

"""
    CubicSpline(nx)

Cubic spline interpolation object for periodic 1D uniform grids.

# Arguments
- `nx::Int`: number of grid points in the periodic domain.

# Notes
This wrapper uses `FastInterpolations.jl`'s `cubic_interp!` under the hood
with periodic boundary conditions. It provides a small adapter exposing the
`PeriodicInterpolation1D`-style `interpolate!` API.

# Example
```julia
using PeriodicInterpolation1D, FastInterpolations

n = 100
u = sin.(2π .* (0:n-1) ./ n)
u_out = zeros(n)

cs = CubicSpline(n)
interpolate!(u_out, cs, u, 0.5)  # shift by 0.5 grid points
```
"""
struct CubicSpline
    nx :: Int
    CubicSpline(nx) = new(nx)
end


"""
    interpolate!(u_out, interp::CubicSpline, u, alpha)

Interpolate array `u` onto `u_out` shifted by `alpha` grid points using cubic
splines. This calls `cubic_interp!` from `FastInterpolations.jl` with periodic
boundary conditions.

# Arguments
- `u_out`: preallocated output array (length `interp.nx`)
- `interp`: `CubicSpline` interpolation object
- `u`: input array (length `interp.nx`)
- `alpha`: fractional shift in units of grid spacing (can be non-integer)
"""
function interpolate!(u_out, interp::CubicSpline, u, alpha)
    xi = 1:interp.nx
    xp = xi .+ alpha
    cubic_interp!(u_out, xi, u, xp, bc=PeriodicBC(endpoint=:exclusive))
    return u_out
end
