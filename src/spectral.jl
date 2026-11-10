export Spectral, interpolate!

"""
    Spectral(nx)

Create a spectral (Fourier) interpolation object for periodic 1D grids.

Spectral interpolation uses the Fourier representation to interpolate with exponential
accuracy for smooth periodic functions. This is the simplest and most accurate method
for very smooth data but assumes periodicity.

# Arguments
- `nx::Int`: Grid size (number of points)

# Fields
- `nx::Int`: Grid size
- `eigvals::Vector{ComplexF64}`: Fourier phase factors (working array)
- `ufft::Vector{ComplexF64}`: FFT workspace

# Example
```julia
spec = Spectral(100)  # spectral interpolation on 100-point grid
```
"""
struct Spectral{order} <: AbstractShiftInterpolation{order}

    nx::Int
    eigvals::Vector{ComplexF64}
    ufft::Vector{ComplexF64}

    function Spectral(nx::Int)

        ufft = zeros(ComplexF64, nx)
        eigvals = zeros(ComplexF64, nx)

        return new{nx - 1}(nx, eigvals, ufft)

    end

end

"""
    interpolate!(u_out, u, alpha, spec::Spectral)

Interpolate array `u` by a displacement `alpha` (in grid units) using pure spectral
(Fourier-based) interpolation.

# Arguments
- `u_out::AbstractVector`: Output array (will be modified in-place)
- `u::AbstractVector`: Input array (periodic boundary conditions assumed)
- `alpha::Real`: Shift amount in grid-point units
- `spec::Spectral`: Spectral interpolation object

# Description
This method uses Fourier space phase shifting:
1. Compute FFT of input
2. Apply phase shift e^(i * 2π * k * alpha / nx) to each frequency mode
3. Inverse FFT to obtain interpolated result

Provides exponential accuracy for smooth periodic functions. **Caution**: only use for
smooth data; sharp features or discontinuities will exhibit Gibbs phenomena.

# Example
```julia
spec = Spectral(100)
f = sin.(2π .* (0:99) ./ 100)
f_interp = zeros(100)
interpolate!(f_interp, f, 0.25, spec)  # shift by 0.25 grid points
```
"""
function interpolate!(u_out, u, alpha::Real, interpolant::Spectral)

    nx = interpolant.nx
    alpha = Float64(alpha)
    interpolant.ufft .= u
    fft!(interpolant.ufft)

    interpolant.eigvals[1] = 1.0
    interpolant.eigvals[nx ÷ 2 + 1] = exp(1im * π * alpha)
    for k in 1:(nx ÷ 2 - 1)
        interpolant.eigvals[k + 1] = exp(1im * 2pi * k * alpha / nx)
        interpolant.eigvals[nx - k + 1] = exp(-1im * 2pi * k * alpha / nx)
    end

    interpolant.ufft .*= interpolant.eigvals

    return u_out .= real(ifft(interpolant.ufft))

end
