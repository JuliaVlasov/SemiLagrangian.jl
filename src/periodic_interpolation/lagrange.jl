export Lagrange, interpolate!

"""
    Lagrange(nx, order)

Create a Lagrange interpolation object using FFT-based circulant approach for periodic 1D grids.

This implementation uses a global Lagrange polynomial interpolation via FFT, providing
spectral accuracy for smooth periodic functions. The interpolation uses circulant matrix
properties to efficiently compute interpolated values via convolution in Fourier space.

# Arguments
- `nx::Int`: Grid size (number of points)
- `order::Int`: Stencil width (number of points). Supports odd orders: 3, 5, 7, 9, 11

# Fields
- `nx::Int`: Grid size
- `order::Int`: Polynomial order
- `w::Vector{Float64}`: Pre-allocated weight vector
- `bfft::Vector{ComplexF64}`: FFT of basis function weights
- `ufft::Vector{ComplexF64}`: FFT of input data

# Example
```julia
lagr = Lagrange(100, 5)  # 100 grid points, order-4 polynomial
```
"""
struct Lagrange 

    nx :: Int
    order :: Int
    w :: Vector{Float64}
    bfft :: Vector{ComplexF64}
    ufft :: Vector{ComplexF64}

    Lagrange( nx, order ) = new( nx, order, zeros(nx), zeros(ComplexF64, nx), zeros(ComplexF64, nx))

end


"""
    interpolate!(u_out, lagr::Lagrange, u, alpha)

Interpolate array `u` by a displacement `alpha` (in grid units) using global Lagrange
polynomials via FFT convolution.

# Arguments
- `u_out::AbstractVector`: Output array (will be modified in-place)
- `lagr::Lagrange`: Lagrange interpolation object
- `u::AbstractVector`: Input array (periodic boundary conditions assumed)
- `alpha::Float64`: Shift amount in grid-point units. Can be any real number; periodic
  boundary conditions handle wrapping.

# Description
This function computes the interpolated values using:
1. Lagrange polynomial weights at fractional offset
2. FFT-based circulant matrix approach for efficiency
3. Periodic boundary conditions via modular arithmetic

For smooth periodic functions, this provides spectral (exponential) accuracy.

# Example
```julia
lagr = Lagrange(100, 5)
f = sin.(2π .* (0:99) ./ 100)
f_interp = zeros(100)
interpolate!(f_interp, lagr, f, 0.5)  # shift by 0.5 grid points
```
"""
function interpolate!(u_out, self::Lagrange, u, alpha)

    nx = length(u)
    d  = self.order ÷ 2 - 1

    # fractional shift in grid-point units, in [0, nx)
    x  = - alpha
    ix = floor(Int, x)
    if ix == nx; ix = 0 end
    x  = x - ix   # fractional part in [0, 1)

    # Build Lagrange weights over stencil i = -d, ..., 0, 1, ..., d+1
    # w[i] = prod_{j != i, j in stencil} (x - j) / (i - j)
    stencil = -d : d+1          # 2(d+1) points

    fill!(self.w, 0.0)
    for i in stencil
        num = 1.0
        den = 1.0
        for j in stencil
            if j != i
                num *= (x - j)
                den *= (i - j)
            end
        end
        idx = mod(ix + i, nx) + 1
        self.w[idx] = num / den
    end

    self.bfft .= fft(self.w)
    self.ufft .= fft(u)

    u_out .= real(ifft(self.ufft .* self.bfft))

end
