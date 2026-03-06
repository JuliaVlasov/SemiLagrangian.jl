function uniform_bsplines_eval_basis(p::Int, x::Float64)
    if p < 0
        return Float64[1.0]
    end

    bspl = zeros(Float64, p + 1)
    bspl[1] = 1.0

    @inbounds for j = 1:p
        xx = -x
        inv_j = 1.0 / j
        saved = 0.0
        for r = 1:j
            xx += 1.0
            temp = bspl[r] * inv_j
            bspl[r] = saved + xx * temp
            saved = (j - xx) * temp
        end
        bspl[j+1] = saved
    end

    return bspl
end

export BSpline, interpolate!

"""
    BSpline(nx, order)

Create a B-spline interpolation object for periodic 1D grids.

B-spline interpolation provides smooth interpolation with controllable accuracy via the
spline order. This implementation uses FFT for efficient periodic boundary handling via
circulant matrix properties.

# Arguments
- `nx::Int`: Grid size (number of points)
- `order::Int`: **Even** spline degree. Common choices: 2, 4, 6, 8
  - Order 2: linear B-splines
  - Order 4: cubic B-splines (smooth C² continuity)
  - Order 6: quintic B-splines (smooth C⁴ continuity)

# Fields
- `nx::Int`: Grid size
- `order::Int`: Spline order
- `eigvals_M::Vector{Float64}`: Eigenvalues of the B-spline mass matrix
- `eigvals_S::Vector{ComplexF64}`: Eigenvalues at displaced points (working array)
- `eikx::Vector{ComplexF64}`: Precomputed phase factors
- `ufft::Vector{ComplexF64}`: FFT workspace

# Errors
- Raises an error if `order` is odd

# Example
```julia
bspl = BSpline(100, 4)  # cubic B-splines on 100-point grid
```
"""
struct BSpline

    nx::Int
    order::Int
    eigvals_M::Vector{Float64}
    eigvals_S::Vector{ComplexF64}
    eikx::Vector{ComplexF64}
    ufft::Vector{ComplexF64}

    function BSpline(nx::Int, order::Int)

        p = order - 1

        isodd(order) && error("Spline interpolation order needs to be even. Order here is: $order")

        ufft = zeros(ComplexF64, nx)
        eigvals_M = zeros(Float64, nx)
        eigvals_S = zeros(ComplexF64, nx)
        eikx = exp.([2π * i * 1im / nx for i in 0:nx-1])

        # compute eigenvalues of degree p b-spline matrix
        biatx = uniform_bsplines_eval_basis(p, 0.0)
        mid = div(p + 1, 2)
        eigvals_M .= biatx[mid]
        
        @inbounds for i in 1:(mid-1)
            coeff = biatx[i + mid] * 2
            for j in 1:nx
                eigvals_M[j] += coeff * cos(2π * i * (j - 1) / nx)
            end
        end

        @inbounds for i in eachindex(eigvals_M)
            eigvals_M[i] = 1.0 / eigvals_M[i]
        end

        new(nx, order, eigvals_M, eigvals_S, eikx, ufft)

    end

end

@inline modulo(a, p) = a - floor(Int, a / p) * p

"""
    interpolate!(u_out, bspl::BSpline, u, alpha)

Interpolate array `u` by a displacement `alpha` (in grid units) using B-spline basis.

# Arguments
- `u_out::AbstractVector`: Output array (will be modified in-place)
- `bspl::BSpline`: B-spline interpolation object
- `u::AbstractVector`: Input array (periodic boundary conditions assumed)
- `alpha::Float64`: Shift amount in grid-point units

# Description
This function computes interpolation using:
1. B-spline basis functions evaluated at the fractional offset
2. FFT for efficient circulant matrix-vector products
3. Periodic boundary conditions

Provides smooth interpolation with accuracy proportional to the spline order.

# Example
```julia
bspl = BSpline(100, 4)
f = sin.(2π .* (0:99) ./ 100)
f_interp = zeros(100)
interpolate!(f_interp, bspl, f, 0.3)  # shift by 0.3 grid points
```
"""
function interpolate!( u_out, interpolant::BSpline, u, alpha::Float64 )

   p = interpolant.order - 1
   nx = interpolant.nx

   interpolant.ufft .= u
   fft!(interpolant.ufft)
    
   # compute eigenvalues of cubic splines evaluated at displaced points
   ishift = floor(Int, alpha)
   beta   = -ishift + alpha
   biatx = uniform_bsplines_eval_basis(p, beta)
   fill!(interpolant.eigvals_S, 0.0im)
   
   mid = div(p + 1, 2)
   half_p = div(p - 1, 2)
   eikx = interpolant.eikx
   eigvals_S = interpolant.eigvals_S
   
   @inbounds for i in -half_p:mid
       coeff = biatx[i + mid]
       idx_offset = ishift + i
       for j in 1:nx
           imode = modulo(idx_offset * (j - 1), nx) + 1
           eigvals_S[j] += coeff * eikx[imode]
       end
   end
          
   # compute interpolating spline using fft and properties of circulant matrices
   @inbounds interpolant.ufft .*= interpolant.eigvals_S .* interpolant.eigvals_M
        
   u_out .= real(ifft(interpolant.ufft))
    

end
