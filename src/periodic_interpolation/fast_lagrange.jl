"""
    FastLagrange interpolation module

    Provides local Lagrange polynomial interpolation over periodic domains using
    fixed-width stencils of 3, 5, 7, 9, or 11 points. Boundary points are handled
    via periodic (circular) wrapping of the input array.

    # ⚠️ Numerical Stability Disclaimer

    This implementation uses **local** (finite-stencil) Lagrange interpolation and is
    **not numerically stable** for large shifts or high-order stencils on coarse grids.

    The key limitations are:

    - **Locality**: Only a small neighbourhood of `2s` points (where `s` is half the
      stencil width) contributes to each interpolated value. If the function varies
      significantly between those points, accuracy degrades rapidly.
    - **Runge's phenomenon**: For large offsets `|p|`, the Lagrange weights grow in
      magnitude and oscillate in sign, leading to catastrophic cancellation.
      Best results are obtained for `p ∈ [-1, 1]` (i.e. interpolating within one
      grid cell of the reference point).
    - **No global spectral accuracy**: Unlike the FFT-based circulant approach
      (`lagrange.jl`), this method does not leverage the global smoothness of the
      function. It is designed for speed, not for maximum accuracy on smooth periodic
      functions.

    Use this module when performance is critical and the shift `p` is small (ideally
    sub-cell). For large shifts or high accuracy requirements on smooth functions,
    prefer the spectral circulant approach.

    # Precomputed inverse constants

    The `inv_*` constants below are reciprocals of the factorial-based denominators
    that appear in the Lagrange weight formulas. Precomputing them avoids repeated
    division at runtime.
"""

# ---------------------------------------------------------------------------
# Precomputed reciprocals of Lagrange denominator factorials
# These are used in the weight (coefficient) functions below to replace
# division by the known integer denominators with a multiply.
# ---------------------------------------------------------------------------
const inv_6       = 1.0 / 6.0        # 1/3!
const inv_12      = 1.0 / 12.0
const inv_24      = 1.0 / 24.0       # 1/4!
const inv_36      = 1.0 / 36.0
const inv_48      = 1.0 / 48.0
const inv_120     = 1.0 / 120.0      # 1/5!
const inv_144     = 1.0 / 144.0
const inv_240     = 1.0 / 240.0
const inv_576     = 1.0 / 576.0
const inv_720     = 1.0 / 720.0      # 1/6!
const inv_1440    = 1.0 / 1440.0
const inv_5040    = 1.0 / 5040.0     # 1/7!
const inv_14400   = 1.0 / 14400.0
const inv_17280   = 1.0 / 17280.0
const inv_30240   = 1.0 / 30240.0
const inv_40320   = 1.0 / 40320.0    # 1/8!
const inv_80640   = 1.0 / 80640.0
const inv_362880  = 1.0 / 362880.0   # 1/9!
const inv_3628800 = 1.0 / 3628800.0  # 1/10!

export FastLagrange

"""
    FastLagrange(stencil)

Interpolation object for local Lagrange interpolation on a periodic 1D grid.

# Arguments
- `stencil::Int`: Number of grid points used per interpolation. Must be one of
  `{3, 5, 7, 9, 11}`. Higher stencil widths give higher polynomial order but
  increase cost and can amplify instability for large offsets `p`.

| stencil | polynomial order | half-width |
|---------|-----------------|------------|
| 3       | 2               | 1          |
| 5       | 4               | 2          |
| 7       | 6               | 3          |
| 9       | 8               | 4          |
| 11      | 10              | 5          |
"""
struct FastLagrange
    stencil::Int
end

# ---------------------------------------------------------------------------
# 3-point (order-2) Lagrange weights and evaluation
# ---------------------------------------------------------------------------

"""
    lagr_3pt_coeff!(pp, p)

Compute the 3-point Lagrange interpolation weights for fractional offset `p`
and store them in the pre-allocated vector `pp` (length ≥ 3).

The stencil nodes are at relative positions {-1, 0, +1}.
The weights satisfy: f(p) ≈ pp[1]*f(-1) + pp[2]*f(0) + pp[3]*f(+1).

# Arguments
- `pp`: Pre-allocated coefficient vector (modified in-place)
- `p`:  Fractional offset in grid-spacing units; best accuracy for `p ∈ [-1, 1]`
"""
@inline function lagr_3pt_coeff!(pp, p)
    pp[1] = 0.5 * p * (p - 1.0)
    pp[2] = -(p + 1.0) * (p - 1.0)
    pp[3] = 0.5 * p * (p + 1.0)
end

"""
    lagr_3pt(fm1, f0, f1, p, pp)

Evaluate the 3-point Lagrange interpolant at offset `p` given function values
at stencil nodes {-1, 0, +1}. Coefficients `pp` must have been filled by
`lagr_3pt_coeff!` beforehand.

# Arguments
- `fm1, f0, f1`: Function values at grid points i-1, i, i+1
- `p`:  Fractional offset
- `pp`: Pre-computed Lagrange weights (from `lagr_3pt_coeff!`)
"""
@inline function lagr_3pt(fm1, f0, f1, p, pp)
    pp[1] * fm1 + pp[2] * f0 + pp[3] * f1
end

"""
    lagr_3pt_vec!(fi, fp, p, pp)

Apply 3-point Lagrange interpolation to the **interior** of array `fi`
(indices 2 to n-1), storing results in `fp`. Boundary points are not filled
here and must be handled separately with periodic wrap-around.

# Arguments
- `fi`: Input array (length n)
- `fp`: Output array (length n, modified in-place for interior points)
- `p`:  Fractional offset
- `pp`: Pre-computed Lagrange weights (from `lagr_3pt_coeff!`)
"""
@inline function lagr_3pt_vec!(fi, fp, p, pp)
    n = length(fi)
    for i = 2:(n-1)
        fp[i] = pp[1] * fi[i-1] + pp[2] * fi[i] + pp[3] * fi[i+1]
    end
end

# ---------------------------------------------------------------------------
# 5-point (order-4) Lagrange weights and evaluation
# ---------------------------------------------------------------------------

"""
    lagr_5pt_coeff!(pp, p)

Compute the 5-point Lagrange interpolation weights for fractional offset `p`.
Stencil nodes are at relative positions {-2, -1, 0, +1, +2}.
"""
@inline function lagr_5pt_coeff!(pp, p)
    pp[1] = p * (p - 1) * (p - 2) * (p + 1) * inv_24
    pp[2] = -p * (p - 1) * (p - 2) * (p + 2) * inv_6
    pp[3] = (p + 1) * (p + 2) * (p - 1) * (p - 2) * 0.25
    pp[4] = -p * (p + 1) * (p + 2) * (p - 2) * inv_6
    pp[5] = p * (p + 1) * (p + 2) * (p - 1) * inv_24
end

"""
    lagr_5pt(fm2, fm1, f0, f1, f2, p, pp)

Evaluate the 5-point Lagrange interpolant at offset `p`.
Stencil: {i-2, i-1, i, i+1, i+2}.
"""
@inline function lagr_5pt(fm2, fm1, f0, f1, f2, p, pp)
    pp[1] * fm2 + pp[2] * fm1 + pp[3] * f0 + pp[4] * f1 + pp[5] * f2
end

"""
    lagr_5pt_vec!(fi, fp, p, pp)

Apply 5-point Lagrange interpolation to interior points (indices 3 to n-2).
"""
@inline function lagr_5pt_vec!(fi, fp, p, pp)
    n = length(fi)
    for i = 3:(n-2)
        fp[i] =
            pp[1] * fi[i-2] +
            pp[2] * fi[i-1] +
            pp[3] * fi[i] +
            pp[4] * fi[i+1] +
            pp[5] * fi[i+2]
    end
end

# ---------------------------------------------------------------------------
# 7-point (order-6) Lagrange weights and evaluation
# ---------------------------------------------------------------------------

"""
    lagr_7pt_coeff!(pp, p)

Compute the 7-point Lagrange interpolation weights for fractional offset `p`.
Stencil nodes at relative positions {-3, -2, -1, 0, +1, +2, +3}.
"""
@inline function lagr_7pt_coeff!(pp, p)
    pp[1] = (p + 2) * (p + 1) * p * (p - 1) * (p - 2) * (p - 3) * inv_720
    pp[2] = -(p + 3) * (p + 1) * p * (p - 1) * (p - 2) * (p - 3) * inv_120
    pp[3] = (p + 3) * (p + 2) * p * (p - 1) * (p - 2) * (p - 3) * inv_48
    pp[4] = -(p + 3) * (p + 2) * (p + 1) * (p - 1) * (p - 2) * (p - 3) * inv_36
    pp[5] = (p + 3) * (p + 2) * (p + 1) * p * (p - 2) * (p - 3) * inv_48
    pp[6] = -(p + 3) * (p + 2) * (p + 1) * p * (p - 1) * (p - 3) * inv_120
    pp[7] = (p + 3) * (p + 2) * (p + 1) * p * (p - 1) * (p - 2) * inv_720
end

"""
    lagr_7pt(fm3, fm2, fm1, f0, f1, f2, f3, p, pp)

Evaluate the 7-point Lagrange interpolant at offset `p`.
Stencil: {i-3, ..., i+3}.
"""
@inline function lagr_7pt(fm3, fm2, fm1, f0, f1, f2, f3, p, pp)
    return pp[1] * fm3 +
           pp[2] * fm2 +
           pp[3] * fm1 +
           pp[4] * f0 +
           pp[5] * f1 +
           pp[6] * f2 +
           pp[7] * f3
end

"""
    lagr_7pt_vec!(fi, fp, p, pp)

Apply 7-point Lagrange interpolation to interior points (indices 4 to n-3).
"""
@inline function lagr_7pt_vec!(fi, fp, p, pp)
    n = length(fi)
    for i = 4:(n-3)
        fp[i] =
            pp[1] * fi[i-3] +
            pp[2] * fi[i-2] +
            pp[3] * fi[i-1] +
            pp[4] * fi[i] +
            pp[5] * fi[i+1] +
            pp[6] * fi[i+2] +
            pp[7] * fi[i+3]
    end
end

# ---------------------------------------------------------------------------
# 9-point (order-8) Lagrange weights and evaluation
# ---------------------------------------------------------------------------

"""
    lagr_9pt_coeff!(pp, p)

Compute the 9-point Lagrange interpolation weights for fractional offset `p`.
Stencil nodes at relative positions {-4, -3, -2, -1, 0, +1, +2, +3, +4}.

The weights are expressed using the factored form
    wᵢ = ∏_{j≠i} (p - j) / (i - j)
with the products over squared differences written as `(p² - k²)` for efficiency.
"""
@inline function lagr_9pt_coeff!(pp, p)
    pp[1] = p * (p - 4) * (p^2 - 9) * (p^2 - 4) * (p^2 - 1) * inv_40320
    pp[2] = -p * (p - 3) * (p^2 - 16) * (p^2 - 4) * (p^2 - 1) * inv_5040
    pp[3] = p * (p - 2) * (p^2 - 16) * (p^2 - 9) * (p^2 - 1) * inv_1440
    pp[4] = -p * (p - 1) * (p^2 - 16) * (p^2 - 9) * (p^2 - 4) * inv_720
    pp[5] = (p^2 - 16) * (p^2 - 9) * (p^2 - 4) * (p^2 - 1) * inv_576
    pp[6] = -(p + 1) * p * (p^2 - 16) * (p^2 - 9) * (p^2 - 4) * inv_720
    pp[7] = (p + 2) * p * (p^2 - 16) * (p^2 - 9) * (p^2 - 1) * inv_1440
    pp[8] = -(p + 3) * p * (p^2 - 16) * (p^2 - 4) * (p^2 - 1) * inv_5040
    pp[9] = (p + 4) * p * (p^2 - 9) * (p^2 - 4) * (p^2 - 1) * inv_40320
end

"""
    lagr_9pt(fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, p, pp)

Evaluate the 9-point Lagrange interpolant at offset `p`.
Stencil: {i-4, ..., i+4}.
"""
@inline function lagr_9pt(fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, p, pp)
    return pp[1] * fm4 +
           pp[2] * fm3 +
           pp[3] * fm2 +
           pp[4] * fm1 +
           pp[5] * f0 +
           pp[6] * f1 +
           pp[7] * f2 +
           pp[8] * f3 +
           pp[9] * f4
end

"""
    lagr_9pt_vec!(fi, fp, p, pp)

Apply 9-point Lagrange interpolation to interior points (indices 5 to n-4).
"""
@inline function lagr_9pt_vec!(fi, fp, p, pp)
    n = length(fi)
    for i = 5:(n-4)
        fp[i] =
            pp[1] * fi[i-4] +
            pp[2] * fi[i-3] +
            pp[3] * fi[i-2] +
            pp[4] * fi[i-1] +
            pp[5] * fi[i] +
            pp[6] * fi[i+1] +
            pp[7] * fi[i+2] +
            pp[8] * fi[i+3] +
            pp[9] * fi[i+4]
    end
    return
end

# ---------------------------------------------------------------------------
# 11-point (order-10) Lagrange weights and evaluation
# ---------------------------------------------------------------------------

"""
    lagr_11pt_coeff!(pp, p)

Compute the 11-point Lagrange interpolation weights for fractional offset `p`.
Stencil nodes at relative positions {-5, ..., 0, ..., +5}.
"""
@inline function lagr_11pt_coeff!(pp, p)
    pp[1] =
        p * (p - 5.0) * (p^2 - 16.0) * (p^2 - 9.0) * (p^2 - 4.0) * (p^2 - 1.0) * inv_3628800
    pp[2] =
        -p * (p - 4.0) * (p^2 - 25.0) * (p^2 - 9.0) * (p^2 - 4.0) * (p^2 - 1.0) * inv_362880
    pp[3] =
        p * (p - 3.0) * (p^2 - 25.0) * (p^2 - 16.0) * (p^2 - 4.0) * (p^2 - 1.0) * inv_80640
    pp[4] =
        -p * (p - 2.0) * (p^2 - 25.0) * (p^2 - 16.0) * (p^2 - 9.0) * (p^2 - 1.0) * inv_30240
    pp[5] =
        p * (p - 1.0) * (p^2 - 25.0) * (p^2 - 16.0) * (p^2 - 9.0) * (p^2 - 4.0) * inv_17280
    pp[6] =
        -(p^2 - 25.0) * (p^2 - 16.0) * (p^2 - 9.0) * (p^2 - 4.0) * (p^2 - 1.0) * inv_14400
    pp[7] =
        (p + 1.0) * p * (p^2 - 25.0) * (p^2 - 16.0) * (p^2 - 9.0) * (p^2 - 4.0) * inv_17280
    pp[8] =
        -(p + 2.0) * p * (p^2 - 25.0) * (p^2 - 16.0) * (p^2 - 9.0) * (p^2 - 1.0) * inv_30240
    pp[9] =
        (p + 3.0) * p * (p^2 - 25.0) * (p^2 - 16.0) * (p^2 - 4.0) * (p^2 - 1.0) * inv_80640
    pp[10] =
        -(p + 4.0) * p * (p^2 - 25.0) * (p^2 - 9.0) * (p^2 - 4.0) * (p^2 - 1.0) * inv_362880
    pp[11] =
        (p + 5.0) * p * (p^2 - 16.0) * (p^2 - 9.0) * (p^2 - 4.0) * (p^2 - 1.0) * inv_3628800
end

"""
    lagr_11pt(fm5, fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, f5, p, pp)

Evaluate the 11-point Lagrange interpolant at offset `p`.
Stencil: {i-5, ..., i+5}.
"""
@inline function lagr_11pt(fm5, fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, f5, p, pp)
    return pp[1] * fm5 +
           pp[2] * fm4 +
           pp[3] * fm3 +
           pp[4] * fm2 +
           pp[5] * fm1 +
           pp[6] * f0 +
           pp[7] * f1 +
           pp[8] * f2 +
           pp[9] * f3 +
           pp[10] * f4 +
           pp[11] * f5
end

"""
    lagr_11pt_vec!(fi, fp, p, pp)

Apply 11-point Lagrange interpolation to interior points (indices 6 to n-5).
"""
@inline function lagr_11pt_vec!(fi, fp, p, pp)
    n = length(fi)
    for i = 6:(n-5)
        fp[i] =
            pp[1] * fi[i-5] +
            pp[2] * fi[i-4] +
            pp[3] * fi[i-3] +
            pp[4] * fi[i-2] +
            pp[5] * fi[i-1] +
            pp[6] * fi[i] +
            pp[7] * fi[i+1] +
            pp[8] * fi[i+2] +
            pp[9] * fi[i+3] +
            pp[10] * fi[i+4] +
            pp[11] * fi[i+5]
    end
    return
end

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

"""
    interpolate!(fp, interpolant::FastLagrange, fi, p)

Interpolate the periodic array `fi` at offset `p` (in grid-spacing units)
and store the result in `fp`.

Boundary points are handled by **periodic wrapping**: the array is treated as
circular, so indices outside `1:n` wrap around modulo `n`.

The interior is processed by the corresponding `lagr_Npt_vec!` kernel; the
`stencil/2 - 1` boundary points at each end are filled individually using
explicitly wrapped index expressions.

# Arguments
- `fp`:           Output array (length n, modified in-place)
- `interpolant`:  `FastLagrange` object specifying the stencil width
- `fi`:           Input array (length n, periodic)
- `p`:            Fractional shift in units of the grid spacing `dx`.
                  Accuracy is best for `p ∈ [-1, 1]`; large values lead to
                  weight blow-up and loss of precision (see module disclaimer).

# Supported stencil widths
`3`, `5`, `7`, `9`, `11`  (polynomial orders 2, 4, 6, 8, 10 respectively)

# Example
```julia
n   = 128
fi  = sin.(2π .* (0:n-1) ./ n)
fp  = similar(fi)
lag = FastLagrange(7)
pp  = zeros(7)
lagr_7pt_coeff!(pp, 0.3)
interpolate!(fp, lag, fi, 0.3)
```
"""
function interpolate!(fp, interpolant::FastLagrange, fi, p)

    n = length(fi)
    stencil = interpolant.stencil
    pp = zeros(stencil)

    return if stencil == 7
        lagr_7pt_coeff!(pp, p)
        fp[1] = lagr_7pt(fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], p, pp)
        fp[2] = lagr_7pt(fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], p, pp)
        fp[3] = lagr_7pt(fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], p, pp)
        lagr_7pt_vec!(fi, fp, p, pp)
        fp[n-2] = lagr_7pt(fi[n-5], fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], p, pp)
        fp[n-1] = lagr_7pt(fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], p, pp)
        fp[n]   = lagr_7pt(fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], p, pp)
    elseif stencil == 5
        lagr_5pt_coeff!(pp, p)
        fp[1] = lagr_5pt(fi[n-1], fi[n], fi[1], fi[2], fi[3], p, pp)
        fp[2] = lagr_5pt(fi[n], fi[1], fi[2], fi[3], fi[4], p, pp)
        lagr_5pt_vec!(fi, fp, p, pp)
        fp[n-1] = lagr_5pt(fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], p, pp)
        fp[n]   = lagr_5pt(fi[n-2], fi[n-1], fi[n], fi[1], fi[2], p, pp)
    elseif stencil == 3
        lagr_3pt_coeff!(pp, p)
        fp[1] = lagr_3pt(fi[n], fi[1], fi[2], p, pp)
        lagr_3pt_vec!(fi, fp, p, pp)
        fp[n] = lagr_3pt(fi[n-1], fi[n], fi[1], p, pp)
    elseif stencil == 9
        lagr_9pt_coeff!(pp, p)
        fp[1] = lagr_9pt(fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], p, pp)
        fp[2] = lagr_9pt(fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], p, pp)
        fp[3] = lagr_9pt(fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], fi[7], p, pp)
        fp[4] = lagr_9pt(fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], fi[7], fi[8], p, pp)
        lagr_9pt_vec!(fi, fp, p, pp)
        fp[n-3] = lagr_9pt(fi[n-7], fi[n-6], fi[n-5], fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], p, pp)
        fp[n-2] = lagr_9pt(fi[n-6], fi[n-5], fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], p, pp)
        fp[n-1] = lagr_9pt(fi[n-5], fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], p, pp)
        fp[n]   = lagr_9pt(fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], p, pp)
    elseif stencil == 11
        lagr_11pt_coeff!(pp, p)
        fp[1] = lagr_11pt(fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], p, pp)
        fp[2] = lagr_11pt(fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], fi[7], p, pp)
        fp[3] = lagr_11pt(fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], fi[7], fi[8], p, pp)
        fp[4] = lagr_11pt(fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], fi[7], fi[8], fi[9], p, pp)
        fp[5] = lagr_11pt(fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], fi[7], fi[8], fi[9], fi[10], p, pp)
        lagr_11pt_vec!(fi, fp, p, pp)
        fp[n-4] = lagr_11pt(fi[n-9], fi[n-8], fi[n-7], fi[n-6], fi[n-5], fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], p, pp)
        fp[n-3] = lagr_11pt(fi[n-8], fi[n-7], fi[n-6], fi[n-5], fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], p, pp)
        fp[n-2] = lagr_11pt(fi[n-7], fi[n-6], fi[n-5], fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], p, pp)
        fp[n-1] = lagr_11pt(fi[n-6], fi[n-5], fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], p, pp)
        fp[n]   = lagr_11pt(fi[n-5], fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], p, pp)
    end
end

