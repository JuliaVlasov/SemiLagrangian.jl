"""
    SplineCache(x; order=3)

Precompute an open B-spline interpolation system for the strictly increasing
grid `x`. `order` is the B-spline degree and must satisfy `1 <= order < length(x)`.
Reuse the cache to interpolate multiple data vectors on the same grid.
"""
struct SplineCache{T<:AbstractFloat,F}
    x::Vector{T}
    knots::Vector{T}
    order::Int
    factor::F
end

function SplineCache(x::AbstractVector{<:Real}; order::Int = 3)
    1 <= order < length(x) ||
        throw(ArgumentError("order must satisfy 1 <= order < length(x)"))
    all(i -> x[i] < x[i + 1], 1:(length(x) - 1)) ||
        throw(ArgumentError("x must be strictly increasing"))

    T = typeof(float(zero(eltype(x))))
    grid = T.(x)
    knots = Vector{T}(undef, length(grid) + order + 1)
    fill!(view(knots, 1:(order + 1)), first(grid))
    fill!(view(knots, (length(knots) - order):length(knots)), last(grid))
    @inbounds for j in 1:(length(grid) - order - 1)
        knots[order + 1 + j] = sum(view(grid, (j + 1):(j + order))) / order
    end

    matrix = zeros(T, length(grid), length(grid))
    @inbounds for row in eachindex(grid)
        first_control, values = _spline_basis(knots, order, grid[row], length(grid))
        for j in eachindex(values)
            matrix[row, first_control + j - 1] = values[j]
        end
    end
    factor = lu(matrix)
    return SplineCache{T,typeof(factor)}(grid, knots, order, factor)
end

"""
    spline_interp(x, y, xq; order=3)
    spline_interp(x, y, x_targets; order=3)

Construct an interpolating open B-spline of degree `order` through `(x, y)`
and evaluate it at one or more target points. `x` must be strictly increasing;
targets outside the grid throw `DomainError`. For repeated interpolation on
the same grid, construct `SplineCache(x; order)` once and reuse it.
"""
function spline_interp(
    x::AbstractVector{<:Real},
    y::AbstractVector,
    xq::Real;
    order::Int = 3,
)
    return spline_interp(SplineCache(x; order), y, xq)
end

function spline_interp(
    x::AbstractVector{<:Real},
    y::AbstractVector,
    x_targets::AbstractArray{<:Real};
    order::Int = 3,
)
    return spline_interp(SplineCache(x; order), y, x_targets)
end

function spline_interp(cache::SplineCache, y::AbstractVector, xq::Real)
    _check_spline_data(cache, y)
    _check_spline_query(cache, xq)
    coefficients = cache.factor \ y
    return _spline_eval(cache, coefficients, xq)
end

function spline_interp(cache::SplineCache, y::AbstractVector, x_targets::AbstractArray{<:Real})
    T = promote_type(eltype(cache.x), eltype(y), eltype(x_targets))
    T = T <: Integer ? float(T) : T
    output = similar(x_targets, T)
    return spline_interp!(output, cache, y, x_targets)
end

"""
    spline_interp!(output, x, y, x_targets; order=3)
    spline_interp!(output, cache, y, x_targets)

Evaluate an interpolating B-spline at each target and store the results in
`output`. The output must have the same shape as `x_targets`.
"""
function spline_interp!(
    output::AbstractArray,
    x::AbstractVector{<:Real},
    y::AbstractVector,
    x_targets::AbstractArray{<:Real};
    order::Int = 3,
)
    return spline_interp!(output, SplineCache(x; order), y, x_targets)
end

function spline_interp!(
    output::AbstractArray,
    cache::SplineCache,
    y::AbstractVector,
    x_targets::AbstractArray{<:Real},
)
    _check_spline_data(cache, y)
    size(output) == size(x_targets) ||
        throw(DimensionMismatch("output and x_targets must have the same size"))
    @inbounds for xq in x_targets
        _check_spline_query(cache, xq)
    end

    coefficients = cache.factor \ y
    @inbounds for i in eachindex(output, x_targets)
        output[i] = _spline_eval(cache, coefficients, x_targets[i])
    end
    return output
end

function _check_spline_data(cache, y)
    length(cache.x) == length(y) ||
        throw(DimensionMismatch("y must have the same length as the cached grid"))
    return nothing
end

function _check_spline_query(cache, xq)
    first(cache.x) <= xq <= last(cache.x) ||
        throw(DomainError(xq, "target must lie within the interpolation grid"))
    return nothing
end

function _spline_eval(cache::SplineCache, coefficients, xq)
    first_control, basis = _spline_basis(cache.knots, cache.order, xq, length(cache.x))
    result = zero(promote_type(eltype(coefficients), eltype(basis)))
    @inbounds for j in eachindex(basis)
        result += basis[j] * coefficients[first_control + j - 1]
    end
    return result
end

function _spline_basis(knots, order, xq, ncontrols)
    span = xq == knots[end] ? ncontrols - 1 : searchsortedlast(knots, xq) - 1
    first_control = span - order + 1
    T = promote_type(eltype(knots), typeof(xq))
    basis = zeros(T, order + 1)
    left = zeros(T, order)
    right = zeros(T, order)
    basis[1] = one(T)

    @inbounds for j in 1:order
        left[j] = xq - knots[span + 2 - j]
        right[j] = knots[span + 1 + j] - xq
        saved = zero(T)
        for r in 0:(j - 1)
            denominator = right[r + 1] + left[j - r]
            term = iszero(denominator) ? zero(T) : basis[r + 1] / denominator
            basis[r + 1] = saved + right[r + 1] * term
            saved = left[j - r] * term
        end
        basis[j + 1] = saved
    end
    return first_control, basis
end
