"""
    HermiteCache(x)

Cache grid geometry for repeated cubic Hermite interpolation on `x`.
The cached interval widths and their reciprocals can be reused with different
value and derivative vectors.
"""
struct HermiteCache{TX<:Real,TH<:Real}
    x::Vector{TX}
    widths::Vector{TH}
    inverse_widths::Vector{TH}
end

function HermiteCache(x::AbstractVector{<:Real})
    length(x) >= 2 || throw(ArgumentError("x must contain at least two points"))
    all(i -> x[i] < x[i + 1], 1:(length(x) - 1)) ||
        throw(ArgumentError("x must be strictly increasing"))

    grid = collect(x)
    TH = eltype(grid) <: Integer ? typeof(float(one(eltype(grid)))) : eltype(grid)
    widths = Vector{TH}(undef, length(grid) - 1)
    inverse_widths = similar(widths)
    @inbounds for i in eachindex(widths)
        widths[i] = TH(grid[i + 1]) - TH(grid[i])
        inverse_widths[i] = inv(widths[i])
    end
    return HermiteCache{eltype(grid),TH}(grid, widths, inverse_widths)
end

"""
    hermite_interp(x, y, dy, xq)
    hermite_interp(x, y, dy, x_targets)

Evaluate the piecewise cubic Hermite interpolant through values `y` and
derivatives `dy` sampled on the strictly increasing grid `x`. Derivatives must
be with respect to the coordinates in `x`. Targets outside the grid domain
throw `DomainError`.

For repeated interpolation on the same grid, create `HermiteCache(x)` once and
reuse it with `hermite_interp(cache, y, dy, x_targets)`.
"""
function hermite_interp(
    x::AbstractVector{<:Real},
    y::AbstractVector,
    dy::AbstractVector,
    xq::Real,
)
    return hermite_interp(HermiteCache(x), y, dy, xq)
end

function hermite_interp(
    x::AbstractVector{<:Real},
    y::AbstractVector,
    dy::AbstractVector,
    x_targets::AbstractArray{<:Real},
)
    return hermite_interp(HermiteCache(x), y, dy, x_targets)
end

function hermite_interp(cache::HermiteCache, y::AbstractVector, dy::AbstractVector, xq::Real)
    _check_hermite_data(cache, y, dy)
    _check_hermite_query(cache, xq)
    return _hermite_eval(cache, y, dy, xq)
end

function hermite_interp(
    cache::HermiteCache,
    y::AbstractVector,
    dy::AbstractVector,
    x_targets::AbstractArray{<:Real},
)
    T = _hermite_output_type(cache, y, dy, eltype(x_targets))
    output = similar(x_targets, T)
    return hermite_interp!(output, cache, y, dy, x_targets)
end

"""
    hermite_interp!(output, x, y, dy, x_targets)
    hermite_interp!(output, cache, y, dy, x_targets)

In-place, batched version of [`hermite_interp`](@ref). `output` must have the
same shape as `x_targets`.
"""
function hermite_interp!(
    output::AbstractArray,
    x::AbstractVector{<:Real},
    y::AbstractVector,
    dy::AbstractVector,
    x_targets::AbstractArray{<:Real},
)
    return hermite_interp!(output, HermiteCache(x), y, dy, x_targets)
end

function hermite_interp!(
    output::AbstractArray,
    cache::HermiteCache,
    y::AbstractVector,
    dy::AbstractVector,
    x_targets::AbstractArray{<:Real},
)
    _check_hermite_data(cache, y, dy)
    size(output) == size(x_targets) ||
        throw(DimensionMismatch("output and x_targets must have the same size"))

    @inbounds for i in eachindex(output, x_targets)
        xq = x_targets[i]
        _check_hermite_query(cache, xq)
        output[i] = _hermite_eval(cache, y, dy, xq)
    end
    return output
end

function _check_hermite_data(cache, y, dy)
    length(cache.x) == length(y) == length(dy) ||
        throw(DimensionMismatch("x, y, and dy must have the same length"))
    return nothing
end

function _check_hermite_query(cache, xq)
    first(cache.x) <= xq <= last(cache.x) ||
        throw(DomainError(xq, "target must lie within the interpolation grid"))
    return nothing
end

function _hermite_output_type(cache, y, dy, Tq)
    T = promote_type(eltype(cache.widths), eltype(y), eltype(dy), Tq)
    return T <: Integer ? float(T) : T
end

function _hermite_eval(cache::HermiteCache, y, dy, xq)
    x = cache.x
    @inbounds xq == x[end] && return y[end]
    interval = min(searchsortedlast(x, xq), length(x) - 1)
    @inbounds xq == x[interval] && return y[interval]

    h = cache.widths[interval]
    t = (xq - x[interval]) * cache.inverse_widths[interval]
    t2 = t * t
    t3 = t2 * t
    h00 = 2t3 - 3t2 + 1
    h10 = t3 - 2t2 + t
    h01 = -2t3 + 3t2
    h11 = t3 - t2
    @inbounds return h00 * y[interval] + h10 * h * dy[interval] +
                   h01 * y[interval + 1] + h11 * h * dy[interval + 1]
end
