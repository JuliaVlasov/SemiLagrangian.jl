"""
    lagrange_interp(x, y, xq; order=3)
    lagrange_interp(x, y, x_targets; order=3)

Evaluate local Lagrange interpolation of degree `order` at one or more target
points. The input grid `x` must be strictly increasing and contain at least
`order + 1` points. Each target uses a contiguous stencil of `order + 1`
samples, centered around the target where possible and shifted at the endpoints.
Targets outside `x`'s domain throw `DomainError`.

For repeated interpolation on the same grid, construct `LagrangeCache(x;
order)` once and call `lagrange_interp(cache, y, targets)` to reuse the
precomputed stencil weights across data vectors and target arrays.

"""
function lagrange_interp(x::AbstractVector{<:Real}, y::AbstractVector, xq::Real; order::Int = 3)
    return lagrange_interp(LagrangeCache(x; order), y, xq)
end

function lagrange_interp(
    x::AbstractVector{<:Real},
    y::AbstractVector,
    x_targets::AbstractArray{<:Real};
    order::Int = 3,
)
    return lagrange_interp(LagrangeCache(x; order), y, x_targets)
end

"""
    lagrange_interp!(output, x, y, x_targets; order=3)

In-place version of [`lagrange_interp`](@ref). `output` must have the same
shape as `x_targets`; targets outside the grid domain throw `DomainError`.
"""
function lagrange_interp!(
    output::AbstractArray,
    x::AbstractVector{<:Real},
    y::AbstractVector,
    x_targets::AbstractArray{<:Real};
    order::Int = 3,
)
    return lagrange_interp!(output, LagrangeCache(x; order), y, x_targets)
end

"""
    LagrangeCache(x; order=3)

Precompute the barycentric node weights for every stencil on the grid `x`.
Reuse the cache with [`lagrange_interp`](@ref) or [`lagrange_interp!`](@ref)
to interpolate multiple data vectors or target arrays on the same grid.
"""
struct LagrangeCache{TX<:Real,TW<:Real}
    x::Vector{TX}
    order::Int
    weights::Matrix{TW}
end

function LagrangeCache(x::AbstractVector{<:Real}; order::Int = 3)
    _check_lagrange_grid(x, order)
    grid = collect(x)
    nstencils = length(grid) - order
    TW = eltype(grid) <: Integer ? typeof(float(one(eltype(grid)))) : eltype(grid)
    weights = Matrix{TW}(undef, nstencils, order + 1)
    @inbounds for first_index = 1:nstencils, i = 0:order
        xi = TW(grid[first_index + i])
        weight = one(TW)
        for j = 0:order
            i == j && continue
            weight /= xi - TW(grid[first_index + j])
        end
        weights[first_index, i + 1] = weight
    end
    return LagrangeCache{eltype(grid),TW}(grid, order, weights)
end

function lagrange_interp(cache::LagrangeCache, y::AbstractVector, xq::Real)
    _check_lagrange_data(cache, y)
    first(cache.x) <= xq <= last(cache.x) ||
        throw(DomainError(xq, "target must lie within the interpolation grid"))
    return _lagrange_eval(cache, y, xq)
end

function lagrange_interp(
    cache::LagrangeCache,
    y::AbstractVector,
    x_targets::AbstractArray{<:Real},
)
    T = promote_type(eltype(cache.x), eltype(y), eltype(x_targets))
    T = T <: Integer ? float(T) : T
    output = similar(x_targets, T)
    return lagrange_interp!(output, cache, y, x_targets)
end

function lagrange_interp!(
    output::AbstractArray,
    cache::LagrangeCache,
    y::AbstractVector,
    x_targets::AbstractArray{<:Real},
)
    _check_lagrange_data(cache, y)
    size(output) == size(x_targets) ||
        throw(DimensionMismatch("output and x_targets must have the same size"))

    @inbounds for i in eachindex(output, x_targets)
        xq = x_targets[i]
        first(cache.x) <= xq <= last(cache.x) ||
            throw(DomainError(xq, "target must lie within the interpolation grid"))
        output[i] = _lagrange_eval(cache, y, xq)
    end
    return output
end

function _check_lagrange_grid(x, order)
    order >= 0 || throw(ArgumentError("order must be nonnegative"))
    length(x) >= order + 1 ||
        throw(ArgumentError("x needs at least order + 1 samples"))
    all(i -> x[i] < x[i + 1], 1:(length(x) - 1)) ||
        throw(ArgumentError("x must be strictly increasing"))
    return nothing
end

function _check_lagrange_data(cache, y)
    length(cache.x) == length(y) ||
        throw(DimensionMismatch("y must have the same length as the cached grid"))
    return nothing
end

function _lagrange_eval(cache::LagrangeCache, y, xq)
    x = cache.x
    order = cache.order
    first_index = clamp(searchsortedlast(x, xq) - order ÷ 2, 1, length(x) - order)
    last_index = first_index + order
    T = promote_type(eltype(cache.weights), eltype(y), typeof(xq))

    @inbounds for i in first_index:last_index
        xq == x[i] && return y[i]
    end

    numerator = zero(T)
    denominator = zero(T)
    @inbounds for i in first_index:last_index
        weight = cache.weights[first_index, i - first_index + 1]
        term = weight / (xq - x[i])
        numerator += term * y[i]
        denominator += term
    end
    return numerator / denominator
end
