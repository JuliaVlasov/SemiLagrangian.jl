using LinearAlgebra
import Base: isless

"""
Abstract supertype for the boundary conditions used by line interpolants.

A boundary condition decides how an interpolation stencil of `order + 1`
points is built when it reaches the ends of the data line. The concrete
subtypes are [`PeriodicBC`](@ref) and [`ExtrapolationBC`](@ref). It is used as the
`bc` type parameter of [`AbstractInterpolation`](@ref) and selected with the
`bc` keyword of the constructors, for example
`Lagrange(5; bc = ExtrapolationBC())`.
"""
abstract type AbstractBoundaryCondition end

"""
    PeriodicBC()

Periodic boundary condition: the data line is considered as one period of a
periodic function, so interpolation stencils wrap around the ends of the line
(the point after the last one is the first one). This is the default boundary
condition of [`Lagrange`](@ref), [`Hermite`](@ref) and [`BSplineLU`](@ref),
and the only one supported by [`BSplineFFT`](@ref).

# Example
```julia
interp = Lagrange(5; bc = PeriodicBC())
```
"""
struct PeriodicBC <: AbstractBoundaryCondition end

"""
    ExtrapolationBC()

Non-periodic boundary condition: interpolation stencils never leave the data
line. Near the ends the stencil is shifted inside the line and the
interpolating polynomial is used to extrapolate, so the result is less
accurate close to the boundaries than in the interior. This is the boundary
condition to use for non-periodic variables (for example a velocity
coordinate).

# Example
```julia
interp = Lagrange(5; bc = ExtrapolationBC())
```
"""
struct ExtrapolationBC <: AbstractBoundaryCondition end

"""
Abstract supertype of all the line interpolation methods
([`Lagrange`](@ref), [`Hermite`](@ref), [`BSplineLU`](@ref), [`BSplineFFT`](@ref)).

An interpolation of order `order` evaluates the function at a shifted position
using a stencil of `order + 1` grid points. The shift is split into an integer
part (in grid cells) and a fractional part in `[0, 1)`; the weights of the
stencil only depend on the fractional part and are given by a table of
polynomial functions.

# Type parameters
- `T` : type of number on which the interpolation works (`Float64`, `BigFloat`, ...)
- `bc` : boundary condition type, [`PeriodicBC`](@ref) or [`ExtrapolationBC`](@ref)
- `order` : order of the interpolation (stencil of `order + 1` points)

# Implementation constraints
- `tabfct::Vector` : this field must be present in the implementation, it is a
  table of `order + 1` functions of the fractional shift that give the weights
  of the stencil.
- Methods needing a preprocessing of the line (such as B-splines, which solve
  a linear system to get the spline coefficients) must specialize `sol`,
  `sol!` and `issolidentity`.
"""
abstract type AbstractInterpolation{T,bc,order} end

"""
Abstract supertype of the periodic interpolants that apply one constant shift to
a whole line (`LagrangeFFT`, `Spectral`, `FastLagrange`, `CubicSpline`).

They are [`AbstractInterpolation`](@ref)s with `Float64` numbers and
[`PeriodicBC`](@ref), but have no stencil table: they implement
`interpolate!(output, interpolant, input, shift)` directly, and the generic
`interpolate!(output, input, shift, interpolant)` forwards to it.
"""
abstract type AbstractShiftInterpolation{order} <: AbstractInterpolation{Float64,PeriodicBC,order} end

"""
   
Return the order of the interpolation `interp`, i.e. the stencil used contains
`order + 1` grid points.

# Example
```julia
get_order(Lagrange(5))  # returns 5
```
"""
get_order(_::AbstractInterpolation{T,bc,order}) where {T,bc,order} = order

"""
Interface method that transforms the data line before the stencil weights are
applied. By default it does nothing and returns the line unchanged. B-spline
interpolations override it to compute the spline coefficients from the values.

# Arguments
- `_::AbstractInterpolation` : interpolation implementation
- `b::AbstractArray` : line (or array, treated line by line) to transform

# Returns
The transformed line.
"""
sol(_::AbstractInterpolation, b::AbstractArray) = b

issolidentity(_::AbstractInterpolation) = true

isbspline(_::AbstractInterpolation) = false

Base.show(io::IO, interp::AbstractInterpolation) = print(io, typeof(interp))

function sol!(
    Y::AbstractArray{T2,N},
    interp_t::AbstractVector{I},
    b::AbstractArray{T2,N},
) where {T,N,I<:AbstractInterpolation{T},T2<:Union{T,OpTuple{N,T}}}
    sz = size(Y)
    p = circshift(1:N, -1)
    perm = p
    bufout = b
    for i = 1:N
        interp = interp_t[i]
        if issolidentity(interp)
            if i == 1
                bufout = b
            else
                perm = perm[p]
            end
        else
            bufin = (i == 1) ? b : permutedims(bufout, perm)
            bufout = zeros(T2, sz)
            for c in CartesianIndices(sz[2:end])
                bo = view(bufout, :, c)
                bi = view(bufin, :, c)
                sol!(bo, interp, bi)
            end
            perm = p
        end
        sz = sz[p]
    end
    permutedims!(Y, bufout, perm)
    return Y
end


"""
Apply the preprocessing `sol` of each interpolation of `interp_t` to the
array `b`, one dimension after the other (`interp_t[i]` acts along dimension
`i`). If all the interpolations are the identity transformation, `b` is
returned without copy, otherwise a new array is returned.
"""
function sol(
    interp_t::AbstractVector{I},
    b::Union{AbstractArray{T,N},AbstractArray{Complex{T},N},AbstractArray{OpTuple{N,T},N}},
) where {T<:Real,N,I<:AbstractInterpolation{T}}
    if all(issolidentity.(interp_t))
        return b
    else
        return sol!(zeros(eltype(b), size(b)), interp_t, b)
    end
end

@inline function getprecal(interp::AbstractInterpolation{T}, decf::T) where {T}
    return [T(fct(decf)) for fct in interp.tabfct]
end

@inline function getprecal!(v::Vector{T}, bsp::AbstractInterpolation{T}, decf::T) where {T}
    return v .= getprecal(bsp, decf)
end

@inline function getprecal!(
    v::Vector{T},
    bsp::Vector{I},
    decf::T,
) where {T,I<:AbstractInterpolation{T}}
    return getprecal!(v, bsp[1], decf)
end

@inline function getprecal!(
    v::Array{T,N},
    bsp::Vector{I},
    decf::NTuple{N,T},
) where {T,N,I<:AbstractInterpolation{T}}
    return v .= dotprod(ntuple(x -> getprecal(bsp[x], decf[x]), N))
end

"""
Return the stencil weights for [`ExtrapolationBC`](@ref) interpolation.

The shift is `decint + decfloat` grid cells, with `decint` an integer and
`decfloat` the fractional part. Near the ends of the line the stencil is
shifted inside the data, so several different sets of weights are needed; this
function returns the vector of all of them (`order + 1` vectors of
`order + 1` weights).
"""
function get_allprecal(
    interp::AbstractInterpolation{T,ExtrapolationBC,order},
    decint::Int,
    decfloat::T,
) where {T,order}
    origin = -div(order, 2)
    indbeg = origin + decint
    indend = indbeg + order
    return [getprecal(interp, decfloat + i) for i = indbeg:indend]
end

"""
Internal helper used to sort complex displacements by integer part and then by
fractional part when inverting a displacement field.
"""
struct ValInv{T}
    a::Int
    b::T
    val::Complex{T}
    ind::CartesianIndex{2}
    function ValInv(ind::CartesianIndex{2}, val::Complex{T}) where {T}
        return new{T}(round(Int, real(val)), imag(val), val, ind)
    end
end

"""
Order two `ValInv` first by their integer part `a` and then by their
fractional part `b`.
"""
function Base.isless(a::ValInv{T}, b::ValInv{T}) where {T}
    return a.a == b.a ? Base.isless(a.b, b.b) : Base.isless(a.a, b.a)
end

function corind(indref::CartesianIndex{2}, indorigin::CartesianIndex{2}, sz::Tuple{Int,Int})
    return CartesianIndex(modone.(indref.I .- indorigin.I, sz))
end

"""
Shift the data `fi` by an offset and store the result in `fp`, for an
interpolation with [`PeriodicBC`](@ref). The offset is given as an integer
number of cells `decint` and the stencil weights `precal` for the fractional
part, computed beforehand with `getprecal(interp, frac)`. Precomputing them
avoids recomputing the weights when the same shift is applied many times.

After the call `fp[i]` approximates `fi` evaluated at `x_i + (decint + frac) dx`.

# Arguments
- `fp::AbstractVector` : output vector
- `fi::AbstractVector` : input vector
- `decint` : offset in units of dx
- `precal::Vector` : vector of length order+1 precompute with getprecal(interp, dec) (dec is the offset)
- `interp::AbstractInterpolation{T, PeriodicBC, order}` : interpolation implementation
- `tabmod=gettabmod(length(fi))` : precompute for "begin at one" modulo

# Returns :
- No return
"""
function interpolate!(
    fp::AbstractVector{T},
    fi::AbstractVector{T},
    decint::Int,
    precal::Vector{T},
    interp::AbstractInterpolation{T,PeriodicBC,order},
    tabmod = gettabmod(length(fi)),
) where {T,order}
    res = sol(interp, fi)
    origin = -div(order, 2)
    lg = length(fi)
    decal = (origin + decint + 5lg) % lg
    @inbounds for i = 1:lg
        indbeg = i + decal
        indend = indbeg + order
        fp[i] = sum(res[tabmod[indbeg:indend]] .* precal)
    end
    return missing
end

"""
Same as the periodic method with precomputed weights, taking a vector of
interpolations and using its first element.
"""
@inline function interpolate!(
    fp::AbstractVector{T},
    fi::AbstractVector{T},
    decint::Int,
    precal::Vector{T},
    tinterp::Vector{I},
    tabmod = gettabmod(length(fi)),
) where {T,I<:AbstractInterpolation{T,PeriodicBC}}
    return interpolate!(fp, fi, decint, precal, tinterp[1], tabmod[1])
end

"""
Multidimensional version of the periodic interpolation with precomputed
weights: `fp[i] = fi(x_i + decint + frac)` for every index `i`, with a tensor
product stencil.

# Arguments
- `fp::AbstractArray{T,N}` : output array
- `fi::AbstractArray{T,N}` : input array (same size as `fp`)
- `decint::NTuple{N,Int}` : integer part of the shift along each dimension, in grid cells
- `precal::Array{T,N}` : tensor product of the stencil weights of each dimension
- `interp::Vector` : the `N` interpolations (periodic boundary condition), one per dimension
- `tabmod` : precomputed index tables used to wrap indices (see `gettabmod`)
"""
function interpolate!(
    fp::AbstractArray{T,N},
    fi::AbstractArray{T,N},
    decint::NTuple{N,Int},
    precal::Array{T,N},
    interp::Vector{I},
    tabmod = gettabmod.(size(fi)),
) where {T,N,I<:AbstractInterpolation{T,PeriodicBC}}
    res = sol(interp, fi)
    sz = size(fi)
    order = get_order.(interp)
    origin = -div.(order, (2,))
    decall = (5 .* sz .+ origin) .% sz .+ sz
    for ind in CartesianIndices(fi)
        deb_i = decint .+ decall .+ ind.I
        end_i = deb_i .+ order
        fp[ind] = sum(res[ntuple(x -> tabmod[x][deb_i[x]:end_i[x]], N)...] .* precal)
    end
    return missing
end

"""
Shift the data `fi` by an offset and store the result in `fp`, for an
interpolation with [`ExtrapolationBC`](@ref). The offset is given as an integer
number of cells `decint` and the sets of stencil weights `allprecal` computed
beforehand with `get_allprecal(interp, decint, frac)`. Several sets are needed
because the stencil is moved inside the line near its ends, where the result
is extrapolated.

# Arguments
- `fp::AbstractVector` : output vector
- `fi::AbstractVector` : input vector
- `decint` : offset in units of dx
- `allprecal::Vector{Vector{T}}` : vector of vector of length order+1 precompute with getprecal(interp, dec) (dec is the offset)
- `interp::AbstractInterpolation{T, ExtrapolationBC, order}` : interpolation implementation
- `tabmod=gettabmod(length(fi))` : precompute for "begin at one" modulo

# Returns :
- No return
"""
function interpolate!(
    fp::AbstractVector{T},
    fi::AbstractVector{T},
    decint::Int,
    allprecal::Vector{Vector{T}},
    interp::AbstractInterpolation{T,ExtrapolationBC,order},
    tabmod = gettabmod(length(fi));
) where {T,order}
    res = sol(interp, fi)
    origin = -div(order, 2)
    lg = length(fi)
    lgp = length(allprecal)
    borne1 = -decint - origin
    borne2 = lg - decint + origin - 1
    for i = 1:borne1
        indbeg = 1
        indend = order + 1
        ind = i
        #        @show 1, i, ind, indbeg, indend, decint, lgp
        fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
    end
    for i = (borne1+1):borne2
        indbeg = i - borne1
        indend = indbeg + order
        ind = borne1 + 1
        #        @show 2, i, ind, indbeg, indend, decint, lgp
        fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
    end
    for i = (borne2+1):lg
        indbeg = lg - order
        indend = lg
        ind = lgp - (lg - i)
        #        @show 3, i, ind, indbeg, indend, decint, lgp
        fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
    end
    return missing
end

"""
Shift the data `fi` by `dec` grid cells and store the result in `fp`:
`fp[i]` approximates the function sampled by `fi` at `x_i + dec * dx`.
This is the simplest entry point: the integer and fractional parts of `dec`
are separated and the stencil weights are computed for you. The boundary
condition is the one of `interp` ([`PeriodicBC`](@ref) or
[`ExtrapolationBC`](@ref)). To apply the same shift many times, use the
variants taking precomputed weights.

# Arguments
- `fp` : output vector of length n
- `fi` : input vector of length n
- `dec` : offset in units of `dx` (grid cells), positive or negative, may be larger than 1
- `interp::AbstractInterpolation` : interpolation implementation

# Returns
Nothing useful; the result is stored in `fp`.

# Example
```julia
n = 64
x = range(0, 1; length = n + 1)[1:n]
fi = sin.(2π .* x)
fp = similar(fi)
interpolate!(fp, fi, 0.3, Lagrange(5))  # fp ≈ sin.(2π .* (x .+ 0.3 / n))
```
"""
function interpolate!(
    fp,
    fi,
    dec,
    interp::AbstractInterpolation{T,BC,order},
) where {T,BC,order}
    decint = convert(Int, floor(dec))
    decfloat = dec - decint
    if BC === PeriodicBC
        return interpolate!(fp, fi, decint, getprecal(interp, decfloat), interp)
    else
        return interpolate!(fp, fi, decint, get_allprecal(interp, decint, decfloat), interp)
    end
end

"""
Cache of the stencil weights of `N` interpolations (one per dimension). The
weights depend only on the fractional part of the shift, so they are
recomputed only when the shift changes between two calls to `getprecal`.
In multithreaded loops, use one cache per thread.

# Fields
- `interps` : interpolations used along each dimension
- `cache_alpha` : last shift for which the weights were computed
- `cache_int` : integer part of `cache_alpha`
- `precal` : tensor product of the stencil weights of each dimension
"""
mutable struct CachePrecal{T,N,I}
    interps::Vector{I}
    cache_alpha::Union{NTuple{N,T},T}
    cache_int::Union{NTuple{N,Int},Int}
    precal::Array{T,N}
    function CachePrecal(
        interps::Vector{I},
        ::T,
    ) where {T<:Real,I<:AbstractInterpolation{T}}
        N = length(interps)
        cache_alpha = (N == 1) ? zero(T) : ntuple(x -> zero(T), N)
        cache_int = (N == 1) ? 0 : ntuple(x -> 0, N)
        sz = totuple(get_order.(interps) .+ 1)
        precal = zeros(T, sz)
        getprecal!(precal, interps, cache_alpha)
        return new{T,N,I}(interps, cache_alpha, cache_int, precal)
    end
    function CachePrecal(
        interps::Vector{I},
        x::Complex{T},
    ) where {T<:Real,I<:AbstractInterpolation{T}}
        return CachePrecal(interps, real(x))
    end
    function CachePrecal(
        interps::Vector{I},
        x::OpTuple{N,T},
    ) where {N,T<:Real,I<:AbstractInterpolation{T}}
        return CachePrecal(interps, x.v[1])
    end
end

"""
Return `(integer_part, weights)` for the shift `alpha` (one value per
dimension), updating the cache only if `alpha` changed.
"""
@inline function getprecal(self::CachePrecal{T,N}, alpha::NTuple{N,T}) where {T,N}
    if alpha != self.cache_alpha
        self.cache_alpha = alpha
        self.cache_int = Int.(floor.(alpha))
        decfloat = alpha .- self.cache_int
        getprecal!(self.precal, self.interps, decfloat)
    end
    return self.cache_int, self.precal
end

"""
Same as above for a two-dimensional shift given as a complex number
`alpha_1 + i alpha_2`.
"""
@inline function getprecal(self::CachePrecal{T,2}, alpha::Complex{T}) where {T<:Real}
    return getprecal(self, reim(alpha))
end

"""
Same as above for a shift given as an `OpTuple`.
"""
@inline function getprecal(self::CachePrecal{T,N}, alpha::OpTuple{N,T}) where {N,T<:Real}
    return getprecal(self, alpha.v)
end

"""
Return `(integer_part, weights)` for the scalar shift `alpha`, updating the
cache only if `alpha` changed.
"""
@inline function getprecal(self::CachePrecal{T,1}, alpha::T) where {T}
    if alpha != self.cache_alpha
        self.cache_alpha = alpha
        self.cache_int = Int(floor(alpha))
        decfloat = alpha - self.cache_int
        getprecal!(self.precal, self.interps[1], decfloat)
    end
    return self.cache_int, self.precal
end

"""
Same as above for a one-dimensional shift given as a 1-tuple.
"""
@inline function getprecal(self::CachePrecal{T,1}, alpha::Tuple{T}) where {T}
    return getprecal(self, alpha[1])
end

"""
Multidimensional periodic interpolation where the shift depends on the grid
point: `dec(ind)` returns the shift (a tuple, a complex number or an
`OpTuple`, in grid cells) to apply at the `CartesianIndex` `ind`.

# Arguments
- `fp` : output array
- `fi` : input array (same size as `fp`)
- `dec::Function` : maps a `CartesianIndex` to the shift at this point
- `interp_t::AbstractVector` : the `N` interpolations, one per dimension

# Keyword arguments
- `tabmod` : precomputed index tables used to wrap indices
- `cache::CachePrecal` : cache of the stencil weights
"""
function interpolate!(
    fp::Union{AbstractArray{T,N},AbstractArray{Complex{T},N}},
    fi::Union{AbstractArray{T,N},AbstractArray{Complex{T},N}},
    dec::Function,
    interp_t::AbstractVector{I};
    tabmod::NTuple{N,Vector{Int}} = gettabmod.(size(fi)),
    cache::CachePrecal{T,N} = CachePrecal(interp_t, one(eltype(fp))),
) where {T,N,I<:AbstractInterpolation{T}}
    N == length(interp_t) || thrown(
        ArgumentError(
            "The number of Interpolation $(length(interp_t)) is different of N=$N",
        ),
    )
    sz = size(fp)
    #    @show "interpolate!", sz
    res = sol(interp_t, fi)

    order = get_order.(interp_t)
    origin = -div.(order, (2,))
    decall = (5 .* sz .+ origin) .% sz .+ sz

    for ind in CartesianIndices(sz)
        dint, tab = getprecal(cache, dec(ind))
        deb_i = dint .+ decall .+ ind.I
        end_i = deb_i + order
        fp[ind] = sum(res[ntuple(x -> tabmod[x][deb_i[x]:end_i[x]], N)...] .* tab)
    end
    #    @show "std2d", diffmax, decminmin, decmaxmax
end

"""
Internal helper: periodic interpolation of `fi` with the displacement field
`bufc` (one array per dimension), the displacements being taken modulo the
mesh length `lgmesh`.
"""
function interpolatemod!(
    fp::AbstractArray{T,N},
    fi::AbstractArray{T,N},
    bufc::Tuple{Array{T,N},Array{T,N}},
    interp_t::AbstractVector{I},
    lgmesh::Union{T,UniformMesh{T}},
) where {T,N,I<:AbstractInterpolation{T,PeriodicBC}}
    (N == 2 && N == length(interp_t)) || thrown(
        ArgumentError(
            "The number of Interpolation $(length(interp_t)) is different of N=$N",
        ),
    )

    sz = size(fp)

    res = sol(interp_t, fi)

    extr = extrema.(bufc)
    decbegin = ntuple(x -> -Int(floor(extr[x][1])), N)
    decend = ntuple(x -> Int(ceil(extr[x][2])), N)

    order = totuple(get_order.(interp_t))
    origin = 0 .- div.(order, (2,))

    dec_b = decbegin .- origin
    dec_e = decend .+ order .+ origin
    fiext = getextarray(res, dec_b, dec_e)

    decall = decbegin

    traitmodbegin!(lgmesh, fiext)

    cache::CachePrecal{T,N} = CachePrecal(interp_t, one(eltype(fp)))
    diffmax = 0
    decminmin = ntuple(x -> Inf, N)
    decmaxmax = ntuple(x -> -Inf, N)
    decabs = ntuple(x -> 0, N)

    for ind in CartesianIndices(sz)
        dint, tab = getprecal(cache, (bufc[1][ind], bufc[2][ind]))
        decabs = max.(decabs, abs.((bufc[1][ind], bufc[2][ind])))
        decmin, decmax = extrema(tab)
        decminmin = min.(decminmin, dint .+ decmin)
        decmaxmax = max.(decmaxmax, dint .+ decmax)

        deb_i = dint .+ decall .+ ind.I
        end_i = deb_i .+ order
        vmin, vmax = extrema(fiext[ntuple(x -> deb_i[x]:end_i[x], N)...])
        diff = vmax - vmin
        diffmax = max(diff, diffmax)
        fp[ind] = sum(fiext[ntuple(x -> deb_i[x]:end_i[x], N)...] .* tab)
    end

    traitmodend!(lgmesh, fp)

    return missing

end

"""
Internal helper: compute by fixed-point iterations the inverse displacement
field of `dec` (one array per dimension), using the periodic interpolation `interp`.
"""
function getinverse(
    dec::NTuple{N,Array{T,N}},
    interp::Vector{I},
) where {T,N,I<:AbstractInterpolation{T,PeriodicBC}}
    sz = size(dec[1])
    res = ntuple(x -> [T(ind.I[x] - 1) for ind in CartesianIndices(sz)], N)
    decinv = ntuple(x -> [T(ind.I[x] - 1) for ind in CartesianIndices(sz)], N)
    buf = zeros(T, sz)
    res2 = ntuple(x -> zeros(T, sz), N)
    res3 = ntuple(x -> zeros(T, sz), N)
    res4 = ntuple(x -> zeros(T, sz), N)
    decfmr = ntuple(x -> zeros(T, sz), N)

    sz_2 = div.(sz, 2)

    for i = 1:N
        interpolatemod!(res2[i], res[i], dec, interp, T(sz[i]))
        interpolatemod!(decfmr[i], -dec[i], dec, interp, T(sz[i]))
        decfmr[i] .= mod.(decfmr[i] .+ sz_2[i], sz[i]) .- sz_2[i]

        res4[i] .= res2[i]

        res3[i] .= mod.(res2[i] .- res[i] .- dec[i] .+ sz_2[i], sz[i]) .- sz_2[i]
    end

    borne = log(eps(T) * prod(sz))

    for z = 1:2
        ind = 1
        note = Inf
        while note > borne
            for i = 1:N
                interpolatemod!(res3[i], res2[i], decfmr, interp, T(sz[i]))
                interpolatemod!(buf, decinv[i], decfmr, interp, T(sz[i]))
                decinv[i] .= buf
            end
            for i = 1:N
                decfmr[i] .= mod.(res[i] .- res3[i] .+ sz_2[i], sz[i]) .- sz_2[i]
                res2[i] .= res3[i]
            end
            note = log(2, norm(decfmr))
            if note > 1
                for i = 1:N
                    decfmr[i] .= decfmr[i] / note
                end
            end
            @show ind, note, borne
            ind += 1
        end
        for i = 1:N
            decinv[i] .= mod.(decinv[i] .- res[i] .+ sz_2[i], sz[i]) .- sz_2[i]
            if z == 1
                decfmr[i] .= decinv[i]
                decinv[i] .= res[i]
                res2[i] .= res4[i]
            end
        end
    end
    return decinv
end


"""
Multidimensional periodic interpolation with the shift of each grid point
stored in the array `bufdec`: `fp[i] = fi(x_i + bufdec[i])`. Returns `true`.

# Arguments
- `fp` : output array
- `fi` : input array (same size as `fp`)
- `bufdec` : array of shifts (in grid cells), one per grid point
- `interp_t::AbstractVector` : the `N` interpolations, one per dimension

# Keyword arguments
- `tabmod` : precomputed index tables used to wrap indices
- `mpid` : unused, kept for compatibility
- `t_split::Tuple` : if given together with `cachethreads`, split of the indices between threads
- `cachethreads::Vector{CachePrecal}` : one weight cache per thread; when given the loop is multithreaded
"""
function interpolate!(
    fp::Union{AbstractArray{T,N},AbstractArray{OpTuple{N,T},N},AbstractArray{Complex{T},N}},
    fi::Union{AbstractArray{T,N},AbstractArray{OpTuple{N,T},N},AbstractArray{Complex{T},N}},
    bufdec::Union{AbstractArray{OpTuple{N,T},N},AbstractArray{Complex{T},N}},
    interp_t::AbstractVector{I};
    tabmod::NTuple{N,Vector{Int}} = gettabmod.(size(fi)),
    mpid = missing,
    t_split::Union{Tuple,Missing} = missing,
    cachethreads::Union{Vector{CachePrecal{T}},Missing} = missing,
    #    itr::AbstractArray=CartesianIndices(fi)
) where {T,N,I<:AbstractInterpolation{T}}
    N == length(interp_t) || thrown(
        ArgumentError(
            "The number of Interpolation $(length(interp_t)) is different of N=$N",
        ),
    )

    isthreads = !ismissing(cachethreads)
    sz = size(fp)
    #    @show "interpolate!", sz
    res = sol(interp_t, fi)

    order = get_order.(interp_t)
    origin = -div.(order, (2,))
    decall = (5 .* sz .+ origin) .% sz .+ sz

    function fct(
        ind::CartesianIndex{N2},
        cache::CachePrecal{T2,N2,I2},
    ) where {N2,T2,I2<:AbstractInterpolation{T2}}
        dint, tab = getprecal(cache, bufdec[ind])
        deb_i = dint .+ decall .+ ind.I
        end_i = deb_i + order
        return fp[ind] = sum(res[ntuple(x -> tabmod[x][deb_i[x]:end_i[x]], N2)...] .* tab)
    end
    ci = CartesianIndices(sz)
    if isthreads
        if ismissing(t_split)
            #            println("SimpleThreads")
            @threads for ind in ci
                local cache = cachethreads[Threads.threadid()]
                fct(ind, cache)
            end
        else
            #            println("SplitThreads")
            @threads for it in t_split
                local cache = cachethreads[Threads.threadid()]
                for ind in ci[it]
                    fct(ind, cache)
                end
            end
        end
    else
        local cache = CachePrecal(interp_t, zero(T))
        for ind in ci
            fct(ind, cache)
        end
    end

    return true
end

"""
Apply `nb` successive interpolations of `from` with the displacement
field `from` itself (`to` receives the result). If `nb < 1`, `to` is a copy of
`from`.
"""
function autointerp!(
    to::Array{OpTuple{N,T},N},
    from::Array{OpTuple{N,T},N},
    nb::Int,
    interp_t::AbstractVector{I};
    mpid = missing,
    t_split::Union{Tuple,Missing} = missing,
    cachethreads::Union{Vector{CachePrecal{T}},Missing} = missing,
) where {N,T,I<:AbstractInterpolation}
    if nb < 1
        to .= from
    end
    fmr = copy(from)
    for i = 1:nb
        interpolate!(
            to,
            from,
            fmr,
            interp_t;
            mpid = mpid,
            t_split = t_split,
            cachethreads = cachethreads,
        )
        if i != nb
            fmr .= to
        end
    end
end


"""
Interpolate in place the `nb` last arrays of `t_buf` with the displacement
field `bufdec`.
"""
function interpbufc!(
    t_buf::Vector{Array{OpTuple{N,T},N}},
    bufdec::Array{OpTuple{N,T},N},
    interp_t::AbstractVector{I},
    nb::Int = length(t_buf);
    mpid = missing,
    t_split::Union{Tuple,Missing} = missing,
    cachethreads::Union{Vector{CachePrecal{T}},Missing} = missing,
) where {N,T,I<:AbstractInterpolation{T}}

    for i = 0:(nb-1)
        buf = t_buf[end-i]
        interpolate!(
            buf,
            copy(buf),
            bufdec,
            interp_t;
            mpid = mpid,
            t_split = t_split,
            cachethreads = cachethreads,
        )
    end

end

interpolate!(fp, fi, dec::Real, interp::AbstractShiftInterpolation) =
    interpolate!(fp, interp, fi, dec)
