using LinearAlgebra
import Base: isless
@enum EdgeType CircEdge = 1 InsideEdge = 2

"""
    AbstractInterpolation{T, edge, order, nd}

Abstract supertype for all interpolation type

# Type parameters
- `T` : type of number on witch interpolation works
- `edge` : type of edge treatment
- `order` : order of interpolation

# Implementation constraint
- `tabfct::Vector` : this attribut must be on the implementation, it is a table of function of size order+1
"""
abstract type AbstractInterpolation{T,edge,order} end

"""
    get_order(_::AbstractInterpolation{T, edge, order}) where{T, edge, order}
Return the order of interpolation implementation       
"""
get_order(_::AbstractInterpolation{T,edge,order}) where {T,edge,order} = order

"""
    sol(_::AbstractInterpolation, line::AbstractVector)

Interface method to transform the treated line, by default this method does nothing

# Arguments :
- `_::AbstractInterpolation` : interpolation implementation
- `line::AbstractVector` : line to transform

# Return :
The transformed line

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

function get_allprecal(
    interp::AbstractInterpolation{T,InsideEdge,order},
    decint::Int,
    decfloat::T,
) where {T,order}
    origin = -div(order, 2)
    indbeg = origin + decint
    indend = indbeg + order
    return [getprecal(interp, decfloat + i) for i = indbeg:indend]
end
struct ValInv{T}
    a::Int
    b::T
    val::Complex{T}
    ind::CartesianIndex{2}
    function ValInv(ind::CartesianIndex{2}, val::Complex{T}) where {T}
        return new{T}(round(Int, real(val)), imag(val), val, ind)
    end
end
function Base.isless(a::ValInv{T}, b::ValInv{T}) where {T}
    return a.a == b.a ? Base.isless(a.b, b.b) : Base.isless(a.a, b.a)
end

function corind(indref::CartesianIndex{2}, indorigin::CartesianIndex{2}, sz::Tuple{Int,Int})
    return CartesianIndex(modone.(indref.I .- indorigin.I, sz))
end

"""
    interpolate!( fp::AbstractVector{T}, 
        fi::AbstractVector{T},
        decint::Int, 
        precal::Vector{T}, 
        interp::AbstractInterpolation{T, CircEdge, order},
        self::AdvectionData{T},
        tabmod=gettabmod(size(fi)) ) where {T, order}

apply an offset to the function fi interpolate by interp struct, the result is in fp vector,
decint and precal are precompute with getprecal method, the TypeEdge is CircEdge

# Arguments
- `fp::AbstractVector` : output vector
- `fi::AbstractVector` : input vector
- `decint` : offset in units of dx
- `precal::Vector` : vector of length order+1 precompute with getprecal(interp, dec) (dec is the offset)
- `interp::AbstractInterpolation{T, CircEdge, order}` : interpolation implementation, note that TypeEdge is CircEdge
- `tabmod=gettabmod(length(fi))` : precompute for "begin at one" modulo

# Returns :
- No return
"""
function interpolate!(
    fp::AbstractVector{T},
    fi::AbstractVector{T},
    decint::Int,
    precal::Vector{T},
    interp::AbstractInterpolation{T,CircEdge,order},
    tabmod = gettabmod(length(fi)),
) where {T,order}
    res = sol(interp, fi)
    # @show typeof(res)
    # @show typeof(tabmod)
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
@inline function interpolate!(
    fp::AbstractVector{T},
    fi::AbstractVector{T},
    decint::Int,
    precal::Vector{T},
    tinterp::Vector{I},
    tabmod = gettabmod(length(fi)),
) where {T,I<:AbstractInterpolation{T,CircEdge}}
    return interpolate!(fp, fi, decint, precal, tinterp[1], tabmod[1])
end

function interpolate!(
    fp::AbstractArray{T,N},
    fi::AbstractArray{T,N},
    decint::NTuple{N,Int},
    precal::Array{T,N},
    interp::Vector{I},
    tabmod = gettabmod.(size(fi)),
) where {T,N,I<:AbstractInterpolation{T,CircEdge}}
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
    interpolate!( 
    fp::AbstractVector{T}, fi::AbstractVector{T}, decint::Int, 
    allprecal::Vector{Vector{T}}, 
    interp::AbstractInterpolation{T, InsideEdge, order},
    tabmod=gettabmod(length(fi))
    ) where {T, order}

apply an offset to the function fi interpolate by interp struct, the result is in fp vector,
decint and precal are precompute with getprecal method, the TypeEdge is InsideEdge, it is a marginal case

# Arguments
- `fp::AbstractVector` : output vector
- `fi::AbstractVector` : input vector
- `decint` : offset in units of dx
- `allprecal::Vector{Vector{T}}` : vector of vector of length order+1 precompute with getprecal(interp, dec) (dec is the offset)
- `interp::AbstractInterpolation{T, InsideEdge, order}` : interpolation implementation, note that TypeEdge is CircEdge
- `tabmod=gettabmod(length(fi))` : precompute for "begin at one" modulo

# Returns :
- No return
"""
function interpolate!(
    fp::AbstractVector{T},
    fi::AbstractVector{T},
    decint::Int,
    allprecal::Vector{Vector{T}},
    interp::AbstractInterpolation{T,InsideEdge,order},
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
    interpolate!( fp, fi, dec, interp)

apply the offset dec to the function fi interpolate by interp struct, the result is in fp Vector

# Arguments
- `fp` : output vector of length n
- `fi` : input vector of length n
- `dec` : offset in units of dx
- `interp::AbstractInterpolation` : interpolation implementation

# Returns :
- No return
"""
function interpolate!(
    fp,
    fi,
    dec,
    interp::AbstractInterpolation{T,edge,order},
) where {T,edge,order}
    decint = convert(Int, floor(dec))
    decfloat = dec - decint
    if edge == CircEdge
        return interpolate!(fp, fi, decint, getprecal(interp, decfloat), interp)
    else
        return interpolate!(fp, fi, decint, get_allprecal(interp, decint, decfloat), interp)
    end
end

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
@inline function getprecal(self::CachePrecal{T,N}, alpha::NTuple{N,T}) where {T,N}
    if alpha != self.cache_alpha
        self.cache_alpha = alpha
        self.cache_int = Int.(floor.(alpha))
        decfloat = alpha .- self.cache_int
        getprecal!(self.precal, self.interps, decfloat)
    end
    return self.cache_int, self.precal
end
@inline function getprecal(self::CachePrecal{T,2}, alpha::Complex{T}) where {T<:Real}
    return getprecal(self, reim(alpha))
end
@inline function getprecal(self::CachePrecal{T,N}, alpha::OpTuple{N,T}) where {N,T<:Real}
    return getprecal(self, alpha.v)
end
@inline function getprecal(self::CachePrecal{T,1}, alpha::T) where {T}
    if alpha != self.cache_alpha
        self.cache_alpha = alpha
        self.cache_int = Int(floor(alpha))
        decfloat = alpha - self.cache_int
        getprecal!(self.precal, self.interps[1], decfloat)
    end
    return self.cache_int, self.precal
end
@inline function getprecal(self::CachePrecal{T,1}, alpha::Tuple{T}) where {T}
    return getprecal(self, alpha[1])
end

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

function interpolatemod!(
    fp::AbstractArray{T,N},
    fi::AbstractArray{T,N},
    bufc::Tuple{Array{T,N},Array{T,N}},
    interp_t::AbstractVector{I},
    lgmesh::Union{T,UniformMesh{T}},
) where {T,N,I<:AbstractInterpolation{T,CircEdge}}
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

    # @show size(fi), decbegin, decend
    # @show order, origin, dec_b, dec_e, size(fiext), decall
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
        # if ind.I[1] <= 10 && ind.I[2] <= 10
        #     indcor = ind +CartesianIndex(decall .- origin)
        #     @show decall, dint, origin
        #     @show ind, fp[ind], indcor, fiext[indcor]
        # end
    end
    # @show diffmax, decminmin, decmaxmax, decabs
    # @show fp[1:10,1:10]
    traitmodend!(lgmesh, fp)
    # @show fp[1:10,1:10]
    return missing
end

function getinverse(
    dec::NTuple{N,Array{T,N}},
    interp::Vector{I},
) where {T,N,I<:AbstractInterpolation{T,CircEdge}}
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

    @show norm(res3)

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


function interpolate!(
    fp::Union{AbstractArray{T,N},AbstractArray{OpTuple{N,T},N},AbstractArray{Complex{T},N}},
    fi::Union{AbstractArray{T,N},AbstractArray{OpTuple{N,T},N},AbstractArray{Complex{T},N}},
    bufdec::Union{AbstractArray{OpTuple{N,T},N},AbstractArray{Complex{T},N}},
    interp_t::AbstractVector{I};
    tabmod::NTuple{N,Vector{Int}} = gettabmod.(size(fi))
) where {T,N,I<:AbstractInterpolation{T}}
    N == length(interp_t) || thrown(
        ArgumentError(
            "The number of Interpolation $(length(interp_t)) is different of N=$N",
        ),
    )

    sz = size(fp)
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
    local cache = CachePrecal(interp_t, zero(T))
    for ind in CartesianIndices(sz)
        fct(ind, cache)
    end

    return true
end

