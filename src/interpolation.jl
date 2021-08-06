using LinearAlgebra
import Base: isless, zero, iterate, +, -, *
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

struct OpTuple{N,T}
    v::NTuple{N,T}
end
Base.show(io::IO, ot::OpTuple) = print(io, ot.v)

(+)(v1::OpTuple,v2::OpTuple) = OpTuple(v1.v .+ v2.v)
(-)(v1::OpTuple,v2::OpTuple) = OpTuple(v1.v .- v2.v)
(*)(a::Number, v::OpTuple) = OpTuple( a .* v.v)
(*)(v::OpTuple, a::Number) = OpTuple( v.v .* a)
Base.zero(::Type{OpTuple{N,T}}) where{N,T}=OpTuple(ntuple(x->zero(T),N))
Base.iterate(v::OpTuple)=Base.iterate(v.v)
Base.iterate(v::OpTuple, state)=Base.iterate(v.v, state)
Base.length(::OpTuple{N}) where N=N

function sol!(
    Y::AbstractArray{T,N},
    interp_t::AbstractVector{I},
    b::AbstractArray{T,N},
) where {T,N,I<:AbstractInterpolation{T}}
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
            # order = get_order(interp)
            # t = getbspline(big(order), 0).(1:order)
            # A = topl(sz[1], t, true)

            #            @show i, sz
            bufin = (i == 1) ? b : permutedims(bufout, perm)
            bufout = zeros(T, sz)
            # diffmax = 0
            for c in CartesianIndices(sz[2:end])
                bo = view(bufout, :, c)
                bi = view(bufin, :, c)
                #                bbi = copy(bi)
                sol!(bo, interp, bi)
                # diff = norm(bbi - A*bo)
                # diffmax = max(diff,diffmax)
            end
            #            @show diffmax
            perm = p
        end
        sz = sz[p]
    end
    permutedims!(Y, bufout, perm)
    return Y
end
function sol(
    interp_t::AbstractVector{I},
    b::Union{AbstractArray{T,N},AbstractArray{Complex{T},N},AbstractArray{OpTuple{N,T},N} }
) where {T<: Real,N,I<:AbstractInterpolation{T}}
    if all(issolidentity.(interp_t))
        return b
    else
        return sol!(zeros(T, size(b)), interp_t, b)
    end
end

@inline function getprecal(interp::AbstractInterpolation{T}, decf::T) where {T}
    return [T(fct(decf)) for fct in interp.tabfct]
end

@inline function getprecal!(v::Vector{T}, bsp::AbstractInterpolation{T}, decf::T) where {T}
    v .= getprecal(bsp, decf)
end
@inline function getprecal!(
    v::Vector{T},
    bsp::Vector{I},
    decf::T,
) where {T,I<:AbstractInterpolation{T}}
    getprecal!(v, bsp[1], decf)
end
@inline function getprecal!(
    v::Array{T,N},
    bsp::Vector{I},
    decf::NTuple{N,T},
) where {T,N,I<:AbstractInterpolation{T}}
    v .= dotprod(ntuple(x -> getprecal(bsp[x], decf[x]), N))
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
#aroundindices()=CartesianIndex.([(1,0),(1,1),(0,1),(-1,1),(-1,0),(-1,-1),(0,-1),(1,-1)])
#aroundcmplx()=[  ]

function corind(indref::CartesianIndex{2}, indorigin::CartesianIndex{2}, sz::Tuple{Int,Int})
    return CartesianIndex(modone.(indref.I .- indorigin.I, sz))
end
cmplxtoindex(v) = CartesianIndex(Int(real(v)), Int(imag(v)))
indextocmplx(ind) = ind.I[1] + ind.I[2] * im

function roundcomplex(v::Complex, sz::Tuple{Int,Int})
    b = div.(sz, (2,))
    return mod(real(v)+b[1],sz[1])-b[1] +im*(mod(imag(v)+b[2],sz[2])-b[2])
end

function calinverse!(
    out::NTuple{N,Array{T,2}},
    tabin::NTuple{N,Array{T,2}},
    order::Int=4
) where {T,N,I<:AbstractInterpolation{T,CircEdge}}

    @assert order % 2 == 0 "order=$order must be even"

    lginterp = order + 1
    origin = div(order, 2)
    sz = size(tabin[1])
    indorigin = CartesianIndex(origin, origin)
    indmiddle = CartesianIndex(origin + 1, origin + 1)
    indlg = CartesianIndex(lginterp,lginterp)

    tab1 = getextarray(tabin[1], (origin, origin), (origin, origin))
    tab2 = getextarray(tabin[2], (origin, origin), (origin, origin))

    largesz = size(tab1)

    traitmodbegin!(T(sz[1]), tab1)
    traitmodbegin!(T(sz[2]), tab2)

    minreal, maxreal = extrema(tab1)
    minimag, maximag = extrema(tab2)
    imagit(v::Int) = minimag+mod(v - minimag, sz[2]):sz[2]:maximag-mod(maximag - v, sz[2])

    @show minreal, maxreal, minimag, maximag

    minintreal = round(Int, minreal)
    maxintreal = round(Int, maxreal)
    realit(v::Int) =
        minintreal+mod(v - minintreal, sz[1]):sz[1]:maxintreal-mod(maxintreal - v, sz[1])

    taball = tab1 + im * tab2
    taballbig = big.(taball)

    # tabinvall est le tableau qui donne à partir de quel indice on va interpoler
    tabvalind = Vector{ValInv}(undef, largesz[1] * largesz[2])
    for ind in CartesianIndices(taball)
        tabvalind[(ind.I[2]-1)*largesz[1]+ind.I[1]] = ValInv(ind, taball[ind])
    end
    sort!(tabvalind)

    rreal = extrema([p.ind.I[1] for p in tabvalind])
    rimag = extrema([p.ind.I[2] for p in tabvalind])

    @show rreal, rimag

    tabresult = [CartesianIndex{2}[] for i in CartesianIndices(sz)]

    nullind = CartesianIndex(0, 0)

    srch(val) = searchsorted(tabvalind, ValInv(nullind, val))
    dist(ind::Int, valsearch) = abs(valsearch - tabvalind[ind].val)
    indend = length(tabvalind) + 1
    for ind in CartesianIndices(sz)
        distmin = Inf
        dec = 0
        indmin::Int = 0
        while distmin >= dec + 0.5
            itdec = dec == 0 ? [0] : [-dec, dec]
            for decdec in itdec
                for i in realit(ind.I[1] + decdec), j in imagit(ind.I[2])
                    valsrch = T(i - decdec) + im * j
                    r = srch(valsrch)
                    for v in [r.start, r.stop]
                        if 0 < v < indend
                            d = dist(v, valsrch)
                            if d < distmin
                                indmin = v
                                distmin = d
                            end
                        end
                    end
                end
            end
            dec += 1
        end
        #        push!(tabresult[tabvalind[indmin].ind-indorigin], ind)
        push!(tabresult[corind(tabvalind[indmin].ind, indorigin, sz)], ind)
    end

    for ind in CartesianIndices(sz)
        p = missing
        moyint = missing
        corind = missing
        for indind in tabresult[ind]
             if ismissing(p)
                tab = taballbig[ind.I[1]:ind.I[1]+2origin, ind.I[2]:ind.I[2]+2origin]
                moyint = round(tab[indmiddle])
                tab .-= moyint
                minre,maxre = extrema(real.(tab))
                minim,maxim = extrema(imag.(tab))
                @show minre,maxre,minim,maxim
                p = getpoly(tab)
                corind = indextocmplx(ind-indlg)
            end
            v = p(roundcomplex(indextocmplx(indind) - moyint, sz)) + corind
            @show indind, taballbig[ind+indorigin]
            @show indind, ind, moyint, corind, v
            out[1][indind] = mod(real(v), sz[1])
            out[2][indind] = mod(imag(v), sz[2])
        end
    end
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
    tabmod = gettabmod(length(fi));
    clockobs::AbstractClockObs = NoClockObs(),
) where {T,order}
    res = sol(interp, fi)
    # @show typeof(res)
    # @show typeof(tabmod)
    origin = -div(order, 2)
    lg = length(fi)
    decal = (origin + decint + 5lg) % lg
    clockbegin(clockobs, 3)
    @inbounds for i = 1:lg
        indbeg = i + decal
        indend = indbeg + order
        fp[i] = sum(res[tabmod[indbeg:indend]] .* precal)
    end
    clockend(clockobs, 3)
    missing
end
@inline function interpolate!(
    fp::AbstractVector{T},
    fi::AbstractVector{T},
    decint::Int,
    precal::Vector{T},
    tinterp::Vector{I},
    tabmod = gettabmod(length(fi));
    clockobs::AbstractClockObs = NoClockObs(),
) where {T,I<:AbstractInterpolation{T,CircEdge}}
    return interpolate!(fp, fi, decint, precal, tinterp[1], tabmod[1], clockobs = clockobs)
end
function interpolate!(
    fp::AbstractArray{T,N},
    fi::AbstractArray{T,N},
    decint::NTuple{N,Int},
    precal::Array{T,N},
    interp::Vector{I},
    tabmod = gettabmod.(size(fi));
    clockobs::AbstractClockObs = NoClockObs(),
) where {T,N,I<:AbstractInterpolation{T,CircEdge}}
    res = sol(interp, fi)
    order = get_order.(interp)
    origin = -div.(order, (2,))
    lg = size(fi)
    decal = (origin .+ decint .+ (5 .* lg)) .% lg
    for i in CartesianIndices(lg)
        indbeg = i.I .+ decal
        indend = indbeg .+ order
        fp[i] = sum(res[ntuple(x -> tabmod[x][indbeg[x]:indend[x]], N)...] .* precal)
    end
    missing
end

@inline function interpolate!(
    fp::AbstractArray{T,1},
    fi::AbstractArray{T,1},
    decint::NTuple{1,Int},
    precal::Array{T,1},
    interp::Vector{I},
    tabmod = gettabmod.(size(fi));
    clockobs::AbstractClockObs = NoClockObs(),
) where {T,N,I<:AbstractInterpolation{T,CircEdge}}
    return interpolate!(
        fp,
        fi,
        decint[1],
        precal,
        interp[1],
        self,
        tabmod[1];
        clockobs = clockobs,
    )
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
    clockobs::AbstractClockObs = NoClockObs(),
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
    for i = borne1+1:borne2
        indbeg = i - borne1
        indend = indbeg + order
        ind = borne1 + 1
        #        @show 2, i, ind, indbeg, indend, decint, lgp
        fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
    end
    for i = borne2+1:lg
        indbeg = lg - order
        indend = lg
        ind = lgp - (lg - i)
        #        @show 3, i, ind, indbeg, indend, decint, lgp
        fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
    end
    missing
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

# function interpolate!(
#     fp::AbstractArray{T,2}, 
#     fi::AbstractArray{T,2}, 
#     dec::NTuple{2, Array{T, 2}}, 
#     interp::AbstractInterpolation{T, edge, order, 2},
#     tabmod=gettabmod.(size(fi))
# ) where {T, edge, order}

#     res = sol(interp, fi)

#     tabfct = interp.tabfct
#     decint1 = Int.(floor.(dec[1]))
#     decint2 = Int.(floor.(dec[2]))
#     decfl1 = dec[1] - decint1
#     decfl2 = dec[2] - decint2

#     origin = -div(order,2)
#     dec1=origin+5size(fp,1)
#     dec2=origin+5size(fp,2)

#     for i=1:size(fp,1), j=1:size(fp,2)
#         deb_i = i + decint1[i,j]+dec1
#         end_i = deb_i+order
#         deb_j = j+decint2[i,j]+dec2
#         end_j = deb_j + order
# #            @show size(tabdec1), size(tabdec2)
# #    tab = tabdec1[:, tabmod[1][i+decint1[i]+size(fp,1)]] .* transpose(tabdec2[:, tabmod[2][j+decint2[j]+size(fp,2)]])
# # tab = tabdec1[:, tabmod[2][j+decint1[j]+size(fp,2)]] .* transpose(tabdec2[:, tabmod[1][i+decint1[i]+size(fp,1)]])
#         fl_i = decfl1[i,j]
#         fl_j = decfl2[i,j]
# #        tab = [f(fl_i) for f in tabfct] .* transpose([f(fl_j) for f in tabfct])
#         tab = dotprod(([f(fl_i) for f in tabfct], [f(fl_j) for f in tabfct]))

# #            @show size(tab), size(fi[tabmod[1][deb_i:end_i], tabmod[2][deb_j:end_j]])
#         fp[i,j] = sum( tab .* fi[tabmod[1][deb_i:end_i], tabmod[2][deb_j:end_j]])
#     end
# end
# OK pour n'importe quel nd

mutable struct CachePrecal{T,N,I}
    interps::Vector{I}
    cache_alpha::Union{NTuple{N,T},T}
    cache_int::Union{NTuple{N,Int},Int}
    precal::Array{T,N}
    function CachePrecal(interps::Vector{I}, x::T) where {T <: Real,I<:AbstractInterpolation{T}}
        N = length(interps)
        cache_alpha = (N == 1) ? zero(T) : ntuple(x -> zero(T), N)
        cache_int = (N == 1) ? 0 : ntuple(x -> 0, N)
        sz = totuple(get_order.(interps) .+ 1)
        precal = zeros(T, sz)
        getprecal!(precal, interps, cache_alpha)
        return new{T,N,I}(interps, cache_alpha, cache_int, precal)
    end
    CachePrecal(interps::Vector{I}, x::Complex{T}) where {T <: Real,I<:AbstractInterpolation{T}}=CachePrecal(interps,real(x))
    CachePrecal(interps::Vector{I}, x::OpTuple{N,T}) where {N,T <: Real,I<:AbstractInterpolation{T}}=CachePrecal(interps,x.v[1])
end
@inline function getprecal(self::CachePrecal{T,N}, alpha::NTuple{N,T}) where {T,N}
    if alpha != self.cache_alpha
        self.cache_alpha = alpha
        self.cache_int = Int.(floor.(alpha))
        decfloat = alpha .- self.cache_int
        #        self.cache_precal = get_precal(getinterp(self),decfloat)
        getprecal!(self.precal, self.interps, decfloat)
    end
    return self.cache_int, self.precal
end
@inline getprecal(self::CachePrecal{T,2}, alpha::Complex{T}) where {T<:Real}=getprecal(self,reim(alpha))
@inline getprecal(self::CachePrecal{T,N}, alpha::OpTuple{N,T}) where {N,T<:Real}=getprecal(self,alpha.v)
@inline function getprecal(self::CachePrecal{T,1}, alpha::T) where {T}
    if alpha != self.cache_alpha
        self.cache_alpha = alpha
        self.cache_int = Int(floor(alpha))
        decfloat = alpha - self.cache_int
        #        self.cache_precal = get_precal(getinterp(self),decfloat)
        getprecal!(self.precal, self.interps[1], decfloat)
    end
    return self.cache_int, self.precal
end
@inline getprecal(self::CachePrecal{T,1}, alpha::Tuple{T}) where {T} =
    getprecal(self, alpha[1])

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
    @show "interpolate!", sz
    res = sol(interp_t, fi)

    order = get_order.(interp_t)
    origin = -div.(order, (2,))
    decall = (5 .* sz .+ origin) .% sz .+ sz

    # diffmax = 0
    # decminmin = ntuple(x -> Inf, N)
    # decmaxmax = ntuple(x -> -Inf, N)
    # decminmin = Inf
    # decmaxmax = -Inf

    for ind in CartesianIndices(sz)
        dint, tab = getprecal(cache, dec(ind))
#        decmin, decmax = extrema(tab)
        # decminmin = min.(decminmin, dint .+ decmin)
        # decmaxmax = max.(decmaxmax, dint .+ decmax)
        # decminmin = min.(decminmin, decmin[1])
        # decmaxmax = max.(decmaxmax, decmax[1])
        deb_i = dint .+ decall .+ ind.I
        end_i = deb_i + order
        # vmin, vmax = extrema(res[ntuple(x -> tabmod[x][deb_i[x]:end_i[x]], N)...])
        # diff = vmax - vmin
        # diffmax = max(diff, diffmax)
        fp[ind] = sum(res[ntuple(x -> tabmod[x][deb_i[x]:end_i[x]], N)...] .* tab)
    end
#    @show "std2d", diffmax, decminmin, decmaxmax
end

function interpolate!(
    fp::Union{AbstractArray{T,N},AbstractArray{OpTuple{N,T},N}},
    fi::Union{AbstractArray{T,N},AbstractArray{OpTuple{N,T},N}},
    bufdec::AbstractArray{OpTuple{N,T},N},
    interp_t::AbstractVector{I};
    tabmod::NTuple{N,Vector{Int}} = gettabmod.(size(fi)),
    cache::CachePrecal{T,N} = CachePrecal(interp_t, zero(eltype(fp))),
    mpid::Union{MPIData, Missing}=missing,
    t_split::Union{Tuple, Missing}=missing
#    itr::AbstractArray=CartesianIndices(fi)
) where {T,N,I <: AbstractInterpolation{T}}
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

    itr = ismissing(mpid) ? CartesianIndices(sz) : CartesianIndices(sz)[t_split[mpid.ind]]

    for ind in itr
        dint, tab = getprecal(cache, bufdec[ind])
        deb_i = dint .+ decall .+ ind.I
        end_i = deb_i + order
        fp[ind] = sum(res[ntuple(x -> tabmod[x][deb_i[x]:end_i[x]], N)...] .* tab)
    end
    if !ismissing(mpid)
        mpibroadcast(mpid, t_split, fp)
    end
    true
#    @show "std2d", diffmax, decminmin, decmaxmax
end
function interpolate2!(
    fp::Union{AbstractArray{T,2},AbstractArray{Complex{T},2}},
    fi::Union{AbstractArray{T,2},AbstractArray{Complex{T},2}},
    bufdec::AbstractArray{Complex{T},2},
    interp_t::AbstractVector{I},
) where {T,N,I<:AbstractInterpolation{T}}

    2 == length(interp_t) || thrown(
        ArgumentError(
            "The number of Interpolation $(length(interp_t)) is different of N=$N",
        ),
    )
    sz = size(fp)
    @show "interpolate2!", sz
    res = sol(interp_t, fi)
    supval(x) = x>0 ? x : 0
    extr = extrema.((real.(bufdec), imag(bufdec)))
    decbegin = supval.(0 .-Int.(floor.((extr[1][1], extr[2][1]))))
    decend = supval.(Int.(ceil.((extr[1][2], extr[2][2]))))

    order = totuple(get_order.(interp_t))
    origin = 0 .- div.(order, (2,))

    dec_b = decbegin .- origin
    dec_e = decend .+ order .+ origin
    fiext = getextarray(res, dec_b, dec_e)

    decall = decbegin
    cache::CachePrecal{T,2} = CachePrecal(interp_t, one(eltype(fp)))

    for ind in CartesianIndices(sz)
        dint, tab = getprecal(cache, bufdec[ind])
        deb_i = dint .+ decall .+ ind.I
        end_i = deb_i .+ order
        fp[ind] = sum(fiext[deb_i[1]:end_i[1], deb_i[2]:end_i[2]] .* tab)
    end
#    @show "std2d", diffmax, decminmin, decmaxmax
end
function interpolatemod!(
    fp::AbstractArray{T,N},
    fi::AbstractArray{T,N},
    dec::Function,
    interp_t::AbstractVector{I},
    lgmesh::Union{T,UniformMesh{T}},
    decbegin::NTuple{N,Int},
    decend::NTuple{N,Int},
) where {T,N,I<:AbstractInterpolation{T,CircEdge}}

    N == length(interp_t) || thrown(
        ArgumentError(
            "The number of Interpolation $(length(interp_t)) is different of N=$N",
        ),
    )

    sz = size(fp)

    res = sol(interp_t, fi)

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
        dint, tab = getprecal(cache, dec(ind))
        decabs = max.(decabs, abs.(dec(ind)))
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
    missing
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
        dint, tab = getprecal(cache, (bufc[1][ind],bufc[2][ind]))
        decabs = max.(decabs, abs.((bufc[1][ind],bufc[2][ind])))
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
    missing
end
# function interpolate!(
#     fp::AbstractVector{T},
#     fi::AbstractVector{T},
#     dec::Function,
#     interp::AbstractInterpolation{T},
#     tabmod::Vector{Int},
#     cache::CachePrecal{T,1},
# ) where {T}
#     lg = length(fi)

#     res::Vector{T} = sol(interp, fi)

#     order = get_order(interp)
#     origin = -div(order,2)
#     decall = 5lg + origin

#     for ind in CartesianIndices((lg,))
#         dint, decfl = getprecal(self, dec(ind)[1])
# #        @show decfl, dint
#         deb_i =  (dint + decall)%lg + ind.I[1]
#         end_i = deb_i + order
#         fp[ind] = sum(res[tabmod[deb_i:end_i]] .* tab)
#     end
#     fp
# end
# function interpolate!(
#     fp::AbstractVector{T},
#     fi::AbstractVector{T},
#     dec::Function,
#     interp_t::AbstractVector{I},
#     tabmod::NTuple{1,Vector{Int}},
#     cache::CachePrecal{T,1},
# ) where {T,I<:AbstractInterpolation{T}}
#     interpolate!(fp, fi, dec, interp_t[1], tabmod[1], cache)
# end
function getinverse(dec::NTuple{N,Array{T,N}}, interp::Vector{I}) where {T,N,I<:AbstractInterpolation{T,CircEdge}}
    sz = size(dec[1])
    res = ntuple(x -> [T(ind.I[x]-1) for ind in CartesianIndices(sz)], N)
    decinv = ntuple(x -> [T(ind.I[x]-1) for ind in CartesianIndices(sz)], N)
    buf = zeros(T,sz)
    res2 = ntuple( x -> zeros(T,sz), N)
    res3 = ntuple( x -> zeros(T,sz), N)
    res4 = ntuple( x -> zeros(T,sz), N)
    decfmr = ntuple( x -> zeros(T,sz), N)

    sz_2 = div.(sz,2)

    for i=1:N
        interpolatemod!(res2[i],res[i],dec,interp, T(sz[i]))
        interpolatemod!(decfmr[i], -dec[i], dec, interp, T(sz[i]))
        decfmr[i] .= mod.(decfmr[i] .+ sz_2[i], sz[i]) .- sz_2[i]

        res4[i] .= res2[i]

        res3[i] .= mod.( res2[i] .- res[i] .- dec[i] .+ sz_2[i], sz[i]) .- sz_2[i]

    end

    @show norm(res3)
    

 
    borne = log(eps(T)*prod(sz))

    for z=1:2
        ind = 1
        note = Inf
        while note > borne
            for i=1:N
                interpolatemod!(res3[i], res2[i], decfmr, interp, T(sz[i]))
                interpolatemod!(buf, decinv[i], decfmr, interp, T(sz[i]))
                decinv[i] .= buf
            end
            for i=1:N
                decfmr[i] .= mod.(res[i] .- res3[i] .+ sz_2[i] , sz[i]) .- sz_2[i]
                res2[i] .= res3[i]
            end
            note = log(2, norm(decfmr))
            if note > 1
                for i=1:N
                    decfmr[i] .= decfmr[i]/note
                end
            end
            @show ind,note,borne
            ind += 1
        end
        for i=1:N
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
 
function autointerp!(
    to::Array{OpTuple{N,T},N},
    from::Array{OpTuple{N,T},N},
    nb::Int,
    interp_t::AbstractVector{I};
     mpid::Union{MPIData, Missing}=missing,
    t_split::Union{Tuple, Missing}=missing
)   where{N,T, I<:AbstractInterpolation}
    if nb < 1
        to .= from
    end
    fmr = copy(from)
    for i=1:nb
        interpolate!(to, from, fmr, interp_t; mpid=mpid, t_split=t_split)
#        @show "autointerp!", i, nb, norm(fmr-to)
        if i != nb
            fmr .= to
        end
    end
end
function interpbufc!(
    t_buf::Vector{Array{OpTuple{N,T}, N}},
    bufdec::Array{OpTuple{N,T}, N},
    interp_t::AbstractVector{I};
    mpid::Union{MPIData, Missing}=missing,
    t_split::Union{Tuple, Missing}=missing
) where {N, T, I <: AbstractInterpolation{T}}

    for buf in t_buf
        interpolate!(buf, copy(buf), bufdec, interp_t, mpid=mpid, t_split=t_split)
    end
end
