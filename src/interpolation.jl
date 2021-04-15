

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
                bo=view(bufout, :, c)
                bi=view(bufin, :, c)
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
function sol(interp_t::AbstractVector{I}, b::AbstractArray{T,N}) where {T,N,I<:AbstractInterpolation{T}}
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
@inline function getprecal!(v::Vector{T}, bsp::Vector{I}, decf::T) where {T,I<:AbstractInterpolation{T}}
    getprecal!(v, bsp[1], decf)
end
@inline function getprecal!(v::Array{T,N}, bsp::Vector{I}, decf::NTuple{N,T}) where {T, N, I <: AbstractInterpolation{T}}
    v .= dotprod(ntuple(x -> getprecal(bsp[x], decf[x]), N))
end

# modulo for "begin to one" array
modone(ind, n) = (n + ind - 1) % n + 1
gettabmod(lg) = modone.(1:10lg, lg) # 
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

"""
    interpolate!( fp::AbstractVector{T}, 
        fi::AbstractVector{T},
        decint::Int, 
        precal::Vector{T}, 
        interp::AbstractInterpolation{T, CircEdge, order},
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
    origin = -div(order, 2)
    lg = length(fi)
    for i = 1:lg
        indbeg = i + origin + decint + 5lg
        indend = indbeg + order
        fp[i] = sum(res[tabmod[indbeg:indend]] .* precal)
    end
    missing
end
function interpolate!(
    fp::AbstractArray{T,N},
    fi::AbstractArray{T,N},
    decint::NTuple{N,Int},
    precal::Array{T,N},
    interp::Vector{I},
    tabmod = gettabmod.(size(fi)),
) where {T, N, I <: AbstractInterpolation{T,CircEdge}}
    res = sol(interp, fi)
    order = get_order.(interp)
    origin = -div.(order, (2,))
    lg = size(fi)
    lg5 = 5 .* lg
    for i in CartesianIndices(lg)
        indbeg = i.I .+ origin .+ decint .+ lg5
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
    tabmod = gettabmod.(size(fi)),
) where {T, N, I <: AbstractInterpolation{T,CircEdge}}
    return interpolate!(fp, fi, decint[1], precal, interp[1],tabmod[1])
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
    tabmod = gettabmod(length(fi)),
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
    cache_alpha::Union{NTuple{N,T}, T}
    cache_int::Union{NTuple{N,Int}, Int}
    precal::Array{T,N}
    function CachePrecal(interps::Vector{I}, x::T) where{T, I<:AbstractInterpolation{T}}
        N=length(interps)
        cache_alpha = (N == 1) ? zero(T) : ntuple(x->zero(T),N)
        cache_int = (N == 1) ? 0 : ntuple(x->0, N)
        sz = totuple(get_order.(interps) .+ 1)
        precal = zeros(T,sz)
        getprecal!(precal, interps, cache_alpha)
        return new{T,N,I}(interps, cache_alpha, cache_int, precal)
    end
end
@inline function getprecal(self::CachePrecal{T,N}, alpha::NTuple{N,T}) where{T, N}
    if alpha != self.cache_alpha
        self.cache_alpha = alpha
        self.cache_int = Int.(floor.(alpha))
        decfloat = alpha .- self.cache_int
       #        self.cache_precal = get_precal(getinterp(self),decfloat)
        getprecal!(self.precal, self.interps, decfloat)
    end
    return self.cache_int, self.precal
end
@inline function getprecal(self::CachePrecal{T,1}, alpha::T) where{T}
    if alpha != self.cache_alpha
        self.cache_alpha = alpha
        self.cache_int = Int(floor(alpha))
        decfloat = alpha - self.cache_int
       #        self.cache_precal = get_precal(getinterp(self),decfloat)
        getprecal!(self.precal, self.interps[1], decfloat)
    end
    return self.cache_int, self.precal
end


function interpolate!(
    fp::AbstractArray{T,N},
    fi::AbstractArray{T,N},
    dec::Function,
    interp_t::AbstractVector{I},
    tabmod::NTuple{N,Vector{Int}},
    cache::CachePrecal{T,N},
) where {T,N,I<:AbstractInterpolation{T}}

    N == length(interp_t) || thrown(ArgumentError("The number of Interpolation $(length(interp_t)) is different of N=$N"))
    sz = size(fp)
 @show sz   
    res = sol(interp_t, fi)
 
    order = get_order.(interp_t)
    origin = -div.(order,(2,))
    decall = 5 .* sz .+ origin

    for ind in CartesianIndices(sz)
        dint, tab = getprecal(cache, dec(ind))
        deb_i = dint .+ decall .+ ind.I
        end_i = deb_i + order
        fp[ind] = sum(res[ntuple(x -> tabmod[x][deb_i[x]:end_i[x]], N)...] .* tab)
    end
end


function interpolate!(
    fp::AbstractVector{T},
    fi::AbstractVector{T},
    dec::Function,
    interp::AbstractInterpolation{T},
    tabmod::Vector{Int},
    cache::CachePrecal{T,1},
) where {T}
    lg = length(fi)
    
    res::Vector{T} = sol(interp, fi)
 
    order = get_order(interp)
    origin = -div(order,2)
    decall = 5lg + origin

    for ind in CartesianIndices((lg,))
        dint, decfl = getprecal(self, dec(ind)[1])
#        @show decfl, dint
        deb_i =  dint + decall + ind.I[1]
        end_i = deb_i + order
        fp[ind] = sum(res[tabmod[deb_i:end_i]] .* tab)
    end
    fp
end
function interpolate!(
    fp::AbstractVector{T},
    fi::AbstractVector{T},
    dec::Function,
    interp_t::AbstractVector{I},
    tabmod::NTuple{1,Vector{Int}},
    cache::CachePrecal{T,1},
) where {T,I<:AbstractInterpolation{T}}
    interpolate!(fp, fi, dec, interp_t[1], tabmod[1], cache)
end
