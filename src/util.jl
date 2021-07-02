
#for debug only
# using CRC32c
# cksum(t::Array)=crc32c(collect(reinterpret(UInt8,t)))
# cksum(t::Array, x::Int32)=crc32c(collect(reinterpret(UInt8,t)), x)

import Base.getindex

# contruct nb iterators thta is the split of 1:lgtot iterator
function splititr(nb, lgtot)
    lg, r = divrem(lgtot, nb)
    return vcat(
        map(x -> ((x-1)*(lg+1)+1):x*(lg+1), 1:r),
        map(x -> ((x-1)*lg+r+1):(x*lg+r), (r+1):nb),
    )
end

function splitvec(nb, v)
    return map(x -> v[x], splititr(nb, length(v)))
end

# construct the permutation of size n that is the transposition of a and b
function transposition(a, b, n)
    p = collect(1:n)
    p[a], p[b] = b, a
    return p
end

# convert a vector to a tuple
totuple(v) = Tuple(x for x in v)

# convert a tuple to a vector
tovector(t) = [x for x in t]

# construct a tuple of size nb, with ones except at index ind the value sz
tupleshape(ind::Int, nb::Int, sz::Int) = Tuple(((x == ind) ? sz : 1) for x = 1:nb)

# construct an array with nb dims, mesh.points on dim=ind, the other dims have a size of one
function tupleshape(ind::Int, nb::Int, v::Vector{T}) where {T}
    return reshape(v, tupleshape(ind, nb, length(v)))
end
function dotprod(t_v::NTuple{N,Vector{T}}) where {N,T}
    #    res = ones(T,totuple(ones(Int,N))) # array of N dimensions with only one one.
    res = ones(T, ntuple(x -> 1, N)) # array of N dimensions with only one one.
    for (ind, v) in enumerate(t_v)
        res = res .* tupleshape(ind, N, v)
    end
    return res
end
dotprod(t_v::Tuple{Vector{T}}) where {T} = t_v[1]
function dotprodother(t_v::NTuple{N,Vector{T}}) where {N,T}
    return prod.(Iterators.product(t_v...))
end

# modulo or div for "begin to one" array
# an other argument for "begin to zero" array ...
@inline modone(ind, n) = mod(ind - 1, n) + 1
@inline modone(ind::CartesianIndex, n) = CartesianIndex(modone.(ind.I, n))
divone(ind, n) = div(ind - 1, n) + 1
gettabmod(lg) = modone.(1:3lg, lg)

# dotprod(v_v::Vector{Vector{T}}) where{T}=dotprod(totuple(v_v))
# modone(x,n)=(x-1)%n+1
# struct CircularArray{T,N} <: AbstractArray{T,N} #inherits from AbstractArray
#     x::AbstractArray{T,N}
#     tabmod::NTuple{N,Vector{Int}}
#     function CircularArray(x::AbstractArray{T,N}) where {T,N} 
# #creates the type only with a vector x
#         sz = size(x)
#         tabmod = ntuple( i -> modone.(1:10*sz[i], sz[i]), N)
#         return new{T,N}(x, tabmod)
#     end
# end
# modtabmod(tm::AbstractVector{Int}, i::Int)=tm[i]

# Base.size(A::CircularArray) = 10 .* size(A.x) #important: use of Base.function

# Base.length(A::CircularArray)=length(A.x)

# function Base.getindex(A::CircularArray, I::Vararg{Int, N}) where N # implements A[I]
#     return Base.getindex(A.x,modtabmod.(A.tabmod, I)...) #this is the magic operation
# end
# function Base.getindex(A::CircularArray, I::Tuple{N, Int}) where N # implements A[I]
#      return Base.getindex(A.x,modtabmod.(A.tabmod, I)...) #this is the magic operation
# end

# Base.getindex(A::CircularArray{T,nd}, CI::CartesianIndex{nd}) where{T,nd}=getindex(A, CI.I)

# Base.getindex(A::CircularArray, I) = (A[i] for i in I) #A[1:5], for example

# function Base.setindex!(A::CircularArray,value,I::Vararg{Int, N}) where N # A[I] = value
#     I2 = size(A)
#     return Base.setindex!(A.x,value,modtabmod.(A.tabmod,I)...)
# end
# function Base.setindex!(A::CircularArray,value,I::Tuple{N, Int}) where N # A[I] = value
#     I2 = size(A)
#     return Base.setindex!(A.x,value,modtabmod.(A.tabmod,I)...)
# end
# Base.setindex!(A::CircularArray,value, CI::CartesianIndex)= setindex!(A,value,CI.I)

# Base.IndexStyle(::Type{CircularArray}) = IndexCartesian()
# a finir

function getextarray(
    tabor::Array{T,N},
    decbeg::NTuple{N,Int},
    decend::NTuple{N,Int},
) where {T,N}
    tab = zeros(T, size(tabor) .+ decbeg .+ decend)
    inddecbeg = CartesianIndex(decbeg)
    sz = size(tabor)
    for ind in CartesianIndices(tab)
        tab[ind] = tabor[modone(ind - inddecbeg, sz)]
    end
    return tab
end

# struct ExtArray{T,N} <: AbstractArray{T,N}
#     tab::Array{T,N}
#     inddecbegin::CartesianIndex{N}
#     function ExtArray(
#         tabor::Array{T,N},
#         decbeg::NTuple{N,Int},
#         decend::NTuple{N,Int},
#     ) where {T,N}
#         tab = zeros(T, size(tabor) .+ decbeg .+ decend)
#         inddecbeg = CartesianIndex(decbeg)
#         sz = size(tabor)
#         for ind in CartesianIndices(tab)
#             tab[ind] = tabor[modone(ind + inddecbeg, sz)]
#         end
#         return new{T,N}(tab, inddecbeg)
#     end
# end
# Base.getindex(A::ExtArray{T,N}, ind::CartesianIndex{N}) where {T,N} =
#     A.tab[ind+A.inddecbegin]
# Base.axes(A::ExtArray{T,N}) where {T,N} =
#     ntuple(x -> (-A.inddecbegin[x]+1):(size(A.tab, x)-A.inddecbegin[x]), N)
# Base.size(A::ExtArray{T,N}) where {T,N} = size(A.tab)
# getarray(A::ExtArray) = A.tab
