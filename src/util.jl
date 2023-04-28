import Base: isless, zero, iterate, getindex, +, -, *, /, ==, !=, mod

struct OpTuple{N,T}
    v::NTuple{N,T}
end

Base.show(io::IO, ot::OpTuple) = print(io, ot.v)

(+)(v1::OpTuple, v2::OpTuple) = OpTuple(v1.v .+ v2.v)
(-)(v1::OpTuple, v2::OpTuple) = OpTuple(v1.v .- v2.v)
(-)(v::OpTuple) = OpTuple(0 .- v.v)
(==)(v1::OpTuple, v2::OpTuple) = v1.v == v2.v
(!=)(v1::OpTuple, v2::OpTuple) = v1.v != v2.v
(*)(a::Number, v::OpTuple) = OpTuple(a .* v.v)
(*)(v::OpTuple, a::Number) = OpTuple(v.v .* a)
(/)(v::OpTuple, a::Number) = OpTuple(v.v ./ a)

Base.zero(::Type{OpTuple{N,T}}) where {N,T} = OpTuple(ntuple(x -> zero(T), N))
Base.iterate(v::OpTuple) = iterate(v.v)
Base.iterate(v::OpTuple, state) = iterate(v.v, state)
Base.length(::OpTuple{N}) where {N} = N
Base.getindex(ot::OpTuple, ind...) = getindex(ot.v, ind...)
Base.mod(v1::OpTuple, v2::OpTuple) = OpTuple(mod.(v1.v, v2.v))

# contruct nb iterators thta is the split of 1:lgtot iterator
function splititr(nb, lgtot)
    lg, r = divrem(lgtot, nb)
    return vcat(
        map(x -> ((x-1)*(lg+1)+1):(x*(lg+1)), 1:r),
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

@inline modone(ind, n) = mod(ind - 1, n) + 1
@inline modone(ind::CartesianIndex, n) = CartesianIndex(modone.(ind.I, n))

divone(ind, n) = div(ind - 1, n) + 1
gettabmod(lg) = modone.(1:(5lg), lg)

function getextarray(
    tabor::Array{T2,N},
    decbeg::NTuple{N,Int},
    decend::NTuple{N,Int},
) where {T2,N}
    tab = zeros(T2, size(tabor) .+ decbeg .+ decend)
    inddecbeg = CartesianIndex(decbeg)
    sz = size(tabor)
    for ind in CartesianIndices(tab)
        tab[ind] = tabor[modone(ind - inddecbeg, sz)]
    end
    return tab
end
