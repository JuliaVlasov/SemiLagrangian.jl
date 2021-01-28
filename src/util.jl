


# contruct nb iterators thta is the split of 1:lgtot iterator
function splititr(nb, lgtot)
    lg, r = divrem(lgtot,nb)
    return vcat( map(x -> ((x-1)*(lg+1)+1):x*(lg+1), 1:r), map(x -> ((x-1)*lg+r+1):(x*lg+r), (r+1):nb) )
end

function splitvec(nb, v)
    return map(x -> v[x], splititr(nb, length(v)))
end

# construct the permutation of size n that is the transposition of a and b
function transperm(a,b,n)
    p = collect(1:n)
    p[a],p[b] = b, a
    return p
end

# convert a vector to a tuple
totuple(v)=Tuple(x for x in v)

# convert a tuple to a vector
tovector(t)=[x for x in t]

# construct a tuple of size nb, with ones except at index ind the value sz
tupleshape(ind::Int, nb::Int, sz::Int)=Tuple(((x==ind) ? sz : 1) for x in 1:nb)

# construct an array with nb dims, mesh.points on dim=ind, the other dims have a size of one
function tupleshape(ind::Int, nb::Int, v::Vector{T}) where{T}
    return reshape(v, tupleshape(ind, nb, length(v)))
end
function dotprod(t_v::NTuple{N, Vector{T}}) where{N,T}
    #    res = ones(T,totuple(ones(Int,N))) # array of N dimensions with only one one.
    res = ones(T,ntuple(x->1,N)) # array of N dimensions with only one one.
    for (ind, v) in enumerate(t_v)
        res = res .* tupleshape(ind, N, v)
    end
    return res
end
function dotprodother(t_v::NTuple{N, Vector{T}}) where{N,T}
    return prod.(Iterators.product(t_v...))
end



# for bspline
function get_kl_ku(order)
    ku = div(order,2)
    kl = order-1-ku
    return kl, ku
end
