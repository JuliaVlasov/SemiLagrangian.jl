
#import Statistics: mean
# export UniformMesh

include("clockobs.jl")

cl_obs = ClockObs(8)


"""

    UniformMesh(start, stop, length)

1D uniform mesh data.

length   : number of points
length-1 : number of cells

If you want remove the last point for periodic domain, set endpoint=false

    - `step`  : size step
    - `points`: Array with node positions
    - `width` : Distance between left and right edges.

"""
struct UniformMesh{T}
    start::T
    stop::T
    length::Int
    step::T
    points::Vector{T}
    width::T
#    pfft
    function UniformMesh(start::T, stop::T, length::Int; 
    endpoint = true,
) where {T}
        cor = endpoint ? 0 : 1
        pdeb = range(start, stop = stop, length = length + cor)
        points = pdeb[1:end-cor]
        step_loc = step(pdeb)
        width = stop - start
        # pfft = if isfft
        #     PrepareFftBig(length; numdims=ndims, dims=ntuple(identity,ndims))
        # else
        #     missing
        # end
        new{T}(start, stop, length, step_loc, points, width)
    end
end

Base.step(mesh::UniformMesh)=mesh.step
Base.length(mesh::UniformMesh)=mesh.length
points(mesh::UniformMesh)=mesh.points
# println("trace121")
function vec_k_fft(mesh::UniformMesh{T}) where{T}
    midx = div(mesh.length,2)
    k = 2T(pi) / (mesh.stop - mesh.start)
    res =  k * vcat(0:midx-1, -midx:-1)
    return res
end

# export compute_charge!

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

# TODO supprimer les fonctions suivantes et les remplacer par les precedentes 
# motivation : ORTHOGONALITE !!!!
function tupleshape(ind::Int, nb::Int, mesh::UniformMesh{T}) where{T}
    return reshape(mesh.points, tupleshape(ind, nb, mesh.length))
end
function dotprod(t_mesh::NTuple{N, UniformMesh{T}}) where{N,T}
#    res = ones(T,totuple(ones(Int,N))) # array of N dimensions with only one one.
    res = ones(T,ntuple(x->1,N)) # array of N dimensions with only one one.
    for (ind, mesh) in enumerate(t_mesh)
        res = res .* tupleshape(ind, N, mesh)
    end
    return res
end
function dotprodother(t_mesh::NTuple{N, UniformMesh{T}}) where{N,T}
    return prod.(Iterators.product(points.(t_mesh)...))
end
    

