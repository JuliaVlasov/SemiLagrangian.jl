
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
    endpoint = false,
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



