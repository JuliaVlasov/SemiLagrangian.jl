
"""

    UniformMesh(start, stop, length)

1D uniform mesh data.

length   : number of points and cells

# Arguments
    - `step`  : size step
    - `points`: Array with node positions
    - `width` : Distance between left and right edges.

"""
struct UniformMesh{T}
    length::Int
    step::T
    points::Vector{T}
    width::T

    function UniformMesh(start::T, stop::T, length::Int) where {T}
        pdeb = range(start, stop = stop, length = length+1)
        points = pdeb[1:end-1]
        step_loc = step(pdeb)
        width = stop - start
        new{T}(length, step_loc, points, width)
    end
end

Base.step(mesh::UniformMesh)=mesh.step
Base.length(mesh::UniformMesh)=mesh.length
points(mesh::UniformMesh)=mesh.points
"""
    vec_k_fft(mesh::UniformMesh{T}) where{T}

Get the fft coefficients of the mesh

# Argument
- `mesh::UniformMesh{T}` : the mesh

# Return
- fft coefficients
"""
function vec_k_fft(mesh::UniformMesh{T}) where{T}
    midx = div(mesh.length,2)
    k = 2T(pi) / (mesh.width)
    res =  k * vcat(0:midx-1, -midx:-1)
    return res
end



