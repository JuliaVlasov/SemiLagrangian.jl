
"""
    UniformMesh{T}
    UniformMesh(start::T, stop::T, length::Int) where {T}

1D uniform mesh data.

# Arguments
- `start::T` : beginning of the mesh
- `stop::T` : end of the mesh
- `length::Int` : number of cells of the mesh

# Implementation
- `step::T`  : size step
- `points::Vector{T}`: Array with node positions
- `width::T` : Distance between left and right edges.

"""
struct UniformMesh{T}
    step::T
    points::Vector{T}
    width::T

    function UniformMesh(start::T, stop::T, length::Int) where {T}
        pdeb = range(start, stop = stop, length = length + 1)
        points = pdeb[1:end-1]
        step_loc = step(pdeb)
        width = stop - start
        new{T}(step_loc, points, width)
    end
end
"""
    Base.step(mesh::UniformMesh)

Get the step of the mesh

# Argument
- `mesh::UniformMesh` : the mesh

# Return
- `step` : the step of the mesh that is the difference between two contiguous points
"""
Base.step(mesh::UniformMesh) = mesh.step
"""
    Base.length(mesh::UniformMesh)

Get the length of the mesh

# Argument
- `mesh::UniformMesh` : the mesh

# Return
- `length` : the length of the mesh that is the number of points and cells
"""
Base.length(mesh::UniformMesh) = length(mesh.points)
"""
    points(mesh::UniformMesh)

Get the points of the mesh

# Argument
- `mesh::UniformMesh` : the mesh

# Return
- `points` : the points of the mesh that is the vector of all points of the mesh except the last
"""
points(mesh::UniformMesh) = mesh.points
"""
    width(mesh::UniformMesh)

Get the width of the mesh

# Argument
- `mesh::UniformMesh` : the mesh

# Return
- `width` : the width that is step*length or distance between left and right edges.
"""
width(mesh::UniformMesh) = mesh.width
"""
    vec_k_fft(mesh::UniformMesh{T}) where{T}

Get the fft coefficients of the mesh

# Argument
- `mesh::UniformMesh{T}` : the mesh

# Return
- fft coefficients
"""
function vec_k_fft(mesh::UniformMesh{T}) where {T}
    midx = div(length(mesh), 2)
    k = 2T(pi) / (width(mesh))
    res = k * vcat(0:midx-1, -midx:-1)
    return res
end
