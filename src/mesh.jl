
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

start(mesh::UniformMesh) = mesh.points[1]
stop(mesh::UniformMesh) = mesh.points[end] + mesh.step
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

# @inline function newval(valref, val, borne, width)
#     diff = val - valref
#     if abs(diff) > borne
#         s = sign(diff)
#         divdiff, moddiff = divrem(diff, width)
#         if abs(moddiff) > borne
#             val -= (divdiff + s) * width
#         else
#             val -= divdiff * width
#         end
#     end
#     return val
# end

@inline function newval(valref, val, borne, width)
    diff = val - valref
    return abs(diff) >= 0.5 ? valref + mod(diff +borne, width) - borne : val
end

function traitmodbegin!(mesh::UniformMesh{T}, f::Array{T,N}) where {T,N}
    sz = size(f)
    wdth = width(mesh)
    borne = wdth / 2
    flagmod = false
    list = [ones(Int, N)]
    flagfait = zeros(Bool, sz)
    incr = map(x -> map(y -> y == x, 1:N), 1:N)
    while !isempty(list)
        indice = popfirst!(list)
#        @show indice
        valref = f[indice...]
        for i = 1:N
            newindice = indice + incr[i]
            if newindice[i] <= sz[i]
                if !flagfait[newindice...]
                    val = f[newindice...]
                    nval = newval(valref, val, borne, wdth)
                    if nval != val
                        f[newindice...] = nval
                        flagmod = true
                    end
                    flagfait[newindice...] = true
                    push!(list, newindice)
                end
            end
        end
    end
    return flagmod
end
function traitmodend!(mesh::UniformMesh{T}, res::Array{T,N}) where {T,N}
    strt = start(mesh)
    stp = stop(mesh)
    wdth = width(mesh)
    for ind in CartesianIndices(res)
        val = res[ind]
        if val < strt || val >= stp
            res[ind] = strt + mod(val - strt, wdth)
        end
    end
end

stdtomesh(mesh::UniformMesh{T}, v) where{T}=mesh.step*v .+ mesh.points[1]
meshtostd(mesh::UniformMesh{T}, v) where{T}=(v .- mesh.points[1]) / mesh.step

function traitmodend!(lg::T, res::Array{T,N}) where {T,N}
    res .= mod.(res, lg)
end
function traitmodbegin!(lg::T, f::Array{T,N}) where {T,N}
    sz = size(f)
    borne = lg / 2
    flagmod = false
    list = [CartesianIndex(ntuple(x->1,N))]
    flagfait = zeros(Bool, sz)
    incr = map(x -> CartesianIndex(ntuple(y -> y == x, N)), 1:N)
    while !isempty(list)
        indice = popfirst!(list)
#        @show indice
        valref = f[indice]
        for i = 1:N
            newindice = indice + incr[i]
            if newindice[i] <= sz[i]
                if !flagfait[newindice]
                    val = f[newindice]
                    nval = newval(valref, val, borne, lg)
                    if nval != val
                        f[newindice] = nval
                        flagmod = true
                    end
                    flagfait[newindice] = true
                    push!(list, newindice)
                end
            end
        end
    end
    return flagmod
end
function gettuple_x(t_mesh::NTuple{N,UniformMesh{T}}) where {T,N}
    sz = length.(t_mesh)
    ret = ntuple(x -> zeros(T, sz), N)
    for ind in CartesianIndices(sz)
        for i = 1:N
            ret[i][ind] = t_mesh[i].points[ind.I[i]]
        end
    end
    return ret
end


