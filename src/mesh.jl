
"""
$(TYPEDEF)

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
        pdeb = range(start; stop = stop, length = length + 1)
        points = pdeb[1:(end-1)]
        step_loc = step(pdeb)
        width = stop - start
        return new{T}(step_loc, points, width)
    end
end

"""
$(SIGNATURES)

    Base.step(mesh::UniformMesh)

Get the step of the mesh

# Argument
- `mesh::UniformMesh` : the mesh

# Return
- `step` : the step of the mesh that is the difference between two contiguous points
"""
Base.step(mesh::UniformMesh) = mesh.step

"""
$(SIGNATURES)
"""
start(mesh::UniformMesh) = mesh.points[1]

"""
$(SIGNATURES)
"""
stop(mesh::UniformMesh) = mesh.points[end] + mesh.step


"""
$(SIGNATURES)

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
    nx = length(mesh)
    k = 2T(pi) / (width(mesh))
    res = k .* collect(fftfreq(nx, nx))
    return res
end

@inline function newval(valref, val, borne, width)
    diff = val - valref
    return abs(diff) >= borne ? valref + mod(diff + borne, width) - borne : val
end

"""
$(SIGNATURES)
"""
function traitmodend!(
    lg::T2,
    f::Array{T2,N},
) where {N,T2<:Union{OpTuple{N,<:Number},Number}} # T2 must OpTuple or Number
    return f .= mod.(f, lg)
end

"""
$(SIGNATURES)
"""
function traitmodbegin!(
    lg::T2,
    f::Array{T2,N},
) where {N,T2<:Union{OpTuple{N,<:Number},Number}} # T2 must OpTuple or Number
    sz = size(f)
    borne = lg / 2
    flagmod = false
    list = [CartesianIndex(ntuple(x -> 1, N))]
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

"""
$(SIGNATURES)
"""
function traitmodbegin!(mesh::UniformMesh{T}, f::Array{T,N}) where {T,N}
    return traitmodbegin!(width(mesh), f)
end

"""
$(SIGNATURES)
"""
function traitmodend!(mesh::UniformMesh{T}, res::Array{T,N}) where {T,N}
    strt = start(mesh)
    res .-= strt
    traitmodend!(width(mesh), res)
    return res .+= strt
end

"""
$(SIGNATURES)
"""
stdtomesh(mesh::UniformMesh{T}, v) where {T} = mesh.step * v .+ mesh.points[1]

"""
$(SIGNATURES)
"""
meshtostd(mesh::UniformMesh{T}, v) where {T} = (v .- mesh.points[1]) / mesh.step
