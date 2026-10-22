
"""
    UniformMesh{T}
    UniformMesh(start::T, stop::T, length::Int) where {T}

One-dimensional uniform mesh over the half-open interval `[start, stop)`.
The mesh stores `length` equally spaced points, omitting the right endpoint,
so `width(mesh) == length(mesh) * step(mesh)`.

# Arguments
- `start::T` : left boundary of the mesh
- `stop::T` : right boundary (not itself a mesh point)
- `length::Int` : number of cells of the mesh

# Fields
- `step::T` : spacing between consecutive points
- `points::Vector{T}` : mesh point coordinates
- `width::T` : distance between the left and right boundaries

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
    Base.step(mesh::UniformMesh)

Return the spacing between consecutive mesh points.

# Argument
- `mesh::UniformMesh` : the mesh

# Return
- `step` : difference between consecutive mesh points
"""
Base.step(mesh::UniformMesh) = mesh.step

"""
    start(mesh::UniformMesh)

Return the left boundary of `mesh`, equal to its first stored point.
"""
start(mesh::UniformMesh) = mesh.points[1]

"""
    stop(mesh::UniformMesh)

Return the right boundary of `mesh`, one step beyond its final stored point.
The right boundary is excluded from `mesh.points`.
"""
stop(mesh::UniformMesh) = mesh.points[end] + mesh.step


"""
    Base.length(mesh::UniformMesh)

Return the number of mesh points (equivalently, cells).

# Argument
- `mesh::UniformMesh` : the mesh

# Return
- `length` : number of points and cells
"""
Base.length(mesh::UniformMesh) = length(mesh.points)
"""
    points(mesh::UniformMesh)

Return the coordinates of the mesh points.

# Argument
- `mesh::UniformMesh` : the mesh

# Return
- `points` : vector of `length(mesh)` coordinates, excluding the right boundary
"""
points(mesh::UniformMesh) = mesh.points
"""
    width(mesh::UniformMesh)

Return the distance between the left and right boundaries.

# Argument
- `mesh::UniformMesh` : the mesh

# Return
- `width` : domain width, equal to `length(mesh) * step(mesh)`
"""
width(mesh::UniformMesh) = mesh.width
"""
    vec_k_fft(mesh::UniformMesh{T}) where{T}

Return the Fourier wavenumbers in the ordering used by `FFTW.fftfreq`.

# Argument
- `mesh::UniformMesh{T}` : the mesh

# Return
- `Vector{T}` of wavenumbers, scaled for the mesh's physical width
"""
function vec_k_fft(mesh::UniformMesh{T}) where {T}
    nx = length(mesh)
    k = 2T(pi) / (width(mesh))
    res = k .* collect(fftfreq(nx, nx))
    return res
end

"""
    newval(valref, val, borne, width)

Choose the periodic representative of `val` nearest to `valref`, using a
period of `width` and a half-period bound `borne`. Used internally when
unwrapping periodic values.
"""
@inline function newval(valref, val, borne, width)
    diff = val - valref
    return abs(diff) >= borne ? valref + mod(diff + borne, width) - borne : val
end

"""
    traitmodend!(period, values)

Wrap every value in `values` modulo `period`, in place, and return `values`.
"""
function traitmodend!(
    lg::T2,
    f::Array{T2,N},
) where {N,T2<:Union{OpTuple{N,<:Number},Number}} # T2 must OpTuple or Number
    return f .= mod.(f, lg)
end

"""
    traitmodbegin!(period, values)

Unwrap periodic values in an array by adjusting adjacent values to the
nearest representative within half a period of their predecessor. Mutates
`values` and returns `true` if any value was changed, otherwise `false`.
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
    traitmodbegin!(mesh::UniformMesh, values)

Unwrap periodic values in `values` using `width(mesh)` as the period. Mutates
`values` and returns whether any value changed.
"""
function traitmodbegin!(mesh::UniformMesh{T}, f::Array{T,N}) where {T,N}
    return traitmodbegin!(width(mesh), f)
end

"""
    traitmodend!(mesh::UniformMesh, values)

Wrap `values` into the interval `[start(mesh), stop(mesh))` in place and
return the array.
"""
function traitmodend!(mesh::UniformMesh{T}, res::Array{T,N}) where {T,N}
    strt = start(mesh)
    res .-= strt
    traitmodend!(width(mesh), res)
    return res .+= strt
end

"""
    stdtomesh(mesh::UniformMesh, values)

Convert standard mesh coordinates (grid units, with zero at the first point)
to physical coordinates. Supports scalars, arrays, and compatible tuples.
"""
stdtomesh(mesh::UniformMesh{T}, v) where {T} = mesh.step * v .+ mesh.points[1]

"""
    meshtostd(mesh::UniformMesh, values)

Convert physical coordinates to standard mesh coordinates (grid units, with
zero at the first point). This is the inverse coordinate transformation of
[`stdtomesh`](@ref).
"""
meshtostd(mesh::UniformMesh{T}, v) where {T} = (v .- mesh.points[1]) / mesh.step
