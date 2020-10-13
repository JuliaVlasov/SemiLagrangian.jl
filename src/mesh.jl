using Test
#import Statistics: mean
export UniformMesh
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
struct UniformMesh{T, ndims}
    start::T
    stop::T
    length::Int
    step::T
    points::Vector{T}
    width::T
    pfft
    function UniformMesh(start::T, stop::T, length::Int; 
    endpoint = true, 
    isfft = false,
    ndims = 1,
) where {T <: AbstractFloat}
        cor = endpoint ? 0 : 1
        pdeb = range(start, stop = stop, length = length + cor)
        points = pdeb[1:end-cor]
        step_loc = step(pdeb)
        width = stop - start
        pfft = isfft ? PrepareFftBig(length; ndims=ndims, numdim=1) : missing
        new{T, ndims}(start, stop, length, step_loc, points, width, pfft)
    end
end
# export compute_charge!

"""
Energie electrique
|| E(t,.) ||_L2 = 
"""
function compute_ee(mesh_x::UniformMesh, elfield::Vector{T}) where {T}
    dx = mesh_x.step
    return dx * sum(elfield.^2)
end
"""
kinetic Energie 
1/2∫∫ v^2 f(x,v,t) dv dx
"""
function compute_ke( mesh_v::UniformMesh, mesh_x::UniformMesh, fvx::Array{T,2}) where {T}
    dx = mesh_x.step
    dv = mesh_v.step
    sum_x = zeros(T,size(fvx,1))
    sum_x .= sum(fvx, dims = 2)[:,1]
    return dx * dv * sum( mesh_v.points .^ 2 .* sum_x)
end

# function compute_etot( mesh_v::UniformMesh, mesh_x::UniformMesh, fvx::Array{T,2}) where {T}
#     return compute_ke( mesh_v, mesh_x, fvx)+compute_ee( mesh_v, mesh_x, fvx)
# end

"""
    compute_charge!( rho, mesh_v, fvx)

 Compute charge density

 ρ(x,t) = ∫ f(x,v,t) dv

"""
function compute_charge!(
    rho::Array{T,ndims},
    meshv::UniformMesh{T,ndims},
    fvx::Array{T,doublendims},
) where {T <: AbstractFloat, ndims, doublendims}
    nx = size(fvx,ndims*2)
    sh = [(nx), (nx, nx)]
    dims_v = ntuple(x -> x, ndims)
#    println("dims_v=$dims_v typeof(fvx)=$(typeof(fvx))")
    dv = meshv.step
    rho .= dv^ndims * reshape(sum(fvx, dims = dims_v),sh[ndims])
    rho .-= sum(rho)/(size(rho,1)^ndims)
    missing
end

# function selectdim(A, d::NTuple{N,Int}, i::NTuple{N,Int}) where{N}
#     if N == 0
#         A
#     else
#         par = Base.selectdim(A, d[end], i[end])
#         selectdim( par, d[1:end-1], i[1:end-1])
#     end
# end
# @inline tuplejoin(x) = x
# @inline tuplejoin(x, y) = (x..., y...)
# @inline tuplejoin(x, y, z...) = tuplejoin(tuplejoin(x, y), z...)

# export compute_e!


# function Base.PermutedDimsArray(pa::Base.PermutedDimsArray{T,N,perm,iperm,AA}, newperm) where {T,N,perm,iperm,AA<:AbstractArray}
#     totuple(v)=Tuple(x for x in v)
#     tovector(t)=[ x for x in t]
#     return PermutedDimsArray(pa.parent, totuple(tovector(newperm)[tovector(perm)]))
# end

"""

    compute_e!( e, mesh, ρ)

    ∇.e = - ρ

Inplace computation of electric field. Fields e and rho are
already allocated.

"""
function compute_elfield!(
    e::Array{T, ndims},
    meshx::UniformMesh{T, ndims},
    rho::Array{T, ndims},
) where {T <: AbstractFloat, ndims}

    nx = meshx.length
    k = 2T(pi) / (meshx.stop - meshx.start)
    modes = Vector{T}(undef, nx)
    modes .= k .* vcat(0:nx÷2-1, -nx÷2:-1)
    modes[1] = one(T)

    if ndims == 1
        e .= real(ifftgen(meshx.pfft, -1im .* fftgen(meshx.pfft, rho) ./ modes))
    elseif ndims == 2
        buf = real(ifftgen(
    meshx.pfft, 
    -1im .* fftgen(meshx.pfft, Base.PermutedDimsArray(rho,(2,1))) ./ modes
    ))
        e .= real(ifftgen(
    meshx.pfft, 
    -1im .* fftgen(meshx.pfft, Base.PermutedDimsArray(buf,(2,1))) ./ modes
    ))
    else
        println("not yet implemented")
    end

end
