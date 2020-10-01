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
struct UniformMesh{T}
    start::T
    stop::T
    length::Int
    step::T
    points::Vector{T}
    width::T
    pfft
    function UniformMesh(start::T, stop::T, length::Int; 
    endpoint = true, 
    isfft = false
) where {T <: AbstractFloat}
        cor = endpoint ? 0 : 1
        pdeb = range(start, stop = stop, length = length + cor)
        points = pdeb[1:end-cor]
        step_loc = T == BigFloat ? step(pdeb) : pdeb.step
        width = stop - start
        pfft = isfft ? PrepareFftBig(length; ndims=1, numdim=1) : missing
        new{T}(start, stop, length, step_loc, points, width, pfft)
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
    rho::Vector{T},
    meshv::UniformMesh,
    fvx::Array{T,2},
) where {T <: AbstractFloat}
    dv = meshv.step
    rho .= dv * vec(sum(fvx, dims = 1))
    rho .-= sum(rho)/size(rho,1)
    missing
end


# export compute_e!

"""

    compute_e!( e, mesh, ρ)

    ∇.e = - ρ

Inplace computation of electric field. Fields e and rho are
already allocated.

"""
function compute_elfield!(
    e::Vector{T},
    meshx::UniformMesh{T},
    rho::Vector{T},
) where {T <: AbstractFloat}

    nx = meshx.length
    k = 2T(pi) / (meshx.stop - meshx.start)
    modes = zeros(T, nx)
    modes .= k .* vcat(0:nx÷2-1, -nx÷2:-1)
    modes[1] = one(T)
    e .= real(ifftgen(meshx.pfft, -1im .* fftgen(meshx.pfft, rho) ./ modes))

end
