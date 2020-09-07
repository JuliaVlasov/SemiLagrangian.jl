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
    function UniformMesh(start::T, stop::T, length::Int; 
    endpoint = true
) where {T <: AbstractFloat}
        cor = endpoint ? 0 : 1
        pdeb = range(start, stop = stop, length = length + cor)
        points = pdeb[1:end-cor]
        step_loc = T == BigFloat ? step(pdeb) : pdeb.step
        width = stop - start
        new{T}(start, stop, length, step_loc, points, width)
    end
end
export compute_charge!
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
    rho .= dv .* vec(sum(fvx, dims = 1))
    rho .-= sum(rho)/size(rho,1)
    missing
end


export compute_e!

"""

    compute_e!( e, mesh, ρ)

    ∇.e = - ρ

Inplace computation of electric field. Fields e and rho are
already allocated.

"""
function compute_e!(
    e::Vector{T},
    meshx::UniformMesh{T},
    rho::Vector{T},
) where {T <: AbstractFloat}

    nx = meshx.length
    k = 2π / (meshx.stop - meshx.start)
    modes = zeros(T, nx)
    modes .= k .* vcat(0:nx÷2-1, -nx÷2:-1)
    modes[1] = 1.0
    e .= real(ifft(-1im .* fft(rho) ./ modes))

end
