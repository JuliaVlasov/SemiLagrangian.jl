
#import Statistics: mean
# export UniformMesh


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
function tupleshape(ind::Int, nb::Int, mesh::UniformMesh{T}) where{T}
    return reshape(mesh.points, tupleshape(ind, nb, mesh.length))
end
# function myproduct(t::NTuple{N,Vector{T}}) where{N,T}
#     res =  Array{Tuple{},N}(undef,totuple(ones(Int,N)))
#     res[1] = ()
#     fct(a,b)=(a...,b)
#     for (ind, tt) in enumerate(t)
#         res = fct.(res,reshape(tt,tupleshape(ind,N,length(tt))))
#     end
#     res
# end
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
"""
Energie electrique
|| E(t,.) ||_L2 = 
"""
function compute_ee(mesh_x::UniformMesh, elfield::Vector{T}) where {T}
    dx = mesh_x.step
    return dx * sum(elfield.^2)
end
function compute_ee(t_mesh_x::NTuple{N,UniformMesh}, elfield::Array{T,N}) where {T,N}
    @assert size(elfield) == length.(t_mesh_x) "size(elfield)=$(size(elfield)) expected:$(length.(t_mesh_x))"
    dx = prod(step, t_mesh_x)
    return dx * sum(elfield.^2)
end
function compute_ee(
    t_mesh_sp::NTuple{N,UniformMesh}, 
    t_elf::NTuple{N,Array{T,N}}
) where {T,N}
    res = zero(T)
    for i=1:N
        res += compute_ee(t_mesh_sp, t_elf[i])
    end
    res
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
function compute_ke_vx( 
    t_mesh_v::NTuple{N, UniformMesh}, 
    t_mesh_x::NTuple{N, UniformMesh}, 
    fvx::Array{T,N2}
) where {T, N, N2}
    @assert N2 == 2N "N2=$N2 must the double of N=$N"
    dx = prod(step, t_mesh_x)
    dv = prod(step, t_mesh_v)
    sum_x = Array{T,N}(undef,size(fvx)[1:N])
    sum_x .= reshape(sum(fvx, dims = totuple((N+1):N2)), size(fvx)[1:N] )
    return dx * dv * sum( dotprod(t_mesh_v) .^ 2 .* sum_x)
end
function compute_ke( 
    t_mesh_sp::NTuple{Nsp, UniformMesh}, 
    t_mesh_v::NTuple{Nv, UniformMesh}, 
    f::Array{T,Nsum}
) where {T, Nsp, Nv, Nsum}
    Nsum == Nsp+Nv || "Nsp=$Nsp, Nv=$Nv, Nsum=$Nsum, we must have Nsum==Nsp+Nv"
    szsp=length.(t_mesh_sp)
    szv=length.(t_mesh_v)
    dsp = prod(step, t_mesh_sp)
    dv = prod(step, t_mesh_v)
    sum_sp = Array{T,Nv}(undef,szv)
    sum_sp .= reshape(sum(f, dims = ntuple(x->x,Nsp)), szv )
    return dsp * dv * sum( dotprod(t_mesh_v) .^ 2 .* sum_sp)
end

# function compute_etot( mesh_v::UniformMesh, mesh_x::UniformMesh, fvx::Array{T,2}) where {T}
#     return compute_ke( mesh_v, mesh_x, fvx)+compute_ee( mesh_v, mesh_x, fvx)
# end

"""
    compute_charge!( rho, mesh_v, fvx)

 Compute charge density

 ρ(x,t) = ∫ f(x,v,t) dv

"""
# function compute_charge!(
#     rho::Vector{T},
#     meshv::UniformMesh,
#     fvx::Array{T,2},
# ) where {T}
#     dv = meshv.step
#     rho .= dv .* vec(sum(fvx, dims = 1))
#     rho .= rho .- mean(rho)
# end
# function compute_charge!(
#     rho::Array{T,N},
#     t_mesh_v::NTuple{N,UniformMesh{T}},
#     f::Array{T,N2},
#     tuple_v
# ) where {T, N, N2}
#     dv = prod(step, t_mesh_v)
#     rho .= dv * reshape(sum(f, dims = tuple_v), size(rho))
#     rho .-= sum(rho)/prod(size(rho))
#     missing
# end
# function compute_charge_vx!(
#     rho::Array{T,N},
#     t_mesh_v::NTuple{N,UniformMesh{T}},
#     fvx::Array{T,N2},
# ) where {T, N, N2}
#     return compute_charge!(rho,t_mesh_v,fvx,ntuple(x -> x, N))
# end
# function compute_charge_xv!(
#     rho::Array{T,N},
#     t_mesh_v::NTuple{N,UniformMesh{T}},
#     fxv,
# ) where {T, N}
#     return compute_charge!(rho,t_mesh_v,fxv,ntuple(x -> N+x, N))
# end
# compute_charge!(rho,t_mesh_v,fvx)=compute_charge_vx!(rho,t_mesh_v,fvx)
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



