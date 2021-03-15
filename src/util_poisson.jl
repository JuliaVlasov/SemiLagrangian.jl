"""
    compute_ke(t_mesh_sp, t_mesh_v, f)

kinetic Energie

∫∫ v^2 f(x,v,t) dv dx

# Arguments
- `t_mesh_sp::NTuple{Nsp, UniformMesh{T}}` : tuple of space meshes
- `t_mesh_v::NTuple{Nv, UniformMesh{T}}` : tuple of velocity meshes
- `f::Array{T,Nsum}` : function data.
"""
function compute_ke( 
    t_mesh_sp::NTuple{Nsp, UniformMesh{T}}, 
    t_mesh_v::NTuple{Nv, UniformMesh{T}}, 
    f::Array{T,Nsum}
) where {T, Nsp, Nv, Nsum}
    Nsum == Nsp+Nv || "Nsp=$Nsp, Nv=$Nv, Nsum=$Nsum, we must have Nsum==Nsp+Nv"
    szv=length.(t_mesh_v)
    dsp = prod(step, t_mesh_sp)
    dv = prod(step, t_mesh_v)
    sum_sp = Array{T,Nv}(undef,szv)
    sum_sp .= reshape(sum(f, dims = ntuple(x->x, Nsp)), szv )
    return  (dsp * dv ) * sum( dotprod(points.(t_mesh_v)) .^ 2 .* sum_sp)
end
"""
    compute_ke(self::AdvectionData)

kinetic Energie

∫∫ v^2 f(x,v,t) dv dx

# Arguments
- `self::AdvectionData` : mutable structure of variables data.
"""
function compute_ke(self::AdvectionData{T, Nsp, Nv, Nsum}) where {T, Nsp, Nv, Nsum}
    Nsum == Nsp+Nv || "Nsp=$Nsp, Nv=$Nv, Nsum=$Nsum, we must have Nsum==Nsp+Nv"
    adv=self.adv
    szv=length.(adv.t_mesh_v)
    dsp = prod(step, adv.t_mesh_sp)
    dv = prod(step, adv.t_mesh_v)
    sum_sp = reshape(sum(getdata(self), dims = ntuple(x->x, Nsp)), szv )
    return (dsp * dv ) * sum(adv.v_square .* sum_sp)
end  

"""
    compute_charge!( rho, mesh_v, fvx)

 Compute charge density

 ρ(x,t) = ∫ f(x,v,t) dv

 # Arguments
 - `rho::Array{T,Nsp}` : result table correctly sized, it is he output
 - `t_mesh_v::NTuple{Nv,UniformMesh{T}}` : velocity meshes
 - `f::Array{T,Nsum}` : input data

"""
function compute_charge!(
    rho::Array{T,Nsp},
    t_mesh_v::NTuple{Nv,UniformMesh{T}},
    f::Array{T,Nsum}
) where {T, Nsp, Nv, Nsum}
    Nsp+Nv == Nsum || thrown(ArgumentError("Nsp=$Nsp Nv=$Nv Nsum=$Nsum we must have Nsp+Nv==Nsum")) 
    dv = prod(step, t_mesh_v)
    rho .= dv * reshape(sum(f, dims = ntuple(x -> Nsp+x, Nv)), size(rho))
    rho .-= sum(rho)/prod(size(rho))
    missing
end

"""
    compute_elfield!(elf::Array{T,1}, mesh::UniformMesh{T}, rho::Array{T,1}) where{T}

computation of electric field
    ∇.e = - ρ

# Argument
 - `self::AdvectionData` : mutable structure of variables data.

"""
# function compute_elfield!( self::AdvectionData{T, Nsp, Nv, Nsum}) where{T, Nsp, Nv, Nsum}
#     pv::PoissonVar{T, Nsp, Nv} = getext(self)

#     sz = size(pv.rho)
#     pfft = pv.pc.pfftbig
#     buf = fftgenall(pfft, pv.rho)
#     # for i=1:Nsp
#     #     size(buf) == size(pv.pc.fctv_k[i]) || thrown(DimensionMismatch("size(buf)=$(size(buf)) size(fctv_k[$i])=$(size(pv.pc.fctv_k[i]))"))
#     # end
#     pv.t_elfield = ntuple(
#     x -> real(ifftgenall(pfft, pv.pc.fctv_k[x] .* buf )),
#     Nsp
# )
#     missing
# end
function compute_elfield(
    t_mesh_x::NTuple{N,UniformMesh{T}},
    rho::Array{T, N},
    pfft
) where {T <: AbstractFloat, N}

    fct_k(v)= im/sum(v.^2)

    v_k = vec_k_fft.(t_mesh_x)
    sz = length.(t_mesh_x)

    buf = fftgenall(pfft, rho) .* fct_k.(collect(Iterators.product(v_k...)))
    buf[1] = 0im

    return ntuple( 
    x -> real(ifftgenall(pfft, reshape(v_k[x],tupleshape(x,N,sz[x])) .* buf )),
    N
)
end
function compute_elfield!(elf::Array{T,1}, mesh::UniformMesh{T}, rho::Array{T,1}) where{T}
    elf .= compute_elfield((mesh,), rho, PrepareFftBig(size(rho), zero(T); numdims=1))[1]
end


"""
    compute_ee(t_mesh_sp, t_elf)

compute electric enegie
|| E(t,.) ||_L2

# Arguments
- `t_mesh_sp::NTuple{N,UniformMesh{T}}` : tuple of space meshes
- `t_elf::NTuple{N,Array{T,N}}` : tuple of electric field
"""
function compute_ee(
    t_mesh_sp::NTuple{N,UniformMesh{T}}, 
    t_elf::NTuple{N,Array{T,N}}
) where {T,N}
    res = zero(T)
    dx = prod(step, t_mesh_sp)
    for i=1:N
        res += sum(t_elf[i].^2)
    end
    dx * res
end

"""
    compute_ee(self::AdvectionData)

compute electric enegie
|| E(t,.) ||_L2

# Argument
- `self::AdvectionData` : veriable advection data structure.

"""
function compute_ee(self::AdvectionData)
    adv = self.adv
    pvar = getext(self)
    dx = prod(step, adv.t_mesh_sp)
    return dx * sum(map(x->sum(x.^2), pvar.t_elfield))
end

