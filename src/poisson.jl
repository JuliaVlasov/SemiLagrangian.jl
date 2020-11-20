
include("fftbig.jl")
include("advection.jl")

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
    compute_charge!( self::AdvectionData)

 Compute charge density

 ρ(x,t) = ∫ f(x,v,t) dv

 # Argument
 - `self::AdvectionData` : mutable structure of variables data.
"""
function compute_charge!(self::AdvectionData{T, Nsp, Nv, Nsum}) where {T, Nsp, Nv, Nsum}
    Nsp+Nv == Nsum || thrown(ArgumentError("Nsp=$Nsp Nv=$Nv Nsum=$Nsum we must have Nsp+Nv==Nsum"))
    pvar = getext(self)
    compute_charge!(pvar.rho, self.adv.t_mesh_v, self.data)
    missing
end
"""
    PoissonConst{T, Nsp, Nv}
    PoissonConst(adv::Advection{T, Nsp, Nv, Nsum}; isfftbig=true)

Constant data for the computation of poisson coefficients

# Arguments
- `adv::Advection{T, Nsp, Nv, Nsum}` : Advection constant data
- `isfftbig=true`: if true compute the fttbig structure

# Implementation
- `adv` : Advection constant data
- `v_k' : vector of vector of fourier coefficents for integration for each space dimension
- `fctv_k` : Array of space dimensions of the inverse of the norm of fourier coefficients
- `pfftbig` : Fourier data for space dimensions
"""
struct PoissonConst{T, Nsp, Nv}
    adv
    v_k
    fctv_k
    pfftbig
    function PoissonConst(
    adv::Advection{T, Nsp, Nv, Nsum}; 
    isfftbig=true
) where{T, Nsp, Nv, Nsum}
        Nsp <= Nv || thrown(ArgumentError("Nsp=$Nsp must less or equal to Nv=$Nv"))
        fct_k(v) = im / sum(v .^ 2)
        v_k = vec_k_fft.(adv.t_mesh_sp)
        sz = length.(adv.t_mesh_sp)
        fctv_k = fct_k.(collect(Iterators.product(v_k...)))
        fctv_k[1] = 0
        pfftbig = if isfftbig
            PrepareFftBig(sz, T; numdims=Nsp, dims=ntuple(x->x,Nsp))
        else
            missing
        end
        return new{T,Nsp,Nv}(adv, v_k, fctv_k, pfftbig)
    end
end
"""
    PoissonVar{T, Nsp, Nv}
    PoissonVar(pc::PoissonConst{T, Nsp, Nv})

mutable structure of variable data for the poisson computation

# Arguments
- `pc::PoissonConst{T, Nsp, Nv}` : poisson constant data

# Implementation
- `pc::PoissonConst{T, Nsp, Nv}` : poisson constant data
- `rho::Array{T, Nsp}` : result of the compute_charge that is the sum along velocity dimensions
- `t_elfield::NTuple{Nsp,Array{Complex{T}, Nsp}}` : electric fields initialized at each beginning of velocity advection subseries
"""
mutable struct PoissonVar{T, Nsp, Nv}
    pc::PoissonConst{T, Nsp, Nv}
    rho::Array{T, Nsp}
    t_elfield
    bufcur
    function PoissonVar(pc::PoissonConst{T, Nsp, Nv}) where{T, Nsp, Nv}
        sz = length.(pc.adv.t_mesh_sp)
        rho = Array{T, Nsp}(undef, sz)
        return new{T, Nsp, Nv}(pc, rho, missing, missing)
    end
end
"""

    compute_elfield!( self:AdvectionData)

computation of electric field
    ∇.e = - ρ



# Argument
 - `self::AdvectionData` : mutable structure of variables data.

"""
function compute_elfield!( self::AdvectionData{T, Nsp, Nv, Nsum}) where{T, Nsp, Nv, Nsum}
    pv::PoissonVar{T, Nsp, Nv} = getext(self)
    sz = size(pv.rho)
    pfft = pv.pc.pfftbig
    buf = fftgenall(pfft, pv.rho) .* pv.pc.fctv_k
    buf[1] = 0im
    pv.t_elfield = ntuple(
    x -> real(ifftgenall(pfft, reshape(pv.pc.v_k[x],tupleshape(x,Nsp,sz[x])) .* buf )),
    Nsp
)
    missing
end

"""
    init!(self::AdvectionData{T, Nsp, Nv, Nsum})

Implementation of the interface function that is called at the begining of each advection
    This is implementation for Vlasov-Poisson equation

"""
# TODO premier parametre de type PoissonVar
function init!(self::AdvectionData{T, Nsp, Nv, Nsum}) where{T, Nsp, Nv, Nsum}
    pv::PoissonVar{T, Nsp, Nv} = getext(self)
    state_dim = getstate_dim(self)
    mesh_v = self.adv.t_mesh_v[state_dim]
    if isvelocitystate(self)
        if self.state_dim == 1
            compute_charge!(self)
            compute_elfield!(self)
        end
#        println("v trace init plus")
        pv.bufcur = (getcur_t(self)/step(mesh_v))*pv.t_elfield[state_dim]
    else
        mesh_sp = self.adv.t_mesh_sp[state_dim]
#        println("sp trace init moins")
       pv.bufcur = (-getcur_t(self)/step(mesh_sp))*mesh_v.points
    end
end
function getpoissonvar(adv::Advection)
    pc = PoissonConst(adv)
    return PoissonVar(pc)
end

"""
    getalpha(self::AdvectionData{T, Nsp, Nv, Nsum}, ind) 

Implementation of the interface function that is called before each interpolation in advection

"""
function getalpha(self::AdvectionData{T, Nsp, Nv, Nsum}, ind) where{T, Nsp, Nv, Nsum}
    pv::PoissonVar{T, Nsp, Nv} = getext(self)
    state_dim = getstate_dim(self)
    alpha = if isvelocitystate(self)
        pv.bufcur[ind[1:Nsp]...]
    else
        pv.bufcur[ind[self.state_dim+Nsp]]
    end
#    println("ind=$ind alpha=$alpha")
    return alpha
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
    compute_ee(t_mesh_sp, t_elf)

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

