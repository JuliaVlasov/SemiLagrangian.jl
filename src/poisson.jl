
include("fftbig.jl")
include("advection.jl")

using MPI


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

function _get_fctv_k(adv::Advection{T, Nsp, Nv, Nsum, timeopt}) where {T, Nsp, Nv, Nsum, timeopt}
    fct_k(v) = im / sum(v .^ 2)
    v_k = vec_k_fft.(adv.t_mesh_sp)
    sz = length.(adv.t_mesh_sp)
    fctv_k_gen = fct_k.(collect(Iterators.product(v_k...)))
    fctv_k_gen[1] = 0
    return ntuple(x->reshape(v_k[x],tupleshape(x,Nsp,sz[x])) .* fctv_k_gen, Nsp)
end

function _get_perm(adv::Advection{T, Nsp, Nv, Nsum, timeopt}, curstate) where {T, Nsp, Nv, Nsum, timeopt}
    if isvelocity(adv, curstate)
        p = transperm(Nsp+1, curstate, Nsum)
        p = p[vcat(Nsp+1:Nsum, 1:Nsp)]
    else
        p = transperm(curstate+Nsp, Nsum, Nsum)
        p = p[transperm(1, curstate, Nsum)]
    end
    return p
end
function _get_split(adv, curstate)
    if adv.nbsplit != 1
        perm = _get_perm(adv,curstate)
        return splititr(adv.nbsplit, sizeall(adv)[perm][end])
    else
        return missing
    end
end
function _get_t_itrfirst(adv::Advection{T,Nsp, Nv, Nsum, timeopt}, t_itr, curid) where{T,Nsp,Nv,Nsum,timeopt}
    
    szitr = sizeitr(adv)
    if adv.nbsplit == 1
        if isvelocity(adv, curid)
#            println("_get_t_itrfirst trace1")
            return Iterators.product(szitr[1:Nsp]...)
        else
#            println("_get_t_itrfirst trace2 return=$(szitr[curid+Nsp])")
            return szitr[curid+Nsp]
        end
    else
        if isvelocity(adv, curid)
#            println("_get_t_itrfirst trace3")
            return ntuple( x -> Iterators.product((szitr[1:Nsp-1]..., (t_itr[x],)...)...), adv.nbsplit)
        else
#            println("_get_t_itrfirst trace4")
            return t_itr
        end
    end
end


"""
    PoissonConst{T, Nsp, Nv}
    PoissonConst(adv::Advection{T, Nsp, Nv, Nsum}; isfftbig=true)

Constant data for the computation of poisson coefficients

# Arguments
- `adv::Advection{T, Nsp, Nv, Nsum, timeopt}` : Advection constant data
- `isfftbig=true`: if true compute the fttbig structure

# Implementation
- `adv` : Advection constant data
- `v_k' : vector of vector of fourier coefficents for integration for each space dimension
- `fctv_k` : Array of space dimensions of the inverse of the norm of fourier coefficients
- `pfftbig` : Fourier data for space dimensions
"""
struct PoissonConst{T, Nsp, Nv}
    adv
    fctv_k
    pfftbig
    t_perms
    tt_split
    t_itrfirst
    t_itrsecond
    function PoissonConst(
    adv::Advection{T, Nsp, Nv, Nsum, timeopt}; 
    isfftbig=true
) where{T, Nsp, Nv, Nsum, timeopt}
        Nsp == Nv || thrown(ArgumentError("Nsp=$Nsp must be equal to Nv=$Nv"))
        fctv_k = _get_fctv_k(adv)
        pfftbig = if isfftbig && T != Float64
            PrepareFftBig(sizeall(adv)[1:Nsp], T; numdims=Nsp, dims=ntuple(x->x,Nsp))
        else
            missing
        end
        t_perms = ntuple(x -> _get_perm(adv, x), Nsum)
        tt_split = ntuple(x-> _get_split(adv, x), Nsum)
        t_itrfirst = ntuple(x -> _get_t_itrfirst(adv, tt_split[x], x), Nsum)
        t_itrsecond = ntuple(x -> Iterators.product(sizeitr(adv)[t_perms[x]][2:Nsum-1]...), Nsum)
        return new{T,Nsp,Nv}(adv, fctv_k, pfftbig, t_perms, tt_split, t_itrfirst, t_itrsecond)
    end
end
getperm(pc,advd)=pc.t_perms[_getcurrentindice(advd)]
gett_split(pc, advd)=pc.tt_split[_getcurrentindice(advd)]

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
    t_itrfirst
    t_itrsecond
    function PoissonVar(pc::PoissonConst{T, Nsp, Nv}) where{T, Nsp, Nv}
        sz = length.(pc.adv.t_mesh_sp)
        rho = Array{T, Nsp}(undef, sz)
        return new{T, Nsp, Nv}(pc, rho, missing, missing, missing, missing)
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

#    ICI Revenir en arriere (au dernier commit) et y aller plus progressivement

    
    sz = size(pv.rho)
    pfft = pv.pc.pfftbig
    buf = fftgenall(pfft, pv.rho)
    # for i=1:Nsp
    #     size(buf) == size(pv.pc.fctv_k[i]) || thrown(DimensionMismatch("size(buf)=$(size(buf)) size(fctv_k[$i])=$(size(pv.pc.fctv_k[i]))"))
    # end
    pv.t_elfield = ntuple(
    x -> real(ifftgenall(pfft, pv.pc.fctv_k[x] .* buf )),
    Nsp
)
    missing
end

# function init!(pv::PoissonVar{T, Nsp, Nv}, advd::AdvectionData{T, Nsp, Nv, Nsum, timeopt}) where{T,Nsp,Nv,Nsum,timeopt}

#     itrfirst = _getitrfirst(pv, self)
#     if timeopt == MPIOpt || timeopt == SplitThreadsOpt
#         nbsplit = timeopt == MPIOpt ? MPI.Comm_size(MPI.COMM_WORLD) : Threads.nthreads()
#         szlast = getsize(self.adv)[getpermdims(self)][Nsum]
#         if isvelocitystate(advd)
#             t_split = vecsplit(1:szlast)
#         end
#     end
# #         t_itrfirst = map(x->)

# # ICI


# end

"""
    initcoef!(pv::PoissonVar{T, Nsp, Nv}, self::AdvectionData{T, Nsp, Nv, Nsum})

Implementation of the interface function that is called at the begining of each advection
    This is implementation for Vlasov-Poisson equation

"""

function initcoef!(pv::PoissonVar{T, Nsp, Nv}, self::AdvectionData{T, Nsp, Nv, Nsum}) where{T, Nsp, Nv, Nsum}
    state_dim = getstate_dim(self)
    mesh_v = self.adv.t_mesh_v[state_dim]
    if isvelocitystate(self)
        if self.state_dim == 1
            global cl_obs
            clockbegin(cl_obs, 5)
            compute_charge!(self)
            compute_elfield!(self)
            clockend(cl_obs, 5)
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

function getdata(pv::PoissonVar{T, Nsp, Nv}, advd::AdvectionData{T, Nsp, Nv, Nsum}) where{T, Nsp, Nv, Nsum}
    p = getperm(pv.pc,advd)
    if p == 1:Nsum
        # the case of identity permutation no copy is needed
        f = self.data
    else
        ptr = pointer(advd.bufdata)
        f = unsafe_wrap(Array, ptr, sizeall(advd.adv)[p])
        permutedims!(f, advd.data, p)
    end
    return f
end
function copydata!(pv::PoissonVar{T, Nsp, Nv}, advd::AdvectionData{T, Nsp, Nv, Nsum, timeopt}, f) where{T, Nsp, Nv, Nsum, timeopt}
    if timeopt == MPIOpt && advd.adv.nbsplit != 1
        t_split = gett_split(pv.pc, advd)
        comm = MPI.COMM_WORLD
        MPI.Barrier(comm)
        if isbitstype(T)
            # for Float64 or Double64 ... per example
            for i=1:advd.adv.nbsplit
                vbcast = selectdim(f, Nsum, t_split[i])
                MPI.Bcast!(vbcast, i-1, comm)
            end
        else
            # for BigFloat ... per example
            for i=0:advd.adv.nbsplit
                vbcast = selectdim(f, Nsum, t_split[i])
                bufr = MPI.bcast(vbcast, i-1, comm)
                if i != id
                    copy!(vbcast, bufr)
                end
            end
        end
        MPI.Barrier(comm)
    end
    p = getperm(pv.pc, advd)
    pinv = invperm(p)

    if f != advd.data
        permutedims!(advd.data, f, pinv)
    end
end





#Obsolete now
"""
    getalpha(self::AdvectionData{T, Nsp, Nv, Nsum}, ind) 

Implementation of the interface function that is called before each interpolation in advection

"""
function getalpha(self::AdvectionData{T, Nsp, Nv, Nsum}, ind, indthread) where{T, Nsp, Nv, Nsum}
    pv::PoissonVar{T, Nsp, Nv} = getext(self)
    state_dim = getstate_dim(self)
    alpha = if isvelocitystate(self)
        if indthread <= 0
            pv.bufcur[ind[1:Nsp]...]
        else
            pv.bufcur[ind[1:(Nsp-1)]...,((indthread-1)*self.adv.szsplit + ind[Nsp],)...]
        end
    else
        pv.bufcur[ind[self.state_dim+Nsp]]
    end
#    println("ind=$ind alpha=$alpha")
    return alpha
end
getalpha(self::AdvectionData, ind)=getalpha(self,ind,0)

function getprecal(pv::PoissonVar{T, Nsp, Nv}, self::AdvectionData{T, Nsp, Nv, Nsum}, ind) where {T, Nsp, Nv, Nsum}
#    alpha = isvelocitystate(self) ? pv.bufcur[ind...] : pv.bufcur[ind]
    alpha = pv.bufcur[ind...]
    decint = convert(Int, floor(alpha))
    decfloat = alpha - decint
    return decint, get_precal(getinterp(self),decfloat)
end


function getitrfirst(pc::PoissonConst, self::AdvectionData{T,Nsp, Nv, Nsum, timeopt}) where{T,Nsp,Nv,Nsum,timeopt}
    itrfirst = pc.t_itrfirst[_getcurrentindice(self)]
    if pc.adv.nbsplit != 1
        ind = timeopt == MPIOpt ? MPI.Comm_rank(MPI.COMM_WORLD)+1 : Threads.threadid()
        return itrfirst[ind]
    else
#        println("trace good itrfirst=$itrfirst")
        return itrfirst
    end
 
end
getitrfirst(pvar::PoissonVar, advd::AdvectionData)=getitrfirst(pvar.pc, advd)
addcolindend(ind::Tuple,tup)=(:,tup...,ind...)
addcolindend(ind::Int,tup)=(:,tup...,ind)

function getitrsecond(pc::PoissonConst, advd::AdvectionData{T,Nsp, Nv, Nsum}, indfirst) where{T,Nsp,Nv,Nsum}
    perm = getperm(pc,advd)
    szitr = sizeitr(advd.adv)[perm]
    tupmid = isvelocitystate(advd) ? szitr[2:Nv] : szitr[2:Nsum-1]
    return addcolindend.((indfirst,), Iterators.product(tupmid...))
end
getitrsecond(pvar::PoissonVar, advd::AdvectionData, indfirst)=getitrsecond(pvar.pc, advd, indfirst)
    
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

