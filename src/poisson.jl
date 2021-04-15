

"""
    compute_charge!( self::AdvectionData)

 Compute charge density

 ρ(x,t) = ∫ f(x,v,t) dv

 # Argument
 - `self::AdvectionData` : mutable structure of variables data.
"""
function compute_charge!(self::AdvectionData{T,N}) where {T,N}
    pvar = getext(self)
    compute_charge!(pvar.rho, self.adv.t_mesh[(div(N,2)+1):N], self.data)
    missing
end

function _get_fctv_k(adv::Advection{T,N,timeopt}) where {T,N,timeopt}
    fct_k(v) = im / sum(v .^ 2)
    Nsp = div(N,2)
    v_k = vec_k_fft.(adv.t_mesh[1:Nsp])
    sz = length.(adv.t_mesh[1:Nsp])
    fctv_k_gen = fct_k.(collect(Iterators.product(v_k...)))
    fctv_k_gen[1] = 0
    return ntuple(x -> reshape(v_k[x], tupleshape(x, Nsp, sz[x])) .* fctv_k_gen, Nsp)
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
struct PoissonConst{T,N, Nsp,Nv}
    adv::Advection
    fctv_k::Any
    v_square::Array{T,Nv}
     pfftbig::Any
    function PoissonConst(
        adv::Advection{T,N,timeopt};
        isfftbig = true,
    ) where {T,N,timeopt}
        N%2 == 0 || thrown(ArgumentError("N=$N must be a multiple of 2"))

        Nsp = Nv = div(N,2)
        fctv_k = _get_fctv_k(adv)
        v_square = dotprod(points.(adv.t_mesh[(Nsp+1):N])) .^ 2 # precompute for ke
        pfftbig = if isfftbig && T != Float64
            PrepareFftBig(sizeall(adv)[1:Nsp], T; numdims = Nsp, dims = ntuple(x -> x, Nsp))
        else
            missing
        end
        return new{T,N, Nsp,Nv}(adv, fctv_k, v_square, pfftbig)
    end
end
getNspNv(pc::PoissonConst{T,N,Nsp,Nv}) where {T,N, Nsp, Nv}= (Nsp, Nv)
"""
    PoissonVar{T, Nsp, Nv} <: AbstractExtDataAdv{T}
    PoissonVar(pc::PoissonConst{T, Nsp, Nv})

mutable structure of variable data for the poisson computation

# Arguments
- `pc::PoissonConst{T, Nsp, Nv}` : poisson constant data

# Implementation
- `pc::PoissonConst{T, Nsp, Nv}` : poisson constant data
- `rho::Array{T, Nsp}` : result of the compute_charge that is the sum along velocity dimensions
- `t_elfield::NTuple{Nsp,Array{Complex{T}, Nsp}}` : electric fields initialized at each beginning of velocity advection subseries
"""
mutable struct PoissonVar{T,N,Nsp,Nv} <: AbstractExtDataAdv{T,N}
    pc::PoissonConst{T,N,Nsp,Nv}
    rho::Array{T,Nsp}
    t_elfield::Union{NTuple{Nsp,Array{T,Nsp}},Missing}
    bufcur
    tupleind

    function PoissonVar(pc::PoissonConst{T,N,Nsp,Nv}) where {T,N,Nsp,Nv}
        sz = length.(pc.adv.t_mesh[1:Nsp])
        rho = Array{T,Nsp}(undef, sz)
        return new{T,N,Nsp,Nv}(pc, rho, missing, missing,missing)
    end
end

getNspNv(_::PoissonVar{T,N,Nsp,Nv}) where {T,N,Nsp,Nv} = (Nsp, Nv)

"""
    compute_elfield!( self:AdvectionData)

computation of electric field
    ∇.e = - ρ

# Argument
 - `self::AdvectionData` : mutable structure of variables data.

"""
function compute_elfield!(self::AdvectionData{T,N}) where {T,N}
    pv = getext(self)
    Nsp, Nv = getNspNv(pv)
    sz = size(pv.rho)
    pfft = pv.pc.pfftbig
    buf = fftgenall(pfft, pv.rho)
    # for i=1:Nsp
    #     size(buf) == size(pv.pc.fctv_k[i]) || thrown(DimensionMismatch("size(buf)=$(size(buf)) size(fctv_k[$i])=$(size(pv.pc.fctv_k[i]))"))
    # end
    pv.t_elfield = ntuple(x -> real(ifftgenall(pfft, pv.pc.fctv_k[x] .* buf)), Nsp)
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
    initcoef!(pv::PoissonVar{T, N,Nsp, Nv}, self::AdvectionData{T, N})

Implementation of the interface function that is called at the begining of each advection
    This is implementation for Vlasov-Poisson equation

"""

function initcoef!(
    pv::PoissonVar{T,N,Nsp,Nv},
    self::AdvectionData{T,N},
) where {T,N,Nsp,Nv}
    st = getst(self)
    adv = self.adv
    if isvelocitystate(self)
        if (Nsp+1) in st.perm[1:st.ndims]
            compute_charge!(self)
            compute_elfield!(self)
        end
        #        println("v trace init plus")
#        pv.bufcur = (getcur_t(self) / step(mesh_v)) * pv.t_elfield[state_dim]
        pv.bufcur = ntuple( x -> (getcur_t(self) / step(adv.t_mesh[st.invp[x]])) * pv.t_elfield[st.invp[x]-Nsp], st.ndims)
#        @show typeof(pv.bufcur)
        #        @show minimum(pv.bufcur),maximum(pv.bufcur) 
    else
#        mesh_sp = self.adv.t_mesh[state_dim]
        #        println("sp trace init moins")
#        pv.bufcur = (-getcur_t(self) / step(mesh_sp)) * mesh_v.points
        pv.tupleind = ntuple( x -> st.perm[st.invp[x]+Nsp]-st.ndims, st.ndims)
        pv.bufcur = ntuple( x -> (-getcur_t(self) / step(adv.t_mesh[st.invp[x]])) * adv.t_mesh[st.invp[x]+Nsp].points, st.ndims)
        #       @show minimum(pv.bufcur),maximum(pv.bufcur) 
#        @show typeof(pv.bufcur)

    end
end
function getpoissonvar(adv::Advection)
    pc = PoissonConst(adv)
    return PoissonVar(pc)
end

function getalpha(
    pv::PoissonVar{T,N,Nsp,Nv},
    self::AdvectionData{T},
    ind,
) where {T,N,Nsp,Nv}
    st = getst(self)
#    @show N,Nsp,Nv, ind.I, pv.tupleind
    if isvelocitystate(self) 
        ntuple(x -> pv.bufcur[x][CartesianIndex(ind.I[end-Nsp+1:end])], size(pv.bufcur,1))
    else
        ntuple(x -> pv.bufcur[x][ind.I[pv.tupleind[x]]], size(pv.bufcur,1))
    end
    #    return isvelocitystate(self) ? pv.bufcur[ind.I[Nsum-Nsp:Nsum-1]...] : pv.bufcur[ind.I[end]]
end
