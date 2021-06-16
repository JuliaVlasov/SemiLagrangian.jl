
@enum TypePoisson StdPoisson=1 StdPoisson2d=2 StdOrder2_1=3 StdOrder2_2=4

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
struct PoissonConst{T,N, Nsp,Nv, type}
    adv::Advection
    fctv_k::Any
    v_square::Array{T,Nv}
    pfftbig::Any
    function PoissonConst(
        adv::Advection{T,N,timeopt};
        isfftbig = true,
        type::TypePoisson=StdPoisson
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
        return new{T,N, Nsp, Nv, type}(adv, fctv_k, v_square, pfftbig)
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
mutable struct PoissonVar{T,N,Nsp,Nv, type} <: AbstractExtDataAdv
    pc::PoissonConst{T,N,Nsp,Nv, type}
    rho::Array{T,Nsp}
    t_elfield::Union{NTuple{Nsp,Array{T,Nsp}},Missing}
    t_elfprec::Union{NTuple{Nsp,Array{T,Nsp}},Missing}
    bufcur_sp
    bufcur_v
    tupleind

    function PoissonVar(pc::PoissonConst{T,N,Nsp,Nv, type}) where {T,N,Nsp,Nv, type}
        sz = length.(pc.adv.t_mesh[1:Nsp])
        rho = Array{T,Nsp}(undef, sz)
        return new{T,N,Nsp,Nv, type}(pc, rho, missing, missing, missing,missing,missing)
    end
end
function getpoissonvar(adv::Advection; type::TypePoisson=StdPoisson)
    pc = PoissonConst(adv, type=type)
    return PoissonVar(pc)
end

getNspNv(_::PoissonVar{T,N,Nsp,Nv}) where {T,N,Nsp,Nv} = (Nsp, Nv)
"""
    compute_charge!( self::AdvectionData)

 Compute charge density

 ρ(x,t) = ∫ f(x,v,t) dv

 # Argument
 - `self::AdvectionData` : mutable structure of variables data.
"""
function compute_charge!(self::AdvectionData{T,N}) where {T,N}
    pvar = getext(self)
    Nsp,Nv = getNspNv(pvar)
    compute_charge!(pvar.rho, self.adv.t_mesh[end-Nv+1:end], self.data)
    missing
end

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
    # if ismissing(pv.t_elfield)
    #     pv.t_elfield = elf
    # else
    #     for x in 1:Nsp
    #         pv.t_elfield[x] .= elf[x]
    #     end
    # end
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
    pv::PoissonVar{T,N,Nsp,Nv, StdPoisson},
    self::AdvectionData{T,N},
) where {T,N,Nsp,Nv}
    st = getst(self)
    adv = self.adv
    if isvelocity(self)
        if (Nsp+1) in st.perm[1:st.ndims]
#            @show cksum(self.data)
            compute_charge!(self)
#            @show cksum(self.parext.rho), size(self.parext.rho)
            compute_elfield!(self)
#            @show cksum.(self.parext.t_elfield)
        end
        #        println("v trace init plus")
#        pv.bufcur = (getcur_t(self) / step(mesh_v)) * pv.t_elfield[state_dim]
#        @show st.ind, st.perm, Nsp
        bc_v = ntuple( x -> (getcur_t(self) / step(adv.t_mesh[st.perm[x]])) * pv.t_elfield[st.perm[x]-Nsp], st.ndims)

        if ismissing(pv.bufcur_v)
            pv.bufcur_v = bc_v
        else
            for x in 1:length(bc_v)
                pv.bufcur_v[x] .= bc_v[x]
            end
        end
#        @show self.state_gen, cksum(pv.bufcur[1]), size(pv.bufcur[1])
#        @show typeof(pv.bufcur)
#        @show minimum(pv.bufcur[1]),maximum(pv.bufcur[1]) 
    else
#        mesh_sp = self.adv.t_mesh[state_dim]
        #        println("sp trace init moins")
#        pv.bufcur = (-getcur_t(self) / step(mesh_sp)) * mesh_v.points
        pv.tupleind = ntuple( x -> st.perm[st.invp[x]+Nsp]-st.ndims, st.ndims)
        bc_sp = ntuple( x -> (-getcur_t(self) / step(adv.t_mesh[st.invp[x]])) * adv.t_mesh[st.invp[x]+Nsp].points, st.ndims)
        if ismissing(pv.bufcur_sp)
            pv.bufcur_sp = bc_sp
        else
            for x in 1:length(bc_sp)
                pv.bufcur_sp[x] .= bc_sp[x]
            end
        end


#        @show self.state_gen, cksum(pv.bufcur[1]), size(pv.bufcur[1])
       #       @show minimum(pv.bufcur),maximum(pv.bufcur) 
#        @show typeof(pv.bufcur)

    end
end
@inline function getalpha(
    pv::PoissonVar{T,N,Nsp,Nv, StdPoisson},
    self::AdvectionData{T},
    ind,
) where {T,N,Nsp,Nv}
    st = getst(self)
#    @show N,Nsp,Nv, ind.I, pv.tupleind, isvelocity(self)
    if isvelocity(self) 
        ntuple(x -> pv.bufcur_v[x][CartesianIndex(ind.I[end-Nsp+1:end])], size(pv.bufcur_v,1))
    else
        ntuple(x -> pv.bufcur_sp[x][ind.I[pv.tupleind[x]]], size(pv.bufcur_sp,1))
    end
    #    return isvelocitystate(self) ? pv.bufcur[ind.I[Nsum-Nsp:Nsum-1]...] : pv.bufcur[ind.I[end]]
end

function initcoef!(
    pv::PoissonVar{T,N,Nsp,Nv, StdPoisson2d},
    self::AdvectionData{T,N},
) where {T,N,Nsp,Nv}
#    st = getst(self)
    adv = self.adv
    compute_charge!(self)
    compute_elfield!(self)

    pv.bufcur_v = (getcur_t(self) / step(adv.t_mesh[2])) * pv.t_elfield[1]
    
    pv.bufcur_sp = (-getcur_t(self) / step(adv.t_mesh[1])) * adv.t_mesh[2].points

end
function initcoef!(
    pv::PoissonVar{T,N,Nsp,Nv, StdOrder2_1},
    self::AdvectionData{T,N},
) where {T,N,Nsp,Nv}
#    st = getst(self)
    adv = self.adv
    compute_charge!(self)
    compute_elfield!(self)
    if ismissing(pv.t_elfprec)
        elf = pv.t_elfield[1]
    else
        elf = (3pv.t_elfield[1] - pv.t_elfprec[1])/2
    end

    pv.t_elfprec = pv.t_elfield
    sz = sizeall(self.adv)
    if ismissing(pv.bufcur_v)
        pv.bufcur_v = zeros(T,sz)
        pv.bufcur_sp = zeros(T,sz)
    end

    bufc_v = (getcur_t(self) / step(adv.t_mesh[2])) * elf
    bufc_sp = (-getcur_t(self) / 2step(adv.t_mesh[1])) * adv.t_mesh[2].points

    interp = adv.t_interp[1]

    lgn = zeros(T,sz[1])

    for i=1:sz[2]
        interpolate!(lgn, bufc_v, bufc_sp[i], interp)
        pv.bufcur_v[:,i] .= lgn
    end

    for i=1:sz[1], j=1:sz[2]
        pv.bufcur_sp[i,j] = (-getcur_t(self)/2step(adv.t_mesh[1])) *( 2adv.t_mesh[2].points[j] + pv.bufcur_v[i,j])
    end
    missing
end
function initcoef!(
    pv::PoissonVar{T,N,Nsp,Nv, StdOrder2_2},
    self::AdvectionData{T,N},
) where {T,N,Nsp,Nv}
#    st = getst(self)
    adv = self.adv
    compute_charge!(self)
    compute_elfield!(self)

    sz = sizeall(self.adv)
 
    bufc_v = (getcur_t(self) / 2step(adv.t_mesh[2])) * pv.t_elfield[1]
    bufc_sp = (-getcur_t(self) / 2step(adv.t_mesh[1])) * adv.t_mesh[2].points

    result = zeros(T,sz)
    interpolate!(result, self.data, ind -> ( bufc_sp[ind.I[2]], bufc_v[ind.I[1]]), adv.t_interp)

    rho = zeros(T,sz[1])

    elf =zeros(T,sz[1])

    compute_charge!(rho, (adv.t_mesh[2],), result)
    compute_elfield!(elf, adv.t_mesh[1], rho)
    pv.bufcur_v = (getcur_t(self) / step(adv.t_mesh[2]))*elf
    pv.bufcur_sp = [(-getcur_t(self) / 2step(adv.t_mesh[1])) * (2*adv.t_mesh[2].points[j] + pv.bufcur_v[i]) for i=1:sz[1],j=1:sz[2] ]
    missing
end

@inline function getalpha(
    pv::PoissonVar{T,N,Nsp,Nv, StdPoisson2d},
    self::AdvectionData{T},
    indext,
    ind
) where {T,N,Nsp,Nv}
    @assert Nsp == Nv == 1 "Nsp=$Nsp Nv=$Nv they must be equal to one"
#    @show ind.I
    return ( pv.bufcur_sp[ind.I[2]], pv.bufcur_v[ind.I[1]])
end
@inline function getalpha(
    pv::PoissonVar{T,N,Nsp,Nv, StdOrder2_2},
    self::AdvectionData{T},
    indext,
    ind
) where {T,N,Nsp,Nv}
    @assert Nsp == Nv == 1 "Nsp=$Nsp Nv=$Nv they must be equal to one"
#    @show ind.I
    return ( pv.bufcur_sp[ind], pv.bufcur_v[ind.I[1]])
end
@inline function getalpha(
    pv::PoissonVar{T,N,Nsp,Nv, StdOrder2_1},
    self::AdvectionData{T},
    indext,
    ind
) where {T,N,Nsp,Nv}
    @assert Nsp == Nv == 1 "Nsp=$Nsp Nv=$Nv they must be equal to one"
#    @show ind.I
    return ( pv.bufcur_sp[ind], pv.bufcur_v[ind])
end
