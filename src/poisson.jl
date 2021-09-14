
@enum TypePoisson StdPoisson = 1 StdPoisson2d = 2  StdABp = 3

function _get_fctv_k(adv::Advection{T,N,timeopt}) where {T,N,timeopt}
    fct_k(v) = im / sum(v .^ 2)
    Nsp = div(N, 2)
    v_k = vec_k_fft.(adv.t_mesh[1:Nsp])
    sz = length.(adv.t_mesh[1:Nsp])
    fctv_k_gen = fct_k.(collect(Iterators.product(v_k...)))
    fctv_k_gen[1] = 0
    return ntuple(x -> reshape(v_k[x], tupleshape(x, Nsp, sz[x])) .* fctv_k_gen, Nsp)
end

# using SHA
# cod(a::Array) = bytes2hex(sha256(reinterpret(UInt8, collect(Iterators.flatten(a)))))



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
struct PoissonConst{T,N,Nsp,Nv,type,typeadd}
    adv::Advection
    fctv_k::Any
    v_square::Array{T,Nv}
    pfftbig::Any
    abcoef::ABcoef
    function PoissonConst(
        adv::Advection{T,N,timeopt};
        isfftbig = true,
        type::TypePoisson = StdPoisson,
        typeadd = 0,
    ) where {T,N,timeopt}
        N % 2 == 0 || thrown(ArgumentError("N=$N must be a multiple of 2"))

        Nsp = Nv = div(N, 2)
        fctv_k = _get_fctv_k(adv)
        v_square = dotprod(points.(adv.t_mesh[(Nsp+1):N])) .^ 2 # precompute for ke
        pfftbig = if isfftbig && T != Float64
            PrepareFftBig(sizeall(adv)[1:Nsp], T; numdims = Nsp, dims = ntuple(x -> x, Nsp))
        else
            missing
        end

        return new{T,N,Nsp,Nv,type,typeadd}(adv, fctv_k, v_square, pfftbig, ABcoef(typeadd))
    end
end
getNspNv(pc::PoissonConst{T,N,Nsp,Nv}) where {T,N,Nsp,Nv} = (Nsp, Nv)




"""
    PoissonVar{T, N, Nsp, Nv} <: AbstractExtDataAdv{T}
    PoissonVar(pc::PoissonConst{T, N, Nsp, Nv})

mutable structure of variable data for the poisson computation

# Arguments
- `pc::PoissonConst{T, N, Nsp, Nv}` : poisson constant data

# Implementation
- `pc::PoissonConst{T, Nsp, Nv}` : poisson constant data
- `rho::Array{T, Nsp}` : result of the compute_charge that is the sum along velocity dimensions
- `t_elfield::NTuple{Nsp,Array{Complex{T}, Nsp}}` : electric fields initialized at each beginning of velocity advection subseries
"""
mutable struct PoissonVar{T,N,Nsp,Nv,type,typeadd} <: AbstractExtDataAdv
    pc::PoissonConst{T,N,Nsp,Nv,type,typeadd}
    rho::Array{T,Nsp}
    t_elfield::Union{NTuple{Nsp,Array{T,Nsp}},Missing}
    bufcur_sp::Any
    bufcur_v::Any
    tupleind::Any

    function PoissonVar(
        pc::PoissonConst{T,N,Nsp,Nv,type,typeadd},
    ) where {T,N,Nsp,Nv,type,typeadd}
        sz = length.(pc.adv.t_mesh[1:Nsp])
        rho = Array{T,Nsp}(undef, sz)
        
        return new{T,N,Nsp,Nv,type,typeadd}(
            pc,
            rho,
            missing,
            missing,
            missing,
            missing
        )
    end
end
function getpoissonvar(adv::Advection; type::TypePoisson = StdPoisson, typeadd = 0)
    pc = PoissonConst(adv, type = type, typeadd = typeadd)
    return PoissonVar(pc)
end

getNspNv(_::PoissonVar{T,N,Nsp,Nv}) where {T,N,Nsp,Nv} = (Nsp, Nv)
"""
    compute_charge!( self::PoissonVar, advd::AdvectionData)

 Compute charge density

 ρ(x,t) = ∫ f(x,v,t) dv

 # Argument
 - `self::AdvectionData` : mutable structure of variables data.
"""
function compute_charge!(self::PoissonVar{T,N, Nsp,Nv}, advd::AdvectionData) where {T,N, Nsp,Nv}
    compute_charge!(self.rho, self.pc.adv.t_mesh[end-Nv+1:end], advd.data)
    missing
end

"""
    compute_elfield!( self:PoissonVar)

computation of electric field
    ∇.e = - ρ

# Argument
 - `self::PoissonVar` : mutable structure of variables data.

"""
function compute_elfield!(self::PoissonVar{T,N, Nsp,Nv}) where {T,N, Nsp,Nv}
    pfft = self.pc.pfftbig
    buf = fftgenall(pfft, self.rho)
    self.t_elfield = ntuple(x -> real(ifftgenall(pfft, self.pc.fctv_k[x] .* buf)), Nsp)
    missing
end



"""
    initcoef!(pv::PoissonVar{T, N,Nsp, Nv}, self::AdvectionData{T, N})

Implementation of the interface function that is called at the begining of each advection
    This is implementation for Vlasov-Poisson equation

"""

function initcoef!(
    pv::PoissonVar{T,N,Nsp,Nv,StdPoisson},
    advd::AdvectionData{T,N},
) where {T,N,Nsp,Nv}
    st = getst(advd)
    adv = advd.adv
    if isvelocity(advd)
        if (Nsp + 1) in st.perm[1:st.ndims]
            compute_charge!(pv, advd)
            compute_elfield!(pv)
        end
        bc_v = ntuple(
            x ->
                (getcur_t(advd) / step(adv.t_mesh[st.perm[x]])) *
                pv.t_elfield[st.perm[x]-Nsp],
            st.ndims,
        )

        if ismissing(pv.bufcur_v)
            pv.bufcur_v = bc_v
        else
            for x = 1:length(bc_v)
                pv.bufcur_v[x] .= bc_v[x]
            end
        end
    else
        pv.tupleind = ntuple(x -> st.perm[st.invp[x]+Nsp] - st.ndims, st.ndims)
        bc_sp = ntuple(
            x ->
                (-getcur_t(advd) / step(adv.t_mesh[st.invp[x]])) *
                adv.t_mesh[st.invp[x]+Nsp].points,
            st.ndims,
        )
        if ismissing(pv.bufcur_sp)
            pv.bufcur_sp = bc_sp
        else
            for x = 1:length(bc_sp)
                pv.bufcur_sp[x] .= bc_sp[x]
            end
        end
    end
end
@inline function getalpha(
    pv::PoissonVar{T,N,Nsp,Nv,StdPoisson},
    advd::AdvectionData{T},
    ind,
) where {T,N,Nsp,Nv}
    st = getst(advd)
    return if isvelocity(advd)
        ntuple(
            x -> pv.bufcur_v[x][CartesianIndex(ind.I[end-Nsp+1:end])],
            size(pv.bufcur_v, 1),
        )
    else
        ntuple(x -> pv.bufcur_sp[x][ind.I[pv.tupleind[x]]], size(pv.bufcur_sp, 1))
    end
end

function initcoef!(
    pv::PoissonVar{T,N,Nsp,Nv,StdPoisson2d},
    advd::AdvectionData{T,N, timeopt},
) where {T,N,Nsp,Nv, timeopt}
    #    st = getst(self)
    adv = advd.adv
    compute_charge!(pv, advd)
    compute_elfield!(pv)

    bufc_v = (getcur_t(advd) / step(adv.t_mesh[2])) * pv.t_elfield[1]
    bufc_sp = (-getcur_t(advd) / step(adv.t_mesh[1])) * adv.t_mesh[2].points

    if ismissing(advd.bufcur)
        advd.bufcur = zeros(OpTuple{N,T}, sizeall(adv))
    end
    for ind in CartesianIndices(sizeall(adv))
        advd.bufcur[ind] = OpTuple((bufc_sp[ind.I[2]], bufc_v[ind.I[1]]))
    end
end
