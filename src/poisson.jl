
@enum TypePoisson StdPoisson = 1 StdPoisson2d = 2 StdOrder2_1 = 3 StdOrder2_2 = 4 StdAB = 5 StdPoisson2dTry =
    6 StdABNew = 7 StdAB2=8

function _get_fctv_k(adv::Advection{T,N,timeopt}) where {T,N,timeopt}
    fct_k(v) = im / sum(v .^ 2)
    Nsp = div(N, 2)
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
    t_elfprec::Union{NTuple{Nsp,Array{T,Nsp}},Missing}
    t_data::Vector{Array{T,N}}
    t_sp::Union{Vector{Array{T,N}},Array{T,N}}
    t_v::Union{Vector{Array{T,N}},Array{T,N}}
    t_bufc_sp::Vector{Array{T,N}}
    t_bufc_v::Vector{Array{T,N}}
    bufcur_sp::Any
    bufcur_v::Any
    tupleind::Any

    function PoissonVar(
        pc::PoissonConst{T,N,Nsp,Nv,type,typeadd},
    ) where {T,N,Nsp,Nv,type,typeadd}
        sz = length.(pc.adv.t_mesh[1:Nsp])
        rho = Array{T,Nsp}(undef, sz)
        sz = sizeall(pc.adv)
        t_data = []
        t_bufc_sp = []
        t_bufc_v = []
        if type == StdPoisson2dTry
            t_sp = zeros(T, sz)
            t_v = zeros(T, sz)
            for ind in CartesianIndices(sz)
                t_sp[ind] = ind.I[1]-1
                t_v[ind] = ind.I[2]-1
            end
        elseif type == StdAB
            t_sp = [zeros(T,sz) for i=0:typeadd]
            t_v = [zeros(T,sz) for i=0:typeadd]
            for ind in CartesianIndices(sz), i=1:typeadd+1
                t_sp[i][ind] = ind.I[1]-1
                t_v[i][ind] = ind.I[2]-1
            end
        else
            t_sp = zeros(T,ntuple(x-> 1,N))
            t_v = zeros(T,ntuple(x-> 1,N))
        end
        return new{T,N,Nsp,Nv,type,typeadd}(
            pc,
            rho,
            missing,
            missing,
            t_data,
            t_sp,
            t_v,
            t_bufc_sp,
            t_bufc_v,
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
    compute_charge!( self::AdvectionData)

 Compute charge density

 ρ(x,t) = ∫ f(x,v,t) dv

 # Argument
 - `self::AdvectionData` : mutable structure of variables data.
"""
function compute_charge!(self::AdvectionData{T,N}) where {T,N}
    pvar = getext(self)
    Nsp, Nv = getNspNv(pvar)
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
    pv::PoissonVar{T,N,Nsp,Nv,StdPoisson},
    self::AdvectionData{T,N},
) where {T,N,Nsp,Nv}
    st = getst(self)
    adv = self.adv
    if isvelocity(self)
        if (Nsp + 1) in st.perm[1:st.ndims]
            #            @show cksum(self.data)
            compute_charge!(self)
            #            @show cksum(self.parext.rho), size(self.parext.rho)
            compute_elfield!(self)
            #            @show cksum.(self.parext.t_elfield)
        end
        #        println("v trace init plus")
        #        pv.bufcur = (getcur_t(self) / step(mesh_v)) * pv.t_elfield[state_dim]
        #        @show st.ind, st.perm, Nsp
        bc_v = ntuple(
            x ->
                (getcur_t(self) / step(adv.t_mesh[st.perm[x]])) *
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
        #        @show self.state_gen, cksum(pv.bufcur[1]), size(pv.bufcur[1])
        #        @show typeof(pv.bufcur)
        #        @show minimum(pv.bufcur[1]),maximum(pv.bufcur[1]) 
    else
        #        mesh_sp = self.adv.t_mesh[state_dim]
        #        println("sp trace init moins")
        #        pv.bufcur = (-getcur_t(self) / step(mesh_sp)) * mesh_v.points
        pv.tupleind = ntuple(x -> st.perm[st.invp[x]+Nsp] - st.ndims, st.ndims)
        bc_sp = ntuple(
            x ->
                (-getcur_t(self) / step(adv.t_mesh[st.invp[x]])) *
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


        #        @show self.state_gen, cksum(pv.bufcur[1]), size(pv.bufcur[1])
        #       @show minimum(pv.bufcur),maximum(pv.bufcur) 
        #        @show typeof(pv.bufcur)

    end
end
@inline function getalpha(
    pv::PoissonVar{T,N,Nsp,Nv,StdPoisson},
    self::AdvectionData{T},
    ind,
) where {T,N,Nsp,Nv}
    st = getst(self)
    #    @show N,Nsp,Nv, ind.I, pv.tupleind, isvelocity(self)
    if isvelocity(self)
        ntuple(
            x -> pv.bufcur_v[x][CartesianIndex(ind.I[end-Nsp+1:end])],
            size(pv.bufcur_v, 1),
        )
    else
        ntuple(x -> pv.bufcur_sp[x][ind.I[pv.tupleind[x]]], size(pv.bufcur_sp, 1))
    end
    #    return isvelocitystate(self) ? pv.bufcur[ind.I[Nsum-Nsp:Nsum-1]...] : pv.bufcur[ind.I[end]]
end

function initcoef!(
    pv::PoissonVar{T,N,Nsp,Nv,StdPoisson2d},
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
    pv::PoissonVar{T,N,Nsp,Nv,StdAB2},
    self::AdvectionData{T,N},
) where {T,N,Nsp,Nv}
    #    st = getst(self)
    adv = self.adv
    sz = sizeall(adv)
    compute_charge!(self)
    compute_elfield!(self)

    vecbufc_v = (getcur_t(self) / step(adv.t_mesh[2])) * pv.t_elfield[1]
    vecbufc_sp = (-getcur_t(self) / step(adv.t_mesh[1])) * adv.t_mesh[2].points
    bufc_v = zeros(T, sz)
    bufc_sp = zeros(T, sz)
    for ind in CartesianIndices(sz)
        bufc_sp[ind] = vecbufc_sp[ind.I[2]]
        bufc_v[ind] = vecbufc_v[ind.I[1]]
    end
    pushfirst!(pv.t_bufc_sp,bufc_sp)
    pushfirst!(pv.t_bufc_v, bufc_v)

    if length(pv.t_bufc_sp) == 3
        deleteat!(pv.t_bufc_sp, 3)
        deleteat!(pv.t_bufc_v, 3)
    end
    fmrbuf = zeros(T, sz)
    for i=1:length(pv.t_bufc_sp)
        interpolate!(fmrbuf,pv.t_bufc_sp[i], ind -> (bufc_sp[ind], bufc_v[ind]), adv.t_interp )
        pv.t_bufc_sp[i] .= fmrbuf
        interpolate!(fmrbuf,pv.t_bufc_v[i], ind -> (bufc_sp[ind], bufc_v[ind]), adv.t_interp )
        pv.t_bufc_v[i] .= fmrbuf
    end


    if ismissing(pv.bufcur_sp)
        pv.bufcur_sp = zeros(T,sz)
        pv.bufcur_v = zeros(T,sz)
    end

    if length(pv.t_bufc_sp) == 2
        pv.bufcur_sp .= (3pv.t_bufc_sp[1] - pv.t_bufc_sp[2])/2
        pv.bufcur_v .= (3pv.t_bufc_v[1] - pv.t_bufc_v[2])/2
    else
        pv.bufcur_sp .= pv.t_bufc_sp[1]
        pv.bufcur_v .= pv.t_bufc_v[1]
    end

end
@inline function getalpha(
    pv::PoissonVar{T,N,Nsp,Nv,StdAB2},
    self::AdvectionData{T},
    indext,
    ind,
) where {T,N,Nsp,Nv}
    @assert Nsp == Nv == 1 "Nsp=$Nsp Nv=$Nv they must be equal to one"
    #    @show ind.I
    return (pv.bufcur_sp[ind], pv.bufcur_v[ind])
end
@inline function getalpha(
    pv::PoissonVar{T,N,Nsp,Nv,StdPoisson2d},
    self::AdvectionData{T},
    indext,
    ind,
) where {T,N,Nsp,Nv}
    @assert Nsp == Nv == 1 "Nsp=$Nsp Nv=$Nv they must be equal to one"
    #    @show ind.I
    return (pv.bufcur_sp[ind.I[2]], pv.bufcur_v[ind.I[1]])
end
@inline function getalpha(
    pv::PoissonVar{T,N,Nsp,Nv,StdPoisson2dTry},
    self::AdvectionData{T},
    indext,
    ind,
) where {T,N,Nsp,Nv}
    @assert Nsp == Nv == 1 "Nsp=$Nsp Nv=$Nv they must be equal to one"
    #    @show ind.I
    return (pv.bufcur_sp[ind.I[2]], pv.bufcur_v[ind.I[1]])
end
function initcoef!(
    pv::PoissonVar{T,N,Nsp,Nv,StdPoisson2dTry},
    self::AdvectionData{T,N},
) where {T,N,Nsp,Nv}
    #    st = getst(self)
    adv = self.adv
    compute_charge!(self)
    compute_elfield!(self)

    bufc_v = (getcur_t(self) / step(adv.t_mesh[2])) * pv.t_elfield[1]

    bufc_sp = (-getcur_t(self) / step(adv.t_mesh[1])) * adv.t_mesh[2].points

    pv.bufcur_v = bufc_v
    pv.bufcur_sp = bufc_sp

    buf = zeros(T, sizeall(adv))

    vmin, vmax = extrema(pv.bufcur_v)
    spmin, spmax = extrema(pv.bufcur_sp)
    
    cnvmin(x) = x < 0 ? Int(-floor(x)) : 0
    cnvmax(x) = x > 0 ? Int(ceil(x)) : 0

    decbegin = cnvmin.((spmin,vmin))
    decend = cnvmax.((spmax, vmax))

    @show decbegin, decend, typeof(decbegin), typeof(decend)

   interpolatemod!(
        buf,
        pv.t_sp,
        ind -> getalpha(pv, self, 0, ind),
        adv.t_interp,
        T(length(adv.t_mesh[1])),
        decbegin,
        decend
    )

    pv.t_sp .= buf

    interpolatemod!(
        buf,
        pv.t_v,
        ind -> getalpha(pv, self, 0, ind),
        adv.t_interp,
        T(length(adv.t_mesh[2])),
        decbegin,
        decend
    )

    pv.t_v .= buf

end


function initcoef!(
    pv::PoissonVar{T,N,Nsp,Nv,StdOrder2_1},
    self::AdvectionData{T,N},
) where {T,N,Nsp,Nv}
    #    st = getst(self)
    adv = self.adv
    compute_charge!(self)
    compute_elfield!(self)
    if ismissing(pv.t_elfprec)
        elf = pv.t_elfield[1]
    else
        elf = (3pv.t_elfield[1] - pv.t_elfprec[1]) / 2
    end

    pv.t_elfprec = pv.t_elfield
    sz = sizeall(self.adv)
    if ismissing(pv.bufcur_v)
        pv.bufcur_v = zeros(T, sz)
        pv.bufcur_sp = zeros(T, sz)
    end

    bufc_v = (getcur_t(self) / step(adv.t_mesh[2])) * elf
    bufc_sp = (-getcur_t(self) / 2step(adv.t_mesh[1])) * adv.t_mesh[2].points

    interp = adv.t_interp[1]

    lgn = zeros(T, sz[1])

    for i = 1:sz[2]
        interpolate!(lgn, bufc_v, bufc_sp[i], interp)
        pv.bufcur_v[:, i] .= lgn
    end

    for i = 1:sz[1], j = 1:sz[2]
        pv.bufcur_sp[i, j] =
            (-getcur_t(self) / 2step(adv.t_mesh[1])) *
            (2adv.t_mesh[2].points[j] + pv.bufcur_v[i, j])
    end
    missing
end
function initcoef!(
    pv::PoissonVar{T,N,Nsp,Nv,StdOrder2_2},
    self::AdvectionData{T,N},
) where {T,N,Nsp,Nv}
    #    st = getst(self)
    adv = self.adv
    compute_charge!(self)
    compute_elfield!(self)

    sz = sizeall(self.adv)

    bufc_v = (getcur_t(self) / 2step(adv.t_mesh[2])) * pv.t_elfield[1]
    bufc_sp = (-getcur_t(self) / 2step(adv.t_mesh[1])) * adv.t_mesh[2].points

    result = zeros(T, sz)
    interpolate!(
        result,
        self.data,
        ind -> (bufc_sp[ind.I[2]], bufc_v[ind.I[1]]),
        adv.t_interp,
    )

    rho = zeros(T, sz[1])

    elf = zeros(T, sz[1])

    compute_charge!(rho, (adv.t_mesh[2],), result)
    compute_elfield!(elf, adv.t_mesh[1], rho)
    pv.bufcur_v = (getcur_t(self) / step(adv.t_mesh[2])) * elf
    pv.bufcur_sp = [
        (-getcur_t(self) / 2step(adv.t_mesh[1])) *
        (2 * adv.t_mesh[2].points[j] + pv.bufcur_v[i]) for i = 1:sz[1], j = 1:sz[2]
    ]
    missing
end
function firstinitdata(
    pv::PoissonVar{T,N,Nsp,Nv,StdAB,typeadd},
    self::AdvectionData,
) where {T,N,Nsp,Nv,typeadd}
    adv = self.adv
    pc = pv.pc
    sz = sizeall(adv)
    tabref = [zeros(T, sz) for i = 1:2typeadd+1]
    buf = zeros(T, sz)
    revtabref = reverse(tabref)
    ref_i = typeadd + 1

    copyto!(tabref[ref_i], self.data)

    rho = zeros(T, sz[1])

    elf = zeros(T, sz[1])


    for ord = 0:typeadd-1
        tabsens = ord == 0 ? [-1] : [1, -1]
        for sens in tabsens
            tab = sens == -1 ? revtabref : tabref
            borne = sens == -1 ? ord + 1 : ord
            for i = 1:borne
                buf .= sum([c(pc.abcoef, i, ord) * tab[ref_i+borne-1-i] for i = 0:ord])
                compute_charge!(rho, (adv.t_mesh[2],), buf)
                compute_elfield!(elf, adv.t_mesh[1], rho)
                bufc_v = (sens * getcur_t(self) / step(adv.t_mesh[2])) * elf
                bufc_sp =
                    (-sens * getcur_t(self) / step(adv.t_mesh[1])) * adv.t_mesh[2].points
                interpolate!(
                    tab[ref_i+borne],
                    tab[ref_i+borne-1],
                    ind -> (bufc_sp[ind.I[2]], bufc_v[ind.I[1]]),
                    adv.t_interp,
                )
            end
        end
    end
    ret = revtabref[ref_i:end]
    @show length(ret), typeadd
    return ret
end
function initcoef!(
    pv::PoissonVar{T,N,Nsp,Nv,StdAB,typeadd},
    self::AdvectionData{T,N},
) where {T,N,Nsp,Nv,typeadd}
    #    st = getst(self)
    pc = pv.pc
    adv = self.adv
    sz = sizeall(adv)
    # NOTE !!!! a ameliorer pour la gestion de la memoire !!!!
    if length(pv.t_data) == 0
        pv.t_data = firstinitdata(pv, self)
    else
        pv.t_data = circshift(pv.t_data, 1)
        pv.t_data[1] = copy(self.data)
    end

    ord = length(pv.t_data) - 1

    locdata = sum([c(pc.abcoef, i, ord) * pv.t_data[i+1] for i = 0:ord])

    rho = zeros(T, sz[1])

    elf = zeros(T, sz[1])

    compute_charge!(rho, (adv.t_mesh[2],), locdata)
    compute_elfield!(elf, adv.t_mesh[1], rho)
    pv.bufcur_v = (getcur_t(self) / step(adv.t_mesh[2])) * elf
    pv.bufcur_sp = (-getcur_t(self) / step(adv.t_mesh[1])) * adv.t_mesh[2].points
    missing
end

function caldata!( bufdata::Array{T,2}, adv::Advection, buf_sp::Array{T,2}, buf_v::Array{T,2}) where{T}
    coef = 1 / sqrt(2T(pi))
    for ind in CartesianIndices(sz)
        sp = stdtomesh(adv.t_mesh[1],buf_sp[ind])
        v = stdtomesh(adv.t_mesh[2],buf_v[ind])
        bufdata[ind] = coef * exp(T(-0.5) * v^2) * (1 + T(big"0.001") * cos( sp / 2))
    end
end

# function firstinitdata(
#     pv::PoissonVar{T,N,Nsp,Nv,StdABNew,typeadd},
#     self::AdvectionData,
# ) where {T,N,Nsp,Nv,typeadd}
#     adv = self.adv
#     pc = pv.pc
#     sz = sizeall(adv)
#     tabref_bufc_sp = [zeros(T, sz) for i = 1:2typeadd+1]
#     tabref_bufc_v = [zeros(T, sz) for i = 1:2typeadd+1]
#     buf = zeros(T, sz)
#     buf_sp = zeros(T, sz)
#     buf_v = zeros(T, sz)
#     revtabref_bufc_sp = reverse(tabref_bufc_sp)
#     revtabref_bufc_v = reverse(tabref_bufc_v)
#     ref_i = typeadd + 1

#     copyto!(tabref_bufc_sp[ref_i], pv.bufc_sp)
#     copyto!(tabref_bufc_v[ref_i], pv.bufc_v)

#     rho = zeros(T, sz[1])

#     elf = zeros(T, sz[1])


#     for ord = 0:typeadd-1
#         tabsens = ord == 0 ? [-1] : [1, -1]
#         for sens in tabsens
#             tab_sp = sens == -1 ? revtabref_bufc_sp : tabref_bufc_sp
#             tab_v = sens == -1 ? revtabref_buf_c_v : tabref__bufc_v
#             borne = sens == -1 ? ord + 1 : ord
#             for i = 1:borne
#                 bufc_sp .= sum([c(pc.abcoef, i, ord) * tab_sp[ref_i+borne-1-i] for i = 0:ord])
#                 buf_v .= sum([c(pc.abcoef, i, ord) * tab_v[ref_i+borne-1-i] for i = 0:ord])
#                 caldata!(buf, adv, buf_sp, buf_v)
#                 compute_charge!(rho, (adv.t_mesh[2],), buf)
#                 compute_elfield!(elf, adv.t_mesh[1], rho)
#                 bufc_v = (sens * getcur_t(self) / step(adv.t_mesh[2])) * elf
#                 bufc_sp =
#                     (-sens * getcur_t(self) / step(adv.t_mesh[1])) * adv.t_mesh[2].points
#                 # interpolate!(
#                 #     tab[ref_i+borne],
#                 #     tab[ref_i+borne-1],
#                 #     ind -> (bufc_sp[ind.I[2]], bufc_v[ind.I[1]]),
#                 #     adv.t_interp,
#                 # )
#                 vmin, vmax = extrema(bufc_v)
#                 spmin, spmax = extrema(bufc_sp)
                
#                 cnvmin(x) = x < 0 ? Int(-floor(x)) : 0
#                 cnvmax(x) = x > 0 ? Int(ceil(x)) : 0
            
#                 decbegin = cnvmin.((spmin,vmin))
#                 decend = cnvmax.((spmax, vmax))
            
#                 @show decbegin, decend, typeof(decbegin), typeof(decend)
#                 interpolatemod!(
#                     tab_sp[ref_i+borne],
#                     tab_sp[ref_i+borne-1],
#                     ind -> (bufc_sp[ind.I[2]], bufc_v[ind.I[1]]),
#                     adv.t_interp,
#                     T(sz[1]),
#                     decbegin,
#                     decend
#                 )
#                 interpolatemod!(
#                     tab_v[ref_i+borne],
#                     tab_v[ref_i+borne-1],
#                     ind -> (bufc_sp[ind.I[2]], bufc_v[ind.I[1]]),
#                     adv.t_interp,
#                     T(sz[2]),
#                     decbegin,
#                     decend
#                 )
#             end
#         end
#     end
#     ret_sp = revtabref_sp[ref_i:end]
#     ret_v = revtabref_v[ref_i:end]
#     @show length(ret_sp), typeadd
#     return ret_sp, ret_v
# end


# function initcoef!(
#     pv::PoissonVar{T,N,Nsp,Nv,StdABNew,typeadd},
#     self::AdvectionData{T,N},
# ) where {T,N,Nsp,Nv,typeadd}
#     #    st = getst(self)
#     pc = pv.pc
#     adv = self.adv
#     sz = sizeall(adv)
#     # NOTE !!!! a ameliorer pour la gestion de la memoire !!!!
#     if length(pv.t_sp) == 0
#         pv.t_sp, pv.t_v = firstinitdata(pv, self)
#     else
#         pv.t_sp = circshift(pv.t_sp, 1)
#         pv.t_v = circshift(pv.t_v, 1)
#         pv.t_data[1] = copy(self.data)
#     end

#     ord = length(pv.t_data) - 1

#     locdata = sum([c(pc.abcoef, i, ord) * pv.t_data[i+1] for i = 0:ord])

#     rho = zeros(T, sz[1])

#     elf = zeros(T, sz[1])

#     compute_charge!(rho, (adv.t_mesh[2],), locdata)
#     compute_elfield!(elf, adv.t_mesh[1], rho)
#     pv.bufcur_v = (getcur_t(self) / step(adv.t_mesh[2])) * elf
#     pv.bufcur_sp = (-getcur_t(self) / step(adv.t_mesh[1])) * adv.t_mesh[2].points
#     missing
# end

@inline function getalpha(
    pv::PoissonVar{T,N,Nsp,Nv,StdAB},
    self::AdvectionData{T},
    indext,
    ind,
) where {T,N,Nsp,Nv}
    @assert Nsp == Nv == 1 "Nsp=$Nsp Nv=$Nv they must be equal to one"
    #    @show ind.I
    return (pv.bufcur_sp[ind.I[2]], pv.bufcur_v[ind.I[1]])
end
@inline function getalpha(
    pv::PoissonVar{T,N,Nsp,Nv,StdOrder2_2},
    self::AdvectionData{T},
    indext,
    ind,
) where {T,N,Nsp,Nv}
    @assert Nsp == Nv == 1 "Nsp=$Nsp Nv=$Nv they must be equal to one"
    #    @show ind.I
    return (pv.bufcur_sp[ind], pv.bufcur_v[ind.I[1]])
end
@inline function getalpha(
    pv::PoissonVar{T,N,Nsp,Nv,StdOrder2_1},
    self::AdvectionData{T},
    indext,
    ind,
) where {T,N,Nsp,Nv}
    @assert Nsp == Nv == 1 "Nsp=$Nsp Nv=$Nv they must be equal to one"
    #    @show ind.I
    return (pv.bufcur_sp[ind], pv.bufcur_v[ind])
end
