
@enum TypePoisson StdPoisson = 1 StdPoisson2d = 2 StdOrder2_1 = 3 StdOrder2_2 = 4 StdAB = 5 StdPoisson2dTry =
    6 StdABNew = 7 StdAB2 = 8 StdRK4 = 9 StdABinit = 10 StdABp = 11

function _get_fctv_k(adv::Advection{T,N,timeopt}) where {T,N,timeopt}
    fct_k(v) = im / sum(v .^ 2)
    Nsp = div(N, 2)
    v_k = vec_k_fft.(adv.t_mesh[1:Nsp])
    sz = length.(adv.t_mesh[1:Nsp])
    fctv_k_gen = fct_k.(collect(Iterators.product(v_k...)))
    fctv_k_gen[1] = 0
    return ntuple(x -> reshape(v_k[x], tupleshape(x, Nsp, sz[x])) .* fctv_k_gen, Nsp)
end

using SHA
cod(a::Array) = bytes2hex(sha256(reinterpret(UInt8, collect(Iterators.flatten(a)))))



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

struct BufcT{T}
    bufc_sp::Array{T,2}
    bufc_v::Array{T,2}
    t::T
end


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
    t_bufc::Vector{BufcT{T}}
    bufcur_sp::Any
    bufcur_v::Any
    tupleind::Any
    ordcur::Int
    indord::Int
    t_bufc_sp::Any
    t_bufc_v::Any
    flend::Bool
    init_data::Vector{Array{T,N}}

    function PoissonVar(
        pc::PoissonConst{T,N,Nsp,Nv,type,typeadd},
    ) where {T,N,Nsp,Nv,type,typeadd}
        sz = length.(pc.adv.t_mesh[1:Nsp])
        rho = Array{T,Nsp}(undef, sz)
        sz = sizeall(pc.adv)
        t_data = []
        t_bufc = []
        if type == StdPoisson2dTry
            t_sp = zeros(T, sz)
            t_v = zeros(T, sz)
            for ind in CartesianIndices(sz)
                t_sp[ind] = ind.I[1] - 1
                t_v[ind] = ind.I[2] - 1
            end
        else
            t_sp = zeros(T, ntuple(x -> 1, N))
            t_v = zeros(T, ntuple(x -> 1, N))
        end
        return new{T,N,Nsp,Nv,type,typeadd}(
            pc,
            rho,
            missing,
            missing,
            t_data,
            t_sp,
            t_v,
            t_bufc,
            missing,
            missing,
            missing,
            0,
            0,
            [],
            [],
            false,
            map(x -> zeros(T, sz), 1:(typeadd-1)),
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

    decbegin = cnvmin.((spmin, vmin))
    decend = cnvmax.((spmax, vmax))

    @show decbegin, decend, typeof(decbegin), typeof(decend)

    interpolatemod!(
        buf,
        pv.t_sp,
        ind -> getalpha(pv, self, 0, ind),
        adv.t_interp,
        T(length(adv.t_mesh[1])),
        decbegin,
        decend,
    )

    pv.t_sp .= buf

    interpolatemod!(
        buf,
        pv.t_v,
        ind -> getalpha(pv, self, 0, ind),
        adv.t_interp,
        T(length(adv.t_mesh[2])),
        decbegin,
        decend,
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
function initcoef!(
    pv::PoissonVar{T,N,Nsp,Nv,StdAB2},
    self::AdvectionData{T,N},
) where {T,N,Nsp,Nv}
    #    st = getst(self)
    adv = self.adv
    sz = sizeall(adv)
    compute_charge!(self)
    compute_elfield!(self)
    dt = getcur_t(self)
    @show "initcoef AB2", self.time_cur, dt

    vecbufc_v = (getcur_t(self) / step(adv.t_mesh[2])) * pv.t_elfield[1]
    vecbufc_sp = (-getcur_t(self) / step(adv.t_mesh[1])) * adv.t_mesh[2].points

    #    @show vecbufc_v
    @show "AB2", cod(vecbufc_v)

    bufc_v = zeros(T, sz)
    bufc_sp = zeros(T, sz)
    for ind in CartesianIndices(sz)
        bufc_sp[ind] = vecbufc_sp[ind.I[2]]
        bufc_v[ind] = vecbufc_v[ind.I[1]]
    end
    if length(pv.t_bufc_sp) > 0
        @show "TR1", sum(pv.t_bufc_sp[1])
    end

    pushfirst!(pv.t_bufc_sp, bufc_sp)
    pushfirst!(pv.t_bufc_v, bufc_v)

    if length(pv.t_bufc_sp) > 1
        @show "TR2", sum(pv.t_bufc_sp[2])
    end

    if length(pv.t_bufc_sp) == 3
        deleteat!(pv.t_bufc_sp, 3)
        deleteat!(pv.t_bufc_v, 3)
    end
    fmrbuf = zeros(T, sz)
    fmr_sp = copy(bufc_sp)
    fmr_v = copy(bufc_v)
    for i = 1:length(pv.t_bufc_sp)
        @show "TR5", i, pv.t_bufc_sp[i][2:2], fmr_sp[2:2], fmr_v[2:2]
        interpolate!(
            fmrbuf,
            pv.t_bufc_sp[i],
            #            ind -> (bufc_sp[ind], bufc_v[ind]),
            ind -> (fmr_sp[ind], fmr_v[ind]),
            adv.t_interp,
        )
        pv.t_bufc_sp[i] .= fmrbuf
        interpolate!(fmrbuf, pv.t_bufc_v[i], ind -> (fmr_sp[ind], fmr_v[ind]), adv.t_interp)
        pv.t_bufc_v[i] .= fmrbuf
    end
    @show sum(pv.t_bufc_sp[1])
    @show sum(pv.t_bufc_v[1])
    if length(pv.t_bufc_sp) > 1
        @show sum(pv.t_bufc_sp[2])
        @show sum(pv.t_bufc_v[2])
    end


    if ismissing(pv.bufcur_sp)
        pv.bufcur_sp = zeros(T, sz)
        pv.bufcur_v = zeros(T, sz)
    end

    if length(pv.t_bufc_sp) == 2
        pv.bufcur_sp .= (3pv.t_bufc_sp[1] - pv.t_bufc_sp[2]) / 2
        pv.bufcur_v .= (3pv.t_bufc_v[1] - pv.t_bufc_v[2]) / 2
    else
        pv.bufcur_sp .= pv.t_bufc_sp[1]
        pv.bufcur_v .= pv.t_bufc_v[1]
    end
    @show "TR6", sum(pv.bufcur_sp), sum(pv.bufcur_v)
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


# function getcur_t(
#     self::AdvectionData{T,2},
#     pv::PoissonVar{T,2,1,1,StdAB,order},
# ) where {T,order}
#     c_ord = pv.ordcur - (pv.indord == 1 && pv.ordcur != 1)
#     ret = getcur_t(self.adv, self.state_gen)
#     if c_ord == order
#         return ret
#     else
#         return ret / prod(c_ord+1:order)
#     end
# end
# function retns(self::AdvectionData{T,2}, pv::PoissonVar{T,2,1,1,StdAB,order},) where {T,order}
#     @show "retns", pv.ordcur, order
#     return pv.ordcur != order
# end

function sum_ab!(
    pv::PoissonVar{T,2,1,1,StdAB,order},
    self::AdvectionData{T,2},
) where {T,order}
    dt = getcur_t(self)
    if pv.ordcur == pv.indord == order
        coefs = map(i -> dt * c(pv.pc.abcoef, i, order), 1:order)
    else
        tabv = map(i -> pv.t_bufc[i].t, 1:pv.ordcur)
        tabv .-= tabv[1]
        coefs = map(i -> _c(dt, tabv, i, pv.ordcur), 1:pv.ordcur)
    end
    #    coefs = map( i -> c(pv.pc.abcoef, i, pv.ordcur), 1:pv.ordcur)
    ords = map(i -> pv.t_bufc[i].order, 1:pv.ordcur)
    @show " sum_ab ", pv.ordcur, pv.indord, coefs, ords
    pv.bufcur_sp .= sum(map(i -> coefs[i] * pv.t_bufc[i].bufc_sp, 1:pv.ordcur))
    pv.bufcur_v .= sum(map(i -> coefs[i] * pv.t_bufc[i].bufc_v, 1:pv.ordcur))
end

function deletebufc!(pv::PoissonVar{T,2,1,1,StdAB,order}) where {T,order}
    mindiff = Inf
    ind = -1
    for i = 2:length(pv.t_bufc)-1
        dt = pv.t_bufc[i-1].t - pv.t_bufc[i+1].t
        if dt < mindiff && 0 < pv.t_bufc[i].order < pv.ordcur
            mindiff = dt
            ind = i
        end
    end
    if ind > 0
        println("delete order=$(pv.t_bufc[ind].order)")
        deleteat!(pv.t_bufc, ind)
    end
end
function interpbufc!(
    pv::PoissonVar{T,2,1,1,type,typeadd},
    self::AdvectionData{T,2},
    vec_sp::Vector{T},
    vec_v::Vector{T},
) where {T,type,typeadd}
    adv = self.adv
    sz = sizeall(adv)
    fmrbuf = zeros(T, sz)
    for bct in pv.t_bufc
        interpolate!(
            fmrbuf,
            bct.bufc_sp,
            ind -> (vec_sp[ind.I[2]], vec_v[ind.I[1]]),
            adv.t_interp,
        )
        bct.bufc_sp .= fmrbuf
        interpolate!(
            fmrbuf,
            bct.bufc_v,
            ind -> (vec_sp[ind.I[2]], vec_v[ind.I[1]]),
            adv.t_interp,
        )
        bct.bufc_v .= fmrbuf
    end
end
function autointerp!(
    to::NTuple{2,Array{T,2}},
    from::NTuple{2,Array{T,2}},
    nb::Int,
    t_interp::Vector{I},
)   where{T, I<:AbstractInterpolation}

    fmr = copy.(from)

    for i=1:nb
        interpolate!(to[1], from[1], ind ->(fmr[1][ind], fmr[2][ind]), t_interp)
        interpolate!(to[2], from[2], ind ->(fmr[1][ind], fmr[2][ind]), t_interp)
        if i != nb
            fmr[1] .= to[1]
            fmr[2] .= to[2]
        end
    end
end

function getenergyall(adv, buf::Array{T,2}) where {T}
    sz = size(buf, 1)
    rho = zeros(T, sz)
    elf = zeros(T, sz)
    compute_charge!(rho, (adv.t_mesh[2],), buf)
    compute_elfield!(elf, adv.t_mesh[1], rho)

    elenergy = compute_ee((adv.t_mesh[1],), (elf,))
    kinenergy = compute_ke((adv.t_mesh[1],), (adv.t_mesh[2],), buf)
    energyall = elenergy + kinenergy
    return energyall
end

function interpbufc!(
    pv::PoissonVar{T,2,1,1,type,typeadd},
    self::AdvectionData{T,2},
    bufc_sp::Array{T,2},
    bufc_v::Array{T,2},
) where {T,type,typeadd}
    adv = self.adv
    sz = sizeall(adv)
    fmrbuf = zeros(T, sz)
    ind = 1
    for bct in pv.t_bufc
        @show "TR3", ind, sum(bct.bufc_sp)
        @show "TR5", ind, bct.bufc_sp[2:2], bufc_sp[2:2], bufc_v[2:2]
        interpolate!(fmrbuf, bct.bufc_sp, ind -> (bufc_sp[ind], bufc_v[ind]), adv.t_interp)
        bct.bufc_sp .= fmrbuf
        @show "TR4", ind, sum(bct.bufc_sp)
        interpolate!(fmrbuf, bct.bufc_v, ind -> (bufc_sp[ind], bufc_v[ind]), adv.t_interp)
        bct.bufc_v .= fmrbuf
        ind += 1
    end
end
struct BufcFmr{T}
    oriplus_sp::Array{T,2}
    oriplus_v::Array{T,2}
    orimoins_sp::Array{T,2}
    orimoins_v::Array{T,2}
    cur_sp::Array{T,2}
    cur_v::Array{T,2}
    function BufcFmr(T::DataType, sz::Tuple{Int,Int})
        return new{T}(
            zeros(T, sz),
            zeros(T, sz),
            zeros(T, sz),
            zeros(T, sz),
            zeros(T, sz),
            zeros(T, sz),
        )
    end
end
function decall!(
    tab::Vector{BufcFmr{T}},
    adv::Advection,
    sens::Int,
    ref_i::Int,
    nb::Int,
    dec::Int = 0,
) where {T}
    for j = 1:nb
        ptab = tab[ref_i-nb+j+dec]

        dec_sp = sens == -1 ? ptab.orimoins_sp : ptab.oriplus_sp
        dec_v = sens == -1 ? ptab.orimoins_v : ptab.oriplus_v

        ori_sp, ori_v = if j == nb
            dec_sp, dec_v
        else
            copy(ptab.cur_sp), copy(ptab.cur_v)
        end

        interpolate!(ptab.cur_sp, ori_sp, ind -> (dec_sp[ind], dec_v[ind]), adv.t_interp)
        interpolate!(ptab.cur_v, ori_v, ind -> (dec_sp[ind], dec_v[ind]), adv.t_interp)
    end
end

function initfirst(
    pv::PoissonVar{T,2,1,1,StdAB,order},
    self::AdvectionData,
    bufcc_sp::Array{T,2},
    bufcc_v::Array{T,2},
) where {T,order}
    adv = self.adv
    pc = pv.pc
    sz = sizeall(adv)
    tabref = [BufcFmr(T, sz) for i = 1:2order-1]
    buf = copy(self.data)
    buf_sp = zeros(T, sz)
    buf_v = zeros(T, sz)
    revtabref = reverse(tabref)
    ref_i = order

    dt = getcur_t(self)

    ptab = tabref[ref_i]
    copyto!(tabref[ref_i].oriplus_sp, bufcc_sp)
    copyto!(tabref[ref_i].oriplus_sp, bufcc_v)

    interpolate!(
        buf_sp,
        ptab.oriplus_sp,
        ind -> (ptab.oriplus_sp[ind], ptab.oriplus_v[ind]),
        adv.t_interp,
    )
    interpolate!(
        buf_v,
        ptab.oriplus_v,
        ind -> (ptab.oriplus_sp[ind], ptab.oriplus_v[ind]),
        adv.t_interp,
    )


    resinv = getinverse((buf_sp, buf_v), adv.t_interp)
    tabref[ref_i].orimoins_sp .= resinv[1]
    tabref[ref_i].orimoins_v .= resinv[2]

    # tabref[ref_i].orimoins_sp .= - tabref[ref_i].oriplus_sp
    # tabref[ref_i].orimoins_v .= - tabref[ref_i].oriplus_v

    rho = zeros(T, sz[1])

    elf = zeros(T, sz[1])


    refen = getenergyall(self.adv, self.data)

    @show refen


    for ord = 1:order-1
        tabsens = ord == 1 ? [-1] : [1, -1]
        for sens in tabsens
            tab = sens == -1 ? revtabref : tabref
            borne = sens == -1 ? ord : ord - 1
            for j = 1:ord-1
                decall!(tab, adv, sens, ref_i, j)
                # for k = 1:j
                #     ptab = tab[ref_i-ord+k]

                #     dec_sp = sens == -1 ? ptab.orimoins_sp : ptab.oriplus_sp
                #     dec_v = sens == -1 ? ptab.orimoins_v : ptab.oriplus_v

                #     ori_sp, ori_v = if k == j
                #         dec_sp, dec_v
                #     else
                #         copy(ptab.cur_sp), copy(ptab.cur_v)
                #     end
                #     interpolate!( ptab.cur_sp, ori_sp, ind -> ( dec_sp[ind], dec_v[ind] ), adv.t_interp)
                #     interpolate!( ptab.cur_v, ori_v, ind -> ( dec_sp[ind], dec_v[ind] ), adv.t_interp)
                # end
            end
            for i = 1:borne
                decall!(tab, adv, sens, ref_i, ord, i - 1)
                # for j=1:ord
                #     ptab = tab[ref_i-ord+i+j-1]

                #     dec_sp = sens == -1 ? ptab.orimoins_sp : ptab.oriplus_sp
                #     dec_v = sens == -1 ? ptab.orimoins_v : ptab.oriplus_v

                #     ori_sp, ori_v = if j == ord
                #         dec_sp, dec_v
                #     else
                #         copy(ptab.cur_sp), copy(ptab.cur_v)
                #     end

                #     interpolate!( ptab.cur_sp, ori_sp, ind -> ( dec_sp[ind], dec_v[ind] ), adv.t_interp)
                #     interpolate!( ptab.cur_v, ori_v, ind -> ( dec_sp[ind], dec_v[ind] ), adv.t_interp)
                # end

                buf_sp .= sum([c(pc.abcoef, j, ord) * tab[ref_i+i-j].cur_sp for j = 1:ord])
                buf_v .= sum([c(pc.abcoef, j, ord) * tab[ref_i+i-j].cur_v for j = 1:ord])
                interpolate!(buf, copy(buf), ind -> (buf_sp[ind], buf_v[ind]), adv.t_interp)
                @show refen
                diffen = abs(refen - getenergyall(self.adv, buf))
                @show sens, ord, i, diffen
                compute_charge!(rho, (adv.t_mesh[2],), buf)
                compute_elfield!(elf, adv.t_mesh[1], rho)
                ptab = tab[ref_i+i]
                vec_v = (getcur_t(self) / step(adv.t_mesh[2])) * elf
                vec_sp = (-getcur_t(self) / step(adv.t_mesh[1])) * adv.t_mesh[2].points
                for ind in CartesianIndices(sz)
                    ptab.oriplus_sp[ind] = vec_sp[ind.I[2]]
                    ptab.oriplus_v[ind] = vec_v[ind.I[1]]
                end
                # ptab.orimoins_sp .= - ptab.oriplus_sp
                # ptab.orimoins_v .= - ptab.oriplus_v
                #                 ptab
                interpolate!(
                    buf_sp,
                    ptab.oriplus_sp,
                    ind -> (ptab.oriplus_sp[ind], ptab.oriplus_v[ind]),
                    adv.t_interp,
                )
                interpolate!(
                    buf_v,
                    ptab.oriplus_v,
                    ind -> (ptab.oriplus_sp[ind], ptab.oriplus_v[ind]),
                    adv.t_interp,
                )

                #                resinv = getinverse((ptab.oriplus_sp,ptab.oriplus_v), adv.t_interp)
                resinv = getinverse((buf_sp, buf_v), adv.t_interp)
                ptab.orimoins_sp .= resinv[1]
                ptab.orimoins_v .= resinv[2]
                # #                interpolate!(ptab.orimoins_sp, resinv[1], ind ->(resinv[1][ind],resinv[2][ind]), adv.t_interp)
                #                 # interpolate!(ptab.orimoins_v, resinv[2], ind ->(resinv[1][ind],resinv[2][ind]), adv.t_interp)
                #                 interpolate!(ptab.orimoins_sp, resinv[1], ind ->(ptab.oriplus_sp[ind],ptab.oriplus_v[ind]), adv.t_interp)
                #                 interpolate!(ptab.orimoins_v, resinv[2], ind ->(ptab.oriplus_sp[ind],ptab.oriplus_v[ind]), adv.t_interp)
            end
        end
    end
    for i = 1:order-1
        decall!(tabref, adv, 1, ref_i, i)
    end
    return tabref
end

function initcoef!(
    pv::PoissonVar{T,2,1,1,StdAB,order},
    self::AdvectionData{T,2},
) where {T,order}
    adv = self.adv
    sz = sizeall(adv)

    compute_charge!(self)
    compute_elfield!(self)
    dt = getcur_t(self)

    vecbufc_v = (dt / step(adv.t_mesh[2])) * pv.t_elfield[1]
    vecbufc_sp = (-dt / step(adv.t_mesh[1])) * adv.t_mesh[2].points
    bufc_v = zeros(T, sz)
    bufc_sp = zeros(T, sz)
    for ind in CartesianIndices(sz)
        bufc_sp[ind] = vecbufc_sp[ind.I[2]]
        bufc_v[ind] = vecbufc_v[ind.I[1]]
    end
    if ismissing(pv.bufcur_sp)
        pv.bufcur_sp = zeros(T, sz)
        pv.bufcur_v = zeros(T, sz)
        tabref = initfirst(pv, self, copy(bufc_sp), copy(bufc_v))
        for i = 1:order-1
            ptab = tabref[i]
            pushfirst!(pv.t_bufc, BufcT(ptab.cur_sp, ptab.cur_v, (i - order) * dt))
        end
    end

    pushfirst!(pv.t_bufc, BufcT(copy(bufc_sp), copy(bufc_v), self.time_cur))

    bufc_sp .= sum(map(i -> c(pv.pc.abcoef, i, order) * pv.t_bufc[i].bufc_sp, 1:order))
    bufc_v .= sum(map(i -> c(pv.pc.abcoef, i, order) * pv.t_bufc[i].bufc_v, 1:order))

    autointerp!((pv.bufcur_sp, pv.bufcur_v),(bufc_sp, bufc_v), order, adv.t_interp)

    deleteat!(pv.t_bufc, length(pv.t_bufc))

    interpbufc!(pv, self, pv.bufcur_sp, pv.bufcur_v)


    #    @show "TR6", sum(pv.bufcur_sp), sum(pv.bufcur_v)

end
@inline function getalpha(
    pv::PoissonVar{T,2,1,1,StdAB},
    self::AdvectionData{T},
    indext,
    ind,
) where {T}
    return (pv.bufcur_sp[ind], pv.bufcur_v[ind])
end
function initcoef!(
    pv::PoissonVar{T,2,1,1,StdABp,order},
    self::AdvectionData{T,2},
) where {T,order}
    adv = self.adv
    sz = sizeall(adv)

    compute_charge!(self)
    compute_elfield!(self)
    dt = getcur_t(self)

    vecbufc_v = (dt / step(adv.t_mesh[2])) * pv.t_elfield[1]
    vecbufc_sp = (-dt / step(adv.t_mesh[1])) * adv.t_mesh[2].points
    bufc_v = zeros(T, sz)
    bufc_sp = zeros(T, sz)
    for ind in CartesianIndices(sz)
        bufc_sp[ind] = vecbufc_sp[ind.I[2]]
        bufc_v[ind] = vecbufc_v[ind.I[1]]
    end
    if ismissing(pv.bufcur_sp)
        pv.bufcur_sp = zeros(T, sz)
        pv.bufcur_v = zeros(T, sz)
        fmr_v = zeros(T, sz)
        fmr_sp = zeros(T, sz)
    
        for indice = 1:order-1
            pushfirst!(pv.t_bufc, BufcT(copy(bufc_sp), copy(bufc_v), (indice - order) * dt))
            fmr_sp = sum(map(i -> c(pv.pc.abcoef, i, indice) * pv.t_bufc[i].bufc_sp, 1:indice))
            fmr_v = sum(map(i -> c(pv.pc.abcoef, i, indice) * pv.t_bufc[i].bufc_v, 1:indice))
            autointerp!((fmr_sp, fmr_v),(copy(fmr_sp), copy(fmr_v)), indice-1, adv.t_interp)
            interpbufc!(pv, self, fmr_sp, fmr_v)
        end
    end
        
    pushfirst!(pv.t_bufc, BufcT(copy(bufc_sp), copy(bufc_v), self.time_cur))

    bufc_sp .= sum(map(i -> c(pv.pc.abcoef, i, order) * pv.t_bufc[i].bufc_sp, 1:order))
    bufc_v .= sum(map(i -> c(pv.pc.abcoef, i, order) * pv.t_bufc[i].bufc_v, 1:order))

    autointerp!((pv.bufcur_sp, pv.bufcur_v),(bufc_sp, bufc_v), order-1, adv.t_interp)

    deleteat!(pv.t_bufc, length(pv.t_bufc))

    interpbufc!(pv, self, pv.bufcur_sp, pv.bufcur_v)


    #    @show "TR6", sum(pv.bufcur_sp), sum(pv.bufcur_v)

end
@inline function getalpha(
    pv::PoissonVar{T,2,1,1,StdABp},
    self::AdvectionData{T},
    indext,
    ind,
) where {T}
    return (pv.bufcur_sp[ind], pv.bufcur_v[ind])
end



function initcoef!(
    pv::PoissonVar{T,2,1,1,StdABinit,order},
    self::AdvectionData{T,2},
) where {T,order}
    adv = self.adv
    sz = sizeall(adv)

    dt = getcur_t(self)

    indice = Int(div(self.time_cur + dt / 2, dt)) + 1

    flag_zero = indice < order

    compute_charge!(self)
    compute_elfield!(self)

    vecbufc_v = (dt / step(adv.t_mesh[2])) * pv.t_elfield[1]
    vecbufc_sp = (-dt / step(adv.t_mesh[1])) * adv.t_mesh[2].points
    bufc_v = zeros(T, sz)
    bufc_sp = zeros(T, sz)
    for ind in CartesianIndices(sz)
        bufc_sp[ind] = vecbufc_sp[ind.I[2]]
        bufc_v[ind] = vecbufc_v[ind.I[1]]
    end
    if ismissing(pv.bufcur_sp)
        pv.bufcur_sp = zeros(T, sz)
        pv.bufcur_v = zeros(T, sz)
    end

    pushfirst!(pv.t_bufc, BufcT(copy(bufc_sp), copy(bufc_v), self.time_cur))

    #    interpbufc!(pv, self, bufc_sp, bufc_v)



    #    @show "TR6", sum(pv.bufcur_sp), sum(pv.bufcur_v)
    if flag_zero
        self.data .= pv.init_data[indice]
        #        fill!(pv.bufcur_sp,0)
        #        fill!(pv.bufcur_v,0)
        fmr_sp = sum(map(i -> c(pv.pc.abcoef, i, indice) * pv.t_bufc[i].bufc_sp, 1:indice))
        fmr_v = sum(map(i -> c(pv.pc.abcoef, i, indice) * pv.t_bufc[i].bufc_v, 1:indice))
        fmr3_sp = copy(fmr_sp)
        fmr3_v = copy(fmr_v)
        fmr2_sp = zeros(T, sz)
        fmr2_v = zeros(T, sz)
        # resprec = Inf
        # prec_sp = eps(T)
        # prec_v = eps(T)
        # i = 1
        for i = 1:order-1
            interpolate!(fmr2_sp, fmr_sp, ind -> (fmr3_sp[ind], fmr3_v[ind]), adv.t_interp)
            interpolate!(fmr2_v, fmr_v, ind -> (fmr3_sp[ind], fmr3_v[ind]), adv.t_interp)
            resnorm_sp = norm(fmr2_sp - fmr3_sp)
            resnorm_v = norm(fmr2_v - fmr3_v)
            @show i, order, resnorm_sp, resnorm_v
            # res = resnorm_sp/prec_sp + resnorm_v/prec_v
            # if res > resprec || i > order ^3
            #     break
            # end
            # resprec = res
            # prec_sp = resnorm_sp
            # prec_v = resnorm_v
            if i != order - 1
                fmr3_sp .= fmr2_sp
                fmr3_v .= fmr2_v
            end
        end
        interpbufc!(pv, self, fmr2_sp, fmr2_v)
    else
        fmr_sp = sum(map(i -> c(pv.pc.abcoef, i, order) * pv.t_bufc[i].bufc_sp, 1:order))
        fmr_v = sum(map(i -> c(pv.pc.abcoef, i, order) * pv.t_bufc[i].bufc_v, 1:order))
        fmr3_sp = copy(fmr_sp)
        fmr3_v = copy(fmr_v)
        fmr2_sp = zeros(T, sz)
        fmr2_v = zeros(T, sz)
        # resprec = Inf
        # prec_sp = eps(T)
        # prec_v = eps(T)
        # i= 1
        for i = 1:order-1
            interpolate!(fmr2_sp, fmr_sp, ind -> (fmr3_sp[ind], fmr3_v[ind]), adv.t_interp)
            interpolate!(fmr2_v, fmr_v, ind -> (fmr3_sp[ind], fmr3_v[ind]), adv.t_interp)
            resnorm_sp = norm(fmr2_sp - fmr3_sp)
            resnorm_v = norm(fmr2_v - fmr3_v)
            @show i, order, resnorm_sp, resnorm_v
            # res = resnorm_sp/prec_sp + resnorm_v/prec_v
            # if res > resprec || i > order ^3
            #     break
            # end
            # resprec = res
            # prec_sp = resnorm_sp
            # prec_v = resnorm_v
            if i != order - 1
                fmr3_sp .= fmr2_sp
                fmr3_v .= fmr2_v
            end
        end
        interpbufc!(pv, self, fmr2_sp, fmr2_v)
        pv.bufcur_sp .=
            sum(map(i -> c(pv.pc.abcoef, i, order) * pv.t_bufc[i].bufc_sp, 1:order))
        pv.bufcur_v .=
            sum(map(i -> c(pv.pc.abcoef, i, order) * pv.t_bufc[i].bufc_v, 1:order))
        deleteat!(pv.t_bufc, length(pv.t_bufc))
    end

end
@inline function getalpha(
    pv::PoissonVar{T,2,1,1,StdABinit},
    self::AdvectionData{T},
    indext,
    ind,
) where {T}
    return (pv.bufcur_sp[ind], pv.bufcur_v[ind])
end

function initcoef!(
    pv::PoissonVar{T,2,1,1,StdRK4,order},
    self::AdvectionData{T,2},
) where {T,order}
    adv = self.adv
    sz = sizeall(adv)

    rho = zeros(T, sz[1])

    elf = zeros(T, sz[1])


    compute_charge!(self)
    compute_elfield!(self)
    dt = getcur_t(self)
    dt2 = dt / 2

    vecbufc_v = (dt2 / step(adv.t_mesh[2])) * pv.t_elfield[1]
    vecbufc_sp = (-dt2 / step(adv.t_mesh[1])) * adv.t_mesh[2].points
    bufc_v1 = zeros(T, sz)
    bufc_sp1 = zeros(T, sz)
    bufc_v = zeros(T, sz)
    bufc_sp = zeros(T, sz)
    for ind in CartesianIndices(sz)
        bufc_sp[ind] = vecbufc_sp[ind.I[2]]
        bufc_v[ind] = vecbufc_v[ind.I[1]]
    end

    autointerp!((bufc_sp1, bufc_v1), (bufc_sp, bufc_v), 4, adv.t_interp)

    fmrdata = zeros(T, sz)
    interpolate!(fmrdata, self.data, ind -> (bufc_sp1[ind], bufc_v1[ind]), adv.t_interp)

    autointerp!((bufc_sp, bufc_v), (bufc_sp1, bufc_v1), 4, adv.t_interp)
    bufc_sp1 .= bufc_sp
    bufc_v1 .= bufc_v



    compute_charge!(rho, (adv.t_mesh[2],), fmrdata)
    compute_elfield!(elf, adv.t_mesh[1], rho)

    vecbufc_v = (dt2 / step(adv.t_mesh[2])) * elf
    vecbufc_sp = (-dt2 / step(adv.t_mesh[1])) * adv.t_mesh[2].points
    bufc_v2 = zeros(T, sz)
    bufc_sp2 = zeros(T, sz)
    for ind in CartesianIndices(sz)
        bufc_sp[ind] = vecbufc_sp[ind.I[2]]
        bufc_v[ind] = vecbufc_v[ind.I[1]]
    end

    autointerp!((bufc_sp2, bufc_v2), (bufc_sp, bufc_v), 4, adv.t_interp)

    # interpolate!(bufc_sp, bufc_sp1, ind ->(bufc_sp2[ind], bufc_v2[ind]), adv.t_interp)
    # interpolate!(bufc_v, bufc_v1, ind ->(bufc_sp2[ind], bufc_v2[ind]), adv.t_interp)
    # bufc_sp1 .= bufc_sp
    # bufc_v1 .= bufc_v

    interpolate!(fmrdata, self.data, ind -> (bufc_sp[ind], bufc_v[ind]), adv.t_interp)

    compute_charge!(rho, (adv.t_mesh[2],), fmrdata)
    compute_elfield!(elf, adv.t_mesh[1], rho)
    vecbufc_v = (dt2 / step(adv.t_mesh[2])) * elf
    vecbufc_sp = (-dt2 / step(adv.t_mesh[1])) * adv.t_mesh[2].points
    bufc_v3 = zeros(T, sz)
    bufc_sp3 = zeros(T, sz)
    for ind in CartesianIndices(sz)
        bufc_sp[ind] = vecbufc_sp[ind.I[2]]
        bufc_v[ind] = vecbufc_v[ind.I[1]]
    end

    autointerp!((bufc_sp3, bufc_v3), (bufc_sp, bufc_v), 4, adv.t_interp)

 
    interpolate!(fmrdata, self.data, ind -> (2bufc_sp3[ind], 2bufc_v3[ind]), adv.t_interp)

    compute_charge!(rho, (adv.t_mesh[2],), fmrdata)
    compute_elfield!(elf, adv.t_mesh[1], rho)

    vecbufc_v = (dt / step(adv.t_mesh[2])) * elf
    vecbufc_sp = (-dt / step(adv.t_mesh[1])) * adv.t_mesh[2].points
    bufc_v4 = zeros(T, sz)
    bufc_sp4 = zeros(T, sz)
    for ind in CartesianIndices(sz)
        bufc_sp4[ind] = vecbufc_sp[ind.I[2]]
        bufc_v4[ind] = vecbufc_v[ind.I[1]]
    end

    if ismissing(pv.bufcur_sp)
        pv.bufcur_sp = zeros(T, sz)
        pv.bufcur_v = zeros(T, sz)
    end

    pv.bufcur_v .= (2bufc_v1 + 4bufc_v2 + 4bufc_v3 + bufc_v4) / 6
    pv.bufcur_sp .= (2bufc_sp1 + 4bufc_sp2 + 4bufc_sp3 + bufc_sp4) / 6




end

@inline function getalpha(
    pv::PoissonVar{T,2,1,1,StdRK4},
    self::AdvectionData{T},
    indext,
    ind,
) where {T}
    return (pv.bufcur_sp[ind], pv.bufcur_v[ind])
end
function caldata!(
    bufdata::Array{T,2},
    adv::Advection,
    buf_sp::Array{T,2},
    buf_v::Array{T,2},
) where {T}
    coef = 1 / sqrt(2T(pi))
    for ind in CartesianIndices(sz)
        sp = stdtomesh(adv.t_mesh[1], buf_sp[ind])
        v = stdtomesh(adv.t_mesh[2], buf_v[ind])
        bufdata[ind] = coef * exp(T(-0.5) * v^2) * (1 + T(big"0.001") * cos(sp / 2))
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
