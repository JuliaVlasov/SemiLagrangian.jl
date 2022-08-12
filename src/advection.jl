
@enum TimeOptimization NoTimeOpt = 1 SimpleThreadsOpt = 2 SplitThreadsOpt = 3 MPIOpt = 4
@enum TimeAlgorithm NoTimeAlg = 1 ABTimeAlg_ip = 2 ABTimeAlg_new = 3 ABTimeAlg_init = 4

struct StateAdv{N}
    ind::Int          # indice
    perm::Vector{Int} # dimensions permutation
    invp::Vector{Int} # inverse of dimension permutation
    ndims::Int        # count of dimension
    stcoef::Int   # state_coef
    isconstdec::Bool #true if constant dec
    function StateAdv(ind, p, ndims, stc, isconst)
        return new{length(p)}(ind, p, invperm(p), ndims, stc, isconst)
    end
end

nosplit(dt::T) where {T} = dt * [1]
standardsplit(dt::T) where {T} = dt * [1, 1]
strangsplit(dt::T) where {T} = dt * [1 // 2, 1 // 1, 1 // 2]
magicsplit(dt::T) where {T} = [tan(dt / 2), sin(dt), tan(dt / 2)]
function triplejumpsplit(dt::T) where {T}
    c = T(2)^(1 // 3)
    c1 = 1 / (2(2 - c))
    c2 = (1 - c) / (2(2 - c))
    d1 = 1 / (2 - c)
    d2 = -c / (2 - c)
    return dt * [c1, d1, c2, d2, c2, d1, c1]
end

function order6split(dt::T) where {T}
    order6 = [
        0.0414649985182624,
        0.123229775946271,
        0.198128671918067,
        0.290553797799558,
        -0.0400061921041533,
        -0.127049212625417,
        0.0752539843015807,
        -0.246331761062075,
        -0.0115113874206879,
        0.357208872795928,
        0.23666992478693111,
        0.20477705429147008,
        0.23666992478693111,
        0.357208872795928,
        -0.0115113874206879,
        -0.246331761062075,
        0.0752539843015807,
        -0.127049212625417,
        -0.0400061921041533,
        0.290553797799558,
        0.198128671918067,
        0.123229775946271,
        0.0414649985182624,
    ]
    b = convert(Vector{T}, order6)
    b[11] = b[13] = T(1 // 2) - sum(b[1:2:9])
    b[12] = T(1 // 1) - 2sum(b[2:2:10])
    @assert isapprox(sum(b), T(2))
    return dt * b
end

function hamsplit_3_11(dt::T, fltrace = false) where {T<:Number}
    a = [
        big"0.168735950563437422448195173400884809990898960535052167820406",
        big"0.377851589220928303880768408101213978086828975884075336276904",
        big"-0.093175079568731452657927163004197576155455872838255008194621",
    ]
    b = [
        big"0.049086460976116245491441126327891629034401134948561362353776",
        big"0.264177609888976700200146195420764624729303307466109303375304",
        big"0.186735929134907054308412678251343746236295557585329334270919",
    ]
    c = [
        big"-0.0000697287150553050840997049705543302201654369851434306789330",
        big"-0.000625704827430047189169785837050053054170994982227957762075",
        big"-0.00221308512404532556162738103226349162177713815465419673226",
    ]
    d = [
        0,
        big"-2.91660045768984781641978316974322113103756038016905421426e-6",
        big"0.0000304848026170003878867997783299987163079240932335253763140",
    ]
    e = [0, 0, big"4.98554938787506812157863826561771683774617978763837906638e-7"]

    result = zeros(T, 11)

    for j = 1:6
        i = div(j + 1, 2)
        result[j] = if j % 2 == 1
            dt * (b[i] + 2c[i] * dt^2 + 4d[i] * dt^4 - 8e[i] * dt^6)
        else
            dt * a[i]
        end
        result[12-j] = result[j]
    end

    return result
end

function table2split(dt::T) where {T}
    a = [big"1.079852426382430882456991", -big"0.579852426382430882456991", 0]
    b = [
        big"0.359950808794143627485664",
        -big"0.1437147273026540434771131",
        big"0.567527837017020831982899",
    ]
    c = [0, -big"0.0139652542242388403673", -big"0.039247029382345626020"]
    note(a, b, c, true)
    result = zeros(T, 9)
    for j = 1:5
        i = div(j + 1, 2)
        result[j] = if j % 2 == 1
            dt * (b[i] + 2c[i] * dt^2)
        else
            dt * a[i]
        end
        result[10-j] = result[j]
    end
    return result
end

"""
    Advection{T}
    Advection(
        t_mesh::NTuple{N,UniformMesh{T}},
        t_interp::Vector{I}
        dt_base::T;
        tab_coef = [1 // 1],
        tab_fct = missing,
        timeopt::TimeOptimization = NoTimeOpt,
    ) where {T, N, I <: AbstractInterpolation{T}}

Immutable structure that contains constant parameters for multidimensional advection

# Type parameters

- `T::DataType` : type of data
- `N` : number of dimensions
- `I` : commun type of interpolation
- `timeopt::TimeOptimization` : time optimization

# Arguments

- `t_mesh::NTuple{N, UniformMesh{T}}` : tuple of meshes (one for each dimension)
- `t_interp::Vector{I}` : tuple of interpolations (one for each dimension)
- `dt_base::T` : time delta for one advection series

# Keywords

- `tab_coef=[1//2, 1//1, 1//2]` : coefficient table for one advection series, the
    coefficients at odd indexes is for space advection series, the coefficients at even indexes is for velocity advection series
- `tab_fct=[identity, identity, identity]` : function table for one advection series, with the same indexes than tab_coef

# Implementation
- `sizeall` : tuple of the sizes of all dimensions (space before velocity)
- `t_mesh_sp` : tuple of space meshes
- `t_mesh_v` : tuple of velocity meshes
- `t_interp_sp` : tuple of space interpolation types
- `t_interp_v` : tuple of velocity interpolation types
- `dt_base::T` : time unit of an advection series
- `tab_coef` : coefficient table
- `v_square` : precompute for ke
- `nbsplit` : number of slices for split
- `mpiid` : MPI id

# Throws

- `ArgumentError` : `Nsp` must be less or equal to `Nv`.

"""
struct Advection{T,N,I,timeopt,timealg,ordalg}
    sizeall::NTuple{N,Int}
    t_mesh::NTuple{N,UniformMesh{T}}
    t_interp::Vector{I}
    dt_base::T
    states::Vector{StateAdv{N}}
    maxcoef::Int
    nbstates::Int
    tab_coef::Vector{T}
    nbsplit::Int
    mpid::Any
    abcoef::ABcoef
    tabmod::NTuple{N,Vector{Int}}
    function Advection(
        t_mesh::NTuple{N,UniformMesh{T}},
        t_interp::Vector{I},
        dt_base::T,
        states::Vector{Tuple{Vector{Int},Int,Int,Bool,Vararg{Bool,N2}}};
        tab_coef::Vector{T} = strangsplit(dt_base),
        timeopt::TimeOptimization = NoTimeOpt,
        timealg::TimeAlgorithm = NoTimeAlg,
        ordalg::Int = timealg != NoTimeAlg ? 4 : 0,
    ) where {T,N,N2,I<:AbstractInterpolation{T}}
        length(t_interp) == N ||
            throw(ArgumentError("size of vector of Interpolation must be equal to N=$N"))
        sizeall = length.(t_mesh)

        newstates = map(i -> StateAdv(i, states[i]...), 1:length(states))

        maxcoef = maximum(x -> x.stcoef, newstates)

        restcoef = length(tab_coef) % maxcoef

        nbstatesplus = length(filter(x -> x.stcoef in restcoef, newstates))

        nbstates = div(length(tab_coef), maxcoef) * length(states) + nbstatesplus

        mpid = timeopt == MPIOpt ? MPIData() : missing
        nbsplit = if timeopt == MPIOpt
            mpid.nb
        elseif timeopt == SplitThreadsOpt
            Threads.nthreads()
        else
            1
        end
        return new{T,N,I,timeopt,timealg,ordalg}(
            sizeall,
            t_mesh,
            t_interp,
            dt_base,
            newstates,
            maxcoef,
            nbstates,
            tab_coef,
            nbsplit,
            mpid,
            ABcoef(ordalg + 1),
            gettabmod.(sizeall),
        )
    end
end
"""
    sizeall(adv::Advection)

Return a tuple of the sizes of each dimensions

# Argument
- `adv::Advection` : Advection structure.
"""
sizeall(adv::Advection) = adv.sizeall

getst(adv::Advection, x) = adv.states[modone(x, length(adv.states))]

function getstcoef(adv::Advection, x)
    return div(x - 1, length(adv.states)) * adv.maxcoef + getst(adv, x).stcoef
end

getcur_t(adv::Advection, x) = adv.tab_coef[getstcoef(adv, x)]

function getinterp(adv::Advection, x)
    st = getst(adv, x)
    return adv.t_interp[st.perm[1:(st.ndims)]]
end

function getordalg(
    adv::Advection{T,N,I,timeopt,timealg,ordalg},
) where {T,N,I,timeopt,timealg,ordalg}
    return ordalg
end

# Interface of external data
abstract type AbstractExtDataAdv end

function initfmrdata(adv::Advection, bufdata::Vector{T}, state) where {T}
    p = getst(adv, state).perm
    if isbitstype(T)
        ptr = pointer(bufdata)
        # the following call is safe because :
        #       - T type is an bits type
        #  and  - length(bufdata) == prod(sizeall(adv))
        return unsafe_wrap(Array, ptr, sizeall(adv)[p]; own = false)
    else
        return zeros(T, sizeall(adv)[p])
    end
end

"""
    AdvectionData{T,N,timeopt}
    AdvectionData(
    adv::Advection{T,N,timeopt},
    data::Array{T,N},
    parext)

Mutable structure that contains variable parameters of advection series

# Type parameters
- `T::DataType` : type of data
- `N` : number of dimensions
- `timeopt::TimeOptimization` : time optimization

# Arguments
- `adv::Advection{T,N}` : link to the constant data of this advection
- `data::Array{T,Nsum}` : Initial data of this advection
- `parext` : external data of this advection to compute alpha of each interpolations

# Implementation
- `adv::Advection{T,N,timeopt}` : link to the constant data of this advection
- `state_coef::Int` : state that is the index of `tab_coef`, it is from one to lenth(tab_coef)
- `state_dim::Int` : the dimension index, from 1 to Nsp in space states, from one to Nv in velocity state
- `data::Array{T,Nsum}` : it is the working buffer
- `bufdata::Vector{T}` : vector of the same size of the working buffer
- `fmrtabdata::NTuple{Nsum,Array{T,Nsum}}` : tuple of array with the same size than data but with permutated dimensions
- `t_buf::NTuple{Nsum, Array{T,2}}` : tuple of buffer that is used to get the linear data for interpolation, one buffer per thread
- `cache_alpha::Union{T,Nothing}` : cache for precal, the precal is compute only when the alpha or decint values change
- `cache_decint::Int64` : for precal cache
- `cache_precal::Vector{T}` : for precal cache
- `parext::ExtDataAdv` : external data of this advection to compute alpha of each interpolations

# Methods to define
- `initcoef!(parext::AbstractExtDataAdv, self::Advection1dData)` : this method called at the beginning of each advection to initialize parext data. The `self.parext` mutable structure is the only data that initcoef! can modify otherwise it leads to unpredictable behaviour.
- `getalpha(parext::AbstractExtDataAdv, self::Advection1dData, ind)` : return the alpha number that is used for interpolation.
- `getperm(parext::AbstractExtDataAdv, advd::Advection1dData)` : get the permutation of the dimension as a function of the current state, the dimension where advection occurs must be first, the dimensions used to compute alpha must be at the end.

"""
mutable struct AdvectionData{T,N,timeopt,timealg}
    adv::Advection{T,N}
    state_gen::Int # indice of the calls of advection!
    time_cur::T # current time
    data::Array{T,N}
    bufdata::Array{T}
    fmrtabdata::Vector{Array{T,N}}
    t_buf::Vector{Array{T}}
    t_itr::Any
    tt_split::Any
    t_cache::Vector{Vector{CachePrecal{T}}}
    parext::AbstractExtDataAdv
    bufcur::Union{Array{OpTuple{N,T},N},Missing}
    t_bufc::Vector{Array{OpTuple{N,T},N}}
    initdatas::Union{Vector{Array{T,N}},Missing}
    function AdvectionData(
        adv::Advection{T,N,I,timeopt,timealg},
        data::Array{T,N},
        parext::AbstractExtDataAdv;
        initdatas::Union{Vector{Array{T,N}},Missing} = missing,
        time_init::T = zero(T),
    ) where {T,N,I,timeopt,timealg}
        s = size(data)
        s == sizeall(adv) ||
            thrown(ArgumentError("size(data)=$s it must be $(sizeall(adv))"))
        nbst = length(adv.states)
        nbthr = if timeopt == SimpleThreadsOpt || timeopt == SplitThreadsOpt
            Threads.nthreads()
        else
            1
        end
        t_buf = map(
            x -> zeros(T, s[getst(adv, x).perm][1:(getst(adv, x).ndims)]..., nbthr),
            1:nbst,
        )
        datanew = Array{T,N}(undef, s)
        bufdata = Vector{T}(undef, length(data))
        fmrtabdata = map(x -> initfmrdata(adv, bufdata, x), 1:nbst)
        copyto!(datanew, data)
        #        @show adv.nbsplit
        if nbst == 1
            t_itr = (splitvec(adv.nbsplit, CartesianIndices(s)),)
        else
            t_itr = ntuple(
                x -> splitvec(
                    adv.nbsplit,
                    CartesianIndices(s[adv.states[x].perm][(adv.states[x].ndims+1):N]),
                ),
                nbst,
            )
        end
        #        @show t_itr
        t_linind = ntuple(x -> LinearIndices(s[adv.states[x].perm]), nbst)
        cartz(x) = CartesianIndices(s[adv.states[x].perm][1:(adv.states[x].ndims)])
        fbegin(x, y) = t_linind[x][cartz(x)[1], t_itr[x][y][1]]
        fend(x, y) = t_linind[x][cartz(x)[end], t_itr[x][y][end]]
        if nbst == 1
            li = LinearIndices(s)
            it = t_itr[1]
            tt_split = (ntuple(y -> (li[it[y][1]]:li[it[y][end]]), adv.nbsplit),)
        else
            tt_split =
                ntuple(x -> ntuple(y -> (fbegin(x, y):fend(x, y)), adv.nbsplit), nbst)
        end
        t_cache =
            map(x -> map(i -> CachePrecal(getinterp(adv, x), zero(T)), 1:nbthr), 1:nbst)

        return new{T,N,timeopt,timealg}(
            adv,
            1,
            time_init,
            datanew,
            bufdata,
            fmrtabdata,
            t_buf,
            t_itr,
            tt_split,
            t_cache,
            parext,
            missing,
            [],
            initdatas,
        )
    end
end

getst(self::AdvectionData) = getst(self.adv, self.state_gen)
getext(self) = self.parext
getdata(self) = self.data
getnbdims(self::AdvectionData) = getst(self).ndims
getstcoef(self::AdvectionData) = getstcoef(self.adv, self.state_gen)
function getcur_t(self::AdvectionData, extdata::AbstractExtDataAdv)
    return getcur_t(self.adv, self.state_gen)
end
getcur_t(self::AdvectionData) = getcur_t(self, self.parext)
_getcurrentindice(self::AdvectionData) = getst(self).perm[1]
function getindsplit(self::AdvectionData{T,N,timeopt}) where {T,N,timeopt}
    if self.adv.nbsplit != 1
        ind = timeopt == MPIOpt ? self.adv.mpid.ind : Threads.threadid()
    else
        ind = 1
    end
    return ind
end
getinterp(self::AdvectionData) = getinterp(self.adv, self.state_gen)

getitr(self::AdvectionData) = self.t_itr[getst(self).ind][getindsplit(self)]
gett_split(self::AdvectionData) = self.tt_split[getst(self).ind]

"""
    nextstate!(self::AdvectionData{T, N})

Function called at the end of advection function to update internal state of AdvectionData structure

# Argument
- `self::AdvectionData{T, N}` : object to update

# return value
- `ret::Bool` : `true` if the series must continue
                `false` at the end of the series.
"""
function retns(self::AdvectionData, extdata::AbstractExtDataAdv)
    return false
end
function nextstate!(self::AdvectionData)
    if self.state_gen < self.adv.nbstates
        self.state_gen += 1
        return true
    else
        self.state_gen = 1
        self.time_cur += self.adv.dt_base
        return retns(self, self.parext)
    end
end

function getformdata(advd::AdvectionData)
    f = advd.fmrtabdata[getst(advd).ind]
    permutedims!(f, advd.data, getst(advd).perm)
    return f
end

function copydata!(advd::AdvectionData{T,N,timeopt,timealg}, f) where {T,N,timeopt,timealg}
    if timeopt == MPIOpt && advd.adv.nbsplit != 1 && length(advd.adv.states) != 1
        mpibroadcast(advd.adv.mpid, gett_split(advd), f)
    end
    return permutedims!(advd.data, f, invperm(getst(advd).perm))
end

function decbegin!(t_trv, t_cal, t_interp::Vector{I}) where {I<:AbstractInterpolation}
    indice = length(t_trv)
    for i = 1:(indice-1)
        buf = t_cal[end]
        autointerp!(buf, copy(buf), indice - 1, t_interp)
        interpbufc!(t_trv, buf, t_interp, i)
        deleteat!(t_cal, length(t_cal))
    end
end

function initcoef!(self::AdvectionData{T,N,timeopt,timealg}) where {T,N,timeopt,timealg}
    nbtours = 3
    isbegin = ismissing(self.bufcur)
    extdata::AbstractExtDataAdv = getext(self)
    initcoef!(extdata, self)
    cachethreads =
        timeopt in (SimpleThreadsOpt, SplitThreadsOpt) ? self.t_cache[1] : missing
    t_sp = timeopt in (MPIOpt, SplitThreadsOpt) ? self.tt_split[1] : missing

    if timealg == ABTimeAlg_new
        adv = self.adv
        ordalg = getordalg(adv)
        if isbegin
            t_ref = []
            t_cal = []
            push!(t_ref, copy(self.bufcur))
            svdata = copy(self.data)
            svbufcur = copy(self.bufcur)
            sens = ordalg * nbtours % 2 == 1 ? 1 : -1
            f = zeros(T, sizeall(adv))
            for indice = 1:ordalg, nb = 1:(indice == ordalg ? nbtours - 1 : nbtours)
                t_ref = reverse(t_ref)
                t_trv = sens * copy.(t_ref)

                decbegin!(t_trv, t_cal, adv.t_interp)

                t_cal = []
                copy!(self.data, svdata)
                for i = 1:indice
                    fmrdec = sum(map(k -> c(adv.abcoef, k, indice) * t_trv[k], 1:indice))
                    if i != 1
                        push!(t_cal, -fmrdec)
                    end
                    if i != 1 || nb != 1
                        deleteat!(t_trv, length(t_trv))
                        deleteat!(t_ref, length(t_ref))
                    end
                    autointerp!(
                        self.bufcur,
                        fmrdec,
                        indice,
                        adv.t_interp;
                        mpid = adv.mpid,
                        t_split = t_sp,
                        cachethreads = cachethreads,
                    )
                    interpbufc!(
                        self.t_bufc,
                        self.bufcur,
                        adv.t_interp;
                        mpid = adv.mpid,
                        t_split = t_sp,
                        cachethreads = cachethreads,
                    )
                    interpolate!(
                        f,
                        self.data,
                        self.bufcur,
                        adv.t_interp;
                        mpid = adv.mpid,
                        t_split = t_sp,
                        cachethreads = cachethreads,
                    )
                    copy!(self.data, f)
                    initcoef!(extdata, self) # calculate bufcur
                    pushfirst!(t_trv, copy(self.bufcur))
                    pushfirst!(t_ref, copy(self.bufcur))
                end

                fmrdec =
                    -sum(map(i -> c(adv.abcoef, i, indice + 1) * t_trv[i], 1:(indice+1)))
                push!(t_cal, fmrdec)

                sens = -sens
            end
            @assert sens == 1 "sens must be positive at this place"

            t_ref = reverse(t_ref)
            t_trv = sens * copy.(t_ref)
            deleteat!(t_trv, 1)
            decbegin!(t_trv, t_cal, adv.t_interp)
            self.t_bufc = t_trv
            self.data .= svdata
            self.bufcur .= svbufcur
        end
    end

    # old version that works at order 4 for poisson and order 2 for other
    if timealg == ABTimeAlg_ip
        adv = self.adv
        ordalg = getordalg(adv)
        if isbegin
            for indice = 1:(ordalg-1)
                pushfirst!(self.t_bufc, copy(self.bufcur))
                fmrdec = sum(map(i -> c(adv.abcoef, i, indice) * self.t_bufc[i], 1:indice))
                autointerp!(
                    fmrdec,
                    copy(fmrdec),
                    indice - 1,
                    adv.t_interp;
                    mpid = adv.mpid,
                    t_split = t_sp,
                    cachethreads = cachethreads,
                )
                interpbufc!(
                    self.t_bufc,
                    fmrdec,
                    adv.t_interp;
                    mpid = adv.mpid,
                    t_split = t_sp,
                    cachethreads = cachethreads,
                )
            end
        end
    end
    if timealg == ABTimeAlg_init
        adv = self.adv
        ordalg = getordalg(adv)
        if isbegin
            for indice = 1:length(self.initdatas)
                pushfirst!(self.t_bufc, copy(self.bufcur))
                ord = min(indice, ordalg)
                fmrdec = sum(map(i -> c(adv.abcoef, i, ord) * self.t_bufc[i], 1:ord))
                if ord == ordalg
                    deleteat!(self.t_bufc, length(self.t_bufc))
                end
                autointerp!(
                    fmrdec,
                    copy(fmrdec),
                    ordalg - 1,
                    adv.t_interp;
                    mpid = adv.mpid,
                    t_split = t_sp,
                    cachethreads = cachethreads,
                )
                interpbufc!(
                    self.t_bufc,
                    fmrdec,
                    adv.t_interp;
                    mpid = adv.mpid,
                    t_split = t_sp,
                    cachethreads = cachethreads,
                )
                copy!(self.data, self.initdatas[indice])
                self.time_cur += getcur_t(self)
                initcoef!(extdata, self)
            end
        end
    end

    if timealg in (ABTimeAlg_ip, ABTimeAlg_new, ABTimeAlg_init)
        pushfirst!(self.t_bufc, copy(self.bufcur))

        bufc = sum(map(i -> c(adv.abcoef, i, ordalg) * self.t_bufc[i], 1:ordalg))

        autointerp!(
            self.bufcur,
            bufc,
            ordalg - 1,
            adv.t_interp;
            mpid = adv.mpid,
            t_split = t_sp,
            cachethreads = cachethreads,
        )

        deleteat!(self.t_bufc, length(self.t_bufc))

        interpbufc!(
            self.t_bufc,
            self.bufcur,
            adv.t_interp;
            mpid = adv.mpid,
            t_split = t_sp,
            cachethreads = cachethreads,
        )
    end
end

"""
    advection!(self::AdvectionData)

Advection function of a multidimensional function `f` discretized on `mesh`

# Argument
- `self::AdvectionData` : mutable structure of variables data

# Return value
- `true` : means that the advection series must continue
- `false` : means that the advection series is ended.
"""
function advection!(self::AdvectionData{T,N,timeopt,timealg}) where {T,N,timeopt,timealg}
    fltrace = true
    adv = self.adv
    interp = getinterp(self)
    extdata = getext(self)
    initcoef!(self)
    curind = _getcurrentindice(self)
    f = getformdata(self)
    st = getst(self)
    sz = size(f)[1:(st.ndims)]
    tabmod = adv.tabmod[st.perm]  # just for optimization of interpolation!

    coltuple = ntuple(x -> Colon(), st.ndims)

    if length(self.adv.states) == 1
        isthreads = timeopt in (SimpleThreadsOpt, SplitThreadsOpt)
        t_sp = timeopt in (MPIOpt, SplitThreadsOpt) ? self.tt_split[1] : missing
        cachethreads = isthreads ? self.t_cache[st.ind] : missing
        interpolate!(
            f,
            self.data,
            self.bufcur,
            adv.t_interp;
            tabmod = tabmod,
            mpid = adv.mpid,
            t_split = t_sp,
            cachethreads = cachethreads,
        )
    elseif timeopt == NoTimeOpt || timeopt == MPIOpt
        local buf = view(self.t_buf[st.ind], coltuple..., 1)
        local cache = self.t_cache[st.ind][1]
        local itr = getitr(self)
        if st.isconstdec
            for indext in itr
                local decint, precal = getprecal(cache, getalpha(extdata, self, indext))
                local slc = view(f, coltuple..., indext)
                interpolate!(buf, slc, decint, precal, interp, tabmod)
                slc .= buf
            end
        else
            for indext in itr
                local slc = view(f, coltuple..., indext)
                interpolate!(
                    buf,
                    slc,
                    indbuf -> getalpha(extdata, self, indext, indbuf),
                    interp;
                    tabmod = tabmod,
                    cache = cache,
                )
                slc .= buf
            end
        end
    elseif timeopt == SimpleThreadsOpt
        local itr = collect(getitr(self))
        if st.isconstdec
            @threads for indext in itr
                local buf = view(self.t_buf[st.ind], coltuple..., Threads.threadid())
                local cache = self.t_cache[st.ind][Threads.threadid()]
                local decint, precal = getprecal(cache, getalpha(extdata, self, indext))
                local slc = view(f, coltuple..., indext)
                interpolate!(buf, slc, decint, precal, interp, tabmod)
                slc .= buf
            end
        else
            @threads for indext in itr
                local buf = view(self.t_buf[st.ind], coltuple..., Threads.threadid())
                local cache = self.t_cache[st.ind][Threads.threadid()]
                local slc = view(f, coltuple..., indext)
                interpolate!(
                    buf,
                    slc,
                    indbuf -> getalpha(extdata, self, indext, indbuf),
                    interp,
                    tabmod,
                    cache,
                )
                slc .= buf
            end
        end
    elseif timeopt == SplitThreadsOpt
        Threads.@threads for indth = 1:Threads.nthreads()
            local buf = view(self.t_buf[st.ind], coltuple..., Threads.threadid())
            local cache = self.t_cache[st.ind][Threads.threadid()]
            local itr = getitr(self)
            if st.isconstdec
                for indext in itr
                    local decint, precal = getprecal(cache, getalpha(extdata, self, indext))
                    local slc = view(f, coltuple..., indext)
                    interpolate!(buf, slc, decint, precal, interp, tabmod)
                    slc .= buf
                end
            else
                for indext in itr
                    local slc = view(f, coltuple..., indext)
                    interpolate!(
                        buf,
                        slc,
                        indbuf -> getalpha(extdata, self, indext, indbuf),
                        interp,
                        tabmod,
                        cache,
                    )
                    slc .= buf
                end
            end
        end
    end
    copydata!(self, f)
    return nextstate!(self)
end
