
@enum TimeOptimization NoTimeOpt = 1 SimpleThreadsOpt = 2 SplitThreadsOpt = 3 MPIOpt = 4

struct StateAdv{N}
    ind::Int          # indice
    perm::Vector{Int} # dimensions permutation
    invp::Vector{Int} # inverse of dimension permutation
    ndims::Int        # count of dimension
    stcoef::Int   # state_coef
    isconstdec::Bool #true if constant dec
    StateAdv(ind, p, ndims, stc, isconst)=new{length(p)}(ind, p, invperm(p), ndims, stc, isconst)
end


standardsplit(T::DataType)=T.([1,1])
standardsplit()=standardsplit(Float64)
strangsplit(T::DataType)=T.([1//2,1//1,1//2])
strangsplit()=strangsplit(Float64)
function triplejumpsplit(T::DataType)
    c=T(2)^(1//3)
    c1 = 1/(2(2-c))
    c2 = (1-c)/(2(2-c))
    d1 = 1/(2-c)
    d2 = -c/(2-c)
    return [c1, d1, c2, d2, c2, d1, c1]
end
triplejumpsplit()=triplejumpsplit(Float64)

order6split()=[ 
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
    0.0414649985182624
]
function order6split(T::DataType)
    b = convert(Vector{T},copy(order6split()))
    b[11] = b[13] = T(1//2) - sum(b[1:2:9])
    b[12] = T(1//1) - 2sum(b[2:2:10])
    @assert isapprox(sum(b),T(2))
    return b
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
- `tab_fct` : function table
- `v_square` : precompute for ke
- `nbsplit` : number of slices for split
- `mpiid` : MPI id

# Throws

- `ArgumentError` : `Nsp` must be less or equal to `Nv`.

"""
struct Advection{T, N, I, timeopt}
    sizeall::NTuple{N, Int}
    t_mesh::NTuple{N,UniformMesh{T}}
    t_interp::Vector{I}
    dt_base::T
    states::Vector{StateAdv{N}}
    maxcoef::Int
    nbstates::Int
    tab_coef::Vector{T}
    tab_fct::Vector{Function}
    nbsplit::Int
    mpid::Any
    function Advection(
        t_mesh::NTuple{N,UniformMesh{T}},
        t_interp::Vector{I},
        dt_base::T,
        states::Vector{Tuple{Vector{Int}, Int, Int, Bool}};
        tab_coef = [1 // 2, 1 // 1, 1 // 2],
        tab_fct = missing,
        timeopt::TimeOptimization = NoTimeOpt,
    ) where {T, N, I <: AbstractInterpolation{T}}
        length(t_interp) == N || throw(ArgumentError("size of vector of Interpolation must be equal to N=$N"))
        sizeall = length.(t_mesh)

        newstates=map( i -> StateAdv(i, states[i]...), 1:length(states))

        maxcoef = maximum(x->x.stcoef, newstates)

        restcoef = length(tab_coef) % maxcoef

        nbstatesplus = length( filter( x -> x.stcoef in restcoef, newstates))

        nbstates = div(length(tab_coef), maxcoef)*length(states) + nbstatesplus

#        v_square = dotprod(points.(t_mesh_v)) .^ 2 # precompute for ke
        mpid = timeopt == MPIOpt ? MPIData() : missing
        nbsplit = if timeopt == MPIOpt
            mpid.nb
        elseif timeopt == SplitThreadsOpt
            Threads.nthreads()
        else
            1
        end
        nfct = ismissing(tab_fct) ? fill(identity, size(tab_coef)) : tab_fct
        return new{T, N, I, timeopt}(
            sizeall,
            t_mesh,
            t_interp,
            dt_base,
            newstates,
            maxcoef,
            nbstates,
            tab_coef,
            nfct,
#            v_square,
            nbsplit,
            mpid,
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

getst(adv::Advection, x) = adv.states[modone(x,length(adv.states))]

getstcoef(adv::Advection, x) = div(x-1, length(adv.states))*adv.maxcoef+ getst(adv,x).stcoef

function getcur_t(adv::Advection, x)
    state_coef = getstcoef(adv, x)
    return adv.tab_fct[state_coef](adv.tab_coef[state_coef] * adv.dt_base)
end

function getinterp(adv::Advection, x)
    st = getst(adv, x)
    return adv.t_interp[st.perm[1:st.ndims]]
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
        return unsafe_wrap(Array, ptr, sizeall(adv)[p], own=false)
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
mutable struct AdvectionData{T,N,timeopt}
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
    clobs::AbstractClockObs
    function AdvectionData(
        adv::Advection{T,N,I, timeopt},
        data::Array{T,N},
        parext::AbstractExtDataAdv;
        clockobs=false
    ) where {T, N, I, timeopt}
        s = size(data)
        s == sizeall(adv) ||
            thrown(ArgumentError("size(data)=$s it must be $(sizeall(adv))"))
        nbst=length(adv.states)
        nbthr =
            timeopt == SimpleThreadsOpt || timeopt == SplitThreadsOpt ? Threads.nthreads() :
            1
        t_buf = map(x -> zeros(T, s[getst(adv,x).perm][1:getst(adv,x).ndims]...,nbthr), 1:nbst)
        datanew = Array{T,N}(undef, s)
        bufdata = Vector{T}(undef, length(data))
        fmrtabdata = map(x -> initfmrdata(adv, bufdata, x), 1:nbst)
        copyto!(datanew, data)
#        @show adv.nbsplit
        t_itr = ntuple(
            x -> splitvec(adv.nbsplit, CartesianIndices(s[adv.states[x].perm][(adv.states[x].ndims + 1):N])),
            nbst,
        )
#        @show t_itr
        t_linind = ntuple(x -> LinearIndices(s[adv.states[x].perm]), nbst)
        cartz(x)=CartesianIndices(s[adv.states[x].perm][1:adv.states[x].ndims])
        fbegin(x, y) = t_linind[x][cartz(x)[1], t_itr[x][y][1]]
        fend(x, y) = t_linind[x][cartz(x)[end], t_itr[x][y][end]]
        tt_split = ntuple(x -> ntuple(y -> (fbegin(x, y):fend(x, y)), adv.nbsplit), nbst)
        t_cache = map(x->map( i -> CachePrecal(getinterp(adv,x), zero(T)), 1:nbthr), 1:nbst)
        return new{T,N,timeopt}(
            adv,
            1,
            zero(T),
            datanew,
            bufdata,
            fmrtabdata,
            t_buf,
            #    t_itrfirst, t_itrsecond, tt_split,
            t_itr,
            tt_split,
            t_cache,
            parext,
            clockobs ? ClockObs(10) : NoClockObs(),
        )
    end
end

getst(self::AdvectionData) = getst(self.adv,self.state_gen)
getext(self) = self.parext
getdata(self) = self.data
getnbdims(self::AdvectionData)=getst(self).ndims
getstcoef(self::AdvectionData)=getstcoef(self.adv, self.state_gen)
getcur_t(self::AdvectionData) = getcur_t(self.adv, self.state_gen)
isvelocitystate(state_coef::Int) = state_coef % 2 == 0
isvelocity(adv::Advection, curid) = isvelocitystate(getstcoef(adv,curid))
isvelocitystate(self::AdvectionData) = isvelocitystate(getstcoef(self))
_getcurrentindice(self::AdvectionData)=getst(self).perm[1]
function getindsplit(self::AdvectionData{T,N,timeopt}) where {T,N,timeopt}
    if self.adv.nbsplit != 1
        ind = timeopt == MPIOpt ? self.adv.mpid.ind : Threads.threadid()
    else
        ind = 1
    end
    return ind
end
getinterp(self::AdvectionData)=getinterp(self.adv, self.state_gen)

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
function nextstate!(self::AdvectionData)
    if self.state_gen < self.adv.nbstates
        self.state_gen += 1
        return true
    else
        self.state_gen = 1
        printall(self.clobs)
        self.time_cur += self.adv.dt_base
        return false
    end
end
# default function of the interface
# initcoef!(parext::AbstractExtDataAdv, self::AdvectionData) = missing

# # this interface function must always be defined
# function getalpha(parext::AbstractExtDataAdv, self::AdvectionData, indext)
#     throw(error("getalpha undefined for $(typeof(parext))"))
# end
# # If not defined we try with only external index
# function getalpha(parext::AbstractExtDataAdv, self::AdvectionData, indext::CartesianIndex, ind::CartesianIndex)
#     return getalpha(parext, self, indext)
# end


# data formating
function getformdata(advd::AdvectionData)
    # ptr = pointer(advd.bufdata)
    # f = unsafe_wrap(Array, ptr, sizeall(advd.adv)[p], own=false)
    f = advd.fmrtabdata[getst(advd).ind]
    permutedims!(f, advd.data, getst(advd).perm)
    return f
end
function copydata!(
    advd::AdvectionData{T,N,timeopt},
    f,
) where {T,N,timeopt}
    if timeopt == MPIOpt && advd.adv.nbsplit != 1
        mpibroadcast(advd.adv.mpid, gett_split(advd), f)
    end
    permutedims!(advd.data, f, invperm(getst(advd).perm))
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
function advection!(
    self::AdvectionData{T,N,timeopt},
) where {T,N,timeopt}
    fltrace = true
    f = self.data
    interp = getinterp(self)
    extdata = getext(self)
    initcoef!(extdata, self)
    curind = _getcurrentindice(self)
    f = getformdata(self)
    st=getst(self)
    sz = size(f)[1:st.ndims]
#    @show size(f)
    tabmod = gettabmod.(sz) # just for optimization of interpolation!
#    @show sz, timeopt

    coltuple = ntuple( x -> Colon(), st.ndims)

    if timeopt == NoTimeOpt || timeopt == MPIOpt
        #        @inbounds for ind in getitr(self)
        # fmrcpt=1
        local buf = view(self.t_buf[st.ind], coltuple..., 1)
        local cache = self.t_cache[st.ind][1]
        local itr = getitr(self)
        # local colitr = collect(itr)
        # @show length(colitr), colitr[1]
#       @show itr
# clockbegin(self.clobs,1)
        if st.isconstdec
            @inbounds for indext in itr
                local decint, precal = getprecal(cache, getalpha(extdata, self, indext))
                local slc = view(f, coltuple..., indext)
                # if fltrace
                #     @show typeof(buf),typeof(slc),length(buf),length(slc)
                #     fltrace = false
                # end
#                clockbegin(self.clobs,2)
#              interpolate!(buf, slc, decint, precal, interp, tabmod; clockobs=self.clobs)
                interpolate!(buf, slc, decint, precal, interp, tabmod)
                slc .= buf
#                clockend(self.clobs,2)
            end           
        else
            for indext in itr
                local slc = view(f, coltuple..., indext)
    #           @show ind
                interpolate!(buf, slc, indbuf -> getalpha(extdata, self, indext, indbuf), interp; tabmod=tabmod, cache=cache)
                slc .= buf
            end
        end
# clockend(self.clobs,1)
    elseif timeopt == SimpleThreadsOpt
        #        @inbounds begin
        local itr = collect(getitr(self))
        if st.isconstdec
#            @show itr
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
                local buf = 
                local cache = self.t_cache[st.ind][Threads.threadid()]
                local slc = view(f, coltuple..., indext)
                interpolate!(buf, slc, indbuf -> getalpha(extdata, self, indext, indbuf), interp, tabmod, cache)
                slc .= buf
            end
        end
        #        end
    elseif timeopt == SplitThreadsOpt
        #        @inbounds begin
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
        #           @show ind
                    interpolate!(buf, slc, indbuf -> getalpha(extdata, self, indext, indbuf), interp, tabmod, cache)
                    slc .= buf
                end
            end
        end
        #        end
    end
    copydata!(self, f)
    return nextstate!(self)
end
