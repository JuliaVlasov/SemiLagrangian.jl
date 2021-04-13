
@enum TimeOptimization NoTimeOpt = 1 SimpleThreadsOpt = 2 SplitThreadsOpt = 3 MPIOpt = 4

struct StateAdv{N}
    perm::Vector{Int} # dimensions permutation
    invp::Vector{Int} # inverse of dimension permutation
    ndims::Int        # count of dimension
    stcoef::Int   # state_coef
    isconstdec::Bool #true if constant dec
    StateAdv(p, ndims, stc, isconst)=new{length(p)}(p, invperm(p), ndims, stc, isconst)
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
        indices=1
        newstates=map( i -> StateAdv(states[i]...), 1:length(states))

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

getst(adv::Advection, x) = adv.states[x]

# Interface of external data
abstract type AbstractExtDataAdv{T,N} end


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
    data::Array{T,N}
    bufdata::Array{T}
    fmrtabdata::Vector{Array{T,N}}
    t_buf::Vector{Vector{Array{T}}}
    t_itr::Any
    tt_split::Any
    t_cache::Vector{Vector{CachePrecal{T}}}
    parext::AbstractExtDataAdv
    function AdvectionData(
        adv::Advection{T,N,I, timeopt},
        data::Array{T,N},
        parext::AbstractExtDataAdv,
    ) where {T, N, I, timeopt}
        s = size(data)
        s == sizeall(adv) ||
            thrown(ArgumentError("size(data)=$s it must be $(sizeall(adv))"))
        nbst=size(adv.states,1)
        nbthr =
            timeopt == SimpleThreadsOpt || timeopt == SplitThreadsOpt ? Threads.nthreads() :
            1
        t_buf = map( x -> map(y -> zeros(T, s[getst(adv,x).perm][1:getst(adv,x).ndims]), 1:nbthr), 1:nbst)
        datanew = Array{T,N}(undef, s)
        bufdata = Vector{T}(undef, length(data))
        fmrtabdata = map(x -> zeros(T, s[adv.states[x].perm]), 1:nbst)
        copyto!(datanew, data)
        @show adv.nbsplit
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
            datanew,
            bufdata,
            fmrtabdata,
            t_buf,
            #    t_itrfirst, t_itrsecond, tt_split,
            t_itr,
            tt_split,
            t_cache,
            parext,
        )
    end
end


getext(self) = self.parext
getdata(self) = self.data
getnbdims(self::AdvectionData)=self.adv.states[self.state_gen].ndims
getstcoef(self::AdvectionData)=self.adv.states[self.state_gen].stcoef
getcur_t(adv::Advection, state_coef::Int) =
    adv.tab_fct[state_coef](adv.tab_coef[state_coef] * adv.dt_base)
getcur_t(self::AdvectionData) = getcur_t(self.adv, getst(self).stcoef)
isvelocitystate(state_coef::Int) = state_coef % 2 == 0
isvelocity(adv::Advection, curid) = isvelocitystate(adv.states[curid].stcoef)
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

getst(self::AdvectionData)=self.adv.states[self.state_gen]
function getinterp(self::AdvectionData)
    st = getst(self)
    return self.adv.t_interp[st.perm[1:st.ndims]]
end
function getinterp(adv::Advection, x)
    st = adv.states[x]
    return adv.t_interp[st.perm[1:st.ndims]]
end


getitr(self::AdvectionData) = self.t_itr[self.state_gen][getindsplit(self)]
gett_split(self::AdvectionData) = self.tt_split[self.state_gen]



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
    if self.state_gen < length(self.adv.states)
        self.state_gen += 1
        return true
    else
        self.state_gen = 1
        return false
    end
end
# default function of the interface
initcoef!(parext::AbstractExtDataAdv, self::AdvectionData) = missing

# this interface function must always be defined
function getalpha(parext::AbstractExtDataAdv, self::AdvectionData, indext)
    throw(error("getalpha undefined for $(typeof(parext))"))
end
# If not defined we try with only external index
function getalpha(parext::AbstractExtDataAdv, self::AdvectionData, indext, ind)
    return getalpha(parext, self, indext)
end

# data formating
function getformdata(advd::AdvectionData)
    # ptr = pointer(advd.bufdata)
    # f = unsafe_wrap(Array, ptr, sizeall(advd.adv)[p], own=false)
    f = advd.fmrtabdata[advd.state_gen]
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
        local buf = self.t_buf[self.state_gen][1]
        local cache = self.t_cache[self.state_gen][1]
        local itr = getitr(self)

#       @show itr
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
    elseif timeopt == SimpleThreadsOpt
        #        @inbounds begin
        local itr = collect(getitr(self))
        if st.isconstdec
#            @show itr
            @threads for indext in itr
                local buf = self.t_buf[self.state_gen][Threads.threadid()]
                local cache = self.t_cache[self.state_gen][Threads.threadid()]
                local decint, precal = getprecal(cache, getalpha(extdata, self, indext))
                local slc = view(f, coltuple..., indext)
                interpolate!(buf, slc, decint, precal, interp, tabmod)
                slc .= buf
            end
        else
            @threads for indext in itr
                local buf = self.t_buf[self.state_gen][Threads.threadid()]
                local cache = self.t_cache[self.state_gen][Threads.threadid()]
                local slc = view(f, coltuple..., indext)
                interpolate!(buf, slc, indbuf -> getalpha(extdata, self, indext, indbuf), interp, tabmod, cache)
                slc .= buf
            end
        end
        #        end
    elseif timeopt == SplitThreadsOpt
        #        @inbounds begin
        Threads.@threads for indth = 1:Threads.nthreads()
            local buf = self.t_buf[self.state_gen][Threads.threadid()]
            local cache = self.t_cache[self.state_gen][Threads.threadid()]
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
