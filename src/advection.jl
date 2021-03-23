
@enum TimeOptimization NoTimeOpt = 1 SimpleThreadsOpt = 2 SplitThreadsOpt = 3 MPIOpt = 4

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
struct Advection{T, N, I}
    sizeall::NTuple{N, Int}
    t_mesh::NTuple{N,UniformMesh{T}}
    t_interp::Vector{I}
    dt_base::T
    tab_coef::Vector{T}
    tab_fct::Vector{Function}
    nbsplit::Int
    mpid::Any
    function Advection(
        t_mesh::NTuple{N,UniformMesh{T}},
        t_interp::Vector{I},
        dt_base::T;
        tab_coef = [1 // 2, 1 // 1, 1 // 2],
        tab_fct = missing,
        timeopt::TimeOptimization = NoTimeOpt,
    ) where {T, N, I <: AbstractInterpolation{T}}
        length(t_interp) == N || throw(ArgumentError("size of vector of Interpolation must be equal to N=$N"))
        sizeall = length.(t_mesh)
#        v_square = dotprod(points.(t_mesh_v)) .^ 2 # precompute for ke
        mpid = timeopt == MPIOpt ? MPIData() : missing
        nbsplit = if timeopt == MPIOpt
            mpid.nb
        elseif timeopt == SplitThreadsOpt
            Threads.nthreads()
        else
            1
        end
        nfct = ismissing(tab_fct) fill(identity, size(tab_coef)), tab_fct
        return new{T,N,I}(
            sizeall,
            t_mesh,
            t_interp,
            dt_base,
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
sizeall(adv) = adv.sizeall

# Interface of external data
abstract type AbstractExtDataAdv{T,Nsum} end

"""
    Advection1dData{T,Nsp,Nv,Nsum,timeopt}
    Advection1dData(
    adv::Advection1d{T,Nsp,Nv,Nsum,timeopt}, 
    data::Array{T,Nsum},
    parext)

Mutable structure that contains variable parameters of advection series

# Type parameters
- `T::DataType` : type of data
- `Nsp` : number of space dimensions
- `Nv` : number of velocity dimensions
- `Nsum` : the total number of dimensions (Nsum = Nsp + Nv)
- `timeopt::TimeOptimization` : time optimization

# Arguments
- `adv::Advection1d{T,Nsp,Nv,Nsum}` : link to the constant data of this advection
- `data::Array{T,Nsum}` : Initial data of this advection
- `parext` : external data of this advection to compute alpha of each interpolations

# Implementation
- `adv::Advection1d{T,Nsp,Nv,Nsum,timeopt}` : link to the constant data of this advection
- `state_coef::Int` : state that is the index of `tab_coef`, it is from one to lenth(tab_coef)
- `state_dim::Int` : the dimension index, from 1 to Nsp in space states, from one to Nv in velocity state
- `data:Array{T,Nsum}` : it is the working buffer
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
mutable struct Advection1dData{T,Nsp,Nv,Nsum,timeopt}
    adv::Advection1d{T,Nsp,Nv,Nsum,timeopt}
    state_coef::Int # from 1 to length(adv.tab_coef)
    state_dim::Int # from 1 to N
    data::Array{T,Nsum}
    bufdata::Array{T}
    fmrtabdata::NTuple{Nsum,Array{T,Nsum}}
    t_buf::NTuple{Nsum,Array{T,2}}
    t_itr::Any
    tt_split::Any
    cache_alpha::Union{T,Nothing}
    cache_decint::Int64
    cache_precal::Vector{T}
    parext::AbstractExtDataAdv
    function Advection1dData(
        adv::Advection1d{T,Nsp,Nv,Nsum,timeopt},
        data::Array{T,Nsum},
        parext::AbstractExtDataAdv,
    ) where {T,Nsp,Nv,Nsum,timeopt}
        s = size(data)
        s == sizeall(adv) ||
            thrown(ArgumentError("size(data)=$s it must be $(sizeall(adv))"))
        nbthr =
            timeopt == SimpleThreadsOpt || timeopt == SplitThreadsOpt ? Threads.nthreads() :
            1
        t_buf = ntuple(x -> zeros(T, s[x], nbthr), Nsum)
        datanew = Array{T,Nsum}(undef, s)
        bufdata = Vector{T}(undef, length(data))
        fmrtabdata = ntuple(x -> zeros(T, s[getperm(parext, x)]), Nsum)
        copyto!(datanew, data)
        t_itr = ntuple(
            x -> splitvec(adv.nbsplit, CartesianIndices(s[getperm(parext, x)][2:Nsum])),
            Nsum,
        )
        t_linind = ntuple(x -> LinearIndices(s[getperm(parext, x)]), Nsum)
        fbegin(x, y) = t_linind[x][1, t_itr[x][y][1]]
        fend(x, y) = t_linind[x][end, t_itr[x][y][end]]
        tt_split = ntuple(x -> ntuple(y -> (fbegin(x, y):fend(x, y)), adv.nbsplit), Nsum)
        ordmax = maximum(get_order.((adv.t_interp_sp..., adv.t_interp_v...)))
        return new{T,Nsp,Nv,Nsum,timeopt}(
            adv,
            1,
            1,
            datanew,
            bufdata,
            fmrtabdata,
            t_buf,
            #    t_itrfirst, t_itrsecond, tt_split,
            t_itr,
            tt_split,
            nothing,
            0,
            zeros(T, ordmax + 1),
            parext,
        )
    end
end


getext(self) = self.parext
getdata(self) = self.data
getcur_t(adv::Advection1d, state_coef::Int) =
    adv.tab_fct[state_coef](adv.tab_coef[state_coef] * adv.dt_base)
getcur_t(self::Advection1dData) = getcur_t(self.adv, self.state_coef)
getstate_dim(self) = self.state_dim
isvelocity(adv::Advection1d{T,Nsp,Nv,Nsum,timeopt}, curid) where {T,Nsp,Nv,Nsum,timeopt} =
    (curid - 1) % Nsum + 1 > Nsp
isvelocitystate(state_coef::Int) = state_coef % 2 == 0
isvelocitystate(self::Advection1dData) = isvelocitystate(self.state_coef)
function getindsplit(
    self::Advection1dData{T,Nsp,Nv,Nsum,timeopt},
) where {T,Nsp,Nv,Nsum,timeopt}
    if self.adv.nbsplit != 1
        ind = timeopt == MPIOpt ? self.adv.mpid.ind : Threads.threadid()
    else
        ind = 1
    end
    return ind
end
function _getcurrentindice(self::Advection1dData{T,Nsp,Nv,Nsum}) where {T,Nsp,Nv,Nsum}
    return isvelocitystate(self) * Nsp + self.state_dim
end
getbufslgn(self::Advection1dData) = self.t_buf[_getcurrentindice(self)]
function getinterp(self::Advection1dData)
    t = isvelocitystate(self) ? self.adv.t_interp_v : self.adv.t_interp_sp
    return t[self.state_dim]
end
function getprecal(self::Advection1dData, alpha)
    if alpha != self.cache_alpha
        self.cache_alpha = alpha
        decint = convert(Int, floor(alpha))
        decfloat = alpha - decint
        self.cache_decint = decint
        #        self.cache_precal = get_precal(getinterp(self),decfloat)
        get_precal!(self.cache_precal, getinterp(self), decfloat)
    end
    return self.cache_decint, self.cache_precal
end

getitr(self) = self.t_itr[_getcurrentindice(self)][getindsplit(self)]
gett_split(self) = self.tt_split[_getcurrentindice(self)]



"""
    nextstate!(self::Advection1dData{T, Nsp, Nv, Nsum})

Function called at the end of advection function to update internal state of Advection1dData structure

# Argument
- `self::Advection1dData{T, Nsp, Nv, Nsum}` : object to update

# return value
- `ret::Bool` : `true` if the series must continue
                `false` at the end of the series.
"""
function nextstate!(self::Advection1dData{T,Nsp,Nv,Nsum}) where {T,Nsp,Nv,Nsum}
    ret = true
    if self.state_dim == [Nv, Nsp][self.state_coef%2+1]
        self.state_dim = 1
        if self.state_coef == size(self.adv.tab_coef, 1)
            self.state_coef = 1
            ret = false
        else
            self.state_coef += 1
        end
    else
        self.state_dim += 1
    end
    self.cache_alpha = nothing
    return ret
end
# default function of the interface
initcoef!(parext::AbstractExtDataAdv, self::Advection1dData) = missing
function getperm(_::AbstractExtDataAdv{T,Nsum}, curstate::Int) where {T,Nsum}
    return transposition(1, curstate, Nsum)
end
function getperm(
    _::AbstractExtDataAdv,
    advd::Advection1dData{T,Nsp,Nv,Nsum},
) where {T,Nsp,Nv,Nsum}
    return transposition(1, _getcurrentindice(advd), Nsum)
end
# this interface function must always be defined
function getalpha(parext::AbstractExtDataAdv, self::Advection1dData, ind)
    throw(error("getalpha undefined for $(typeof(parext))"))
end

# data formating
function getformdata(advd::Advection1dData{T,Nsp,Nv,Nsum}) where {T,Nsp,Nv,Nsum}
    p = getperm(getext(advd), advd)
    if p == 1:Nsum
        # the case of identity permutation no copy is needed
        f = advd.data
    else
        # ptr = pointer(advd.bufdata)
        # f = unsafe_wrap(Array, ptr, sizeall(advd.adv)[p], own=false)
        f = advd.fmrtabdata[_getcurrentindice(advd)]
        permutedims!(f, advd.data, p)
    end
    return f
end
function copydata!(
    advd::Advection1dData{T,Nsp,Nv,Nsum,timeopt},
    f,
) where {T,Nsp,Nv,Nsum,timeopt}
    if timeopt == MPIOpt && advd.adv.nbsplit != 1
        mpibroadcast(advd.adv.mpid, gett_split(advd), f)
    end
    p = getperm(getext(advd), advd)
    pinv = invperm(p)

    if f != advd.data
        permutedims!(advd.data, f, pinv)
    end
end
"""
    advection!(self::Advection1dData)

Advection1d function of a multidimensional function `f` discretized on `mesh`

# Argument
- `self::Advection1dData` : mutable structure of variables data

# Return value
- `true` : means that the advection series must continue
- `false` : means that the advection series is ended.
"""
function advection!(
    self::Advection1dData{T,Nsp,Nv,Nsum,timeopt},
) where {T,Nsp,Nv,Nsum,timeopt}
    f = self.data
    tabbuf = getbufslgn(self)

    interp = getinterp(self)
    extdata = getext(self)
    initcoef!(extdata, self)
    curind = _getcurrentindice(self)
    f = getformdata(self)
    #    @show size(f)
    tabmod = gettabmod(size(f, 1)) # just for optimization of interpolation!
    if timeopt == NoTimeOpt || timeopt == MPIOpt
        local buf = view(tabbuf, :, 1)
        #        @inbounds for ind in getitr(self)
        # fmrcpt=1
        for ind in getitr(self)
            local decint, precal = getprecal(self, getalpha(extdata, self, ind))
            local lgn = view(f, :, ind)
            interpolate!(buf, lgn, decint, precal, interp, tabmod)
            lgn .= buf
            #    fmrcpt +=1
            #     if fmrcpt == 30
            #         GC.gc()
            #     end
        end
    elseif timeopt == SimpleThreadsOpt
        #        @inbounds begin
        Threads.@threads for ind in collect(getitr(self))
            local buf = view(tabbuf, :, Threads.threadid())
            local decint, precal = getprecal(self, getalpha(extdata, self, ind))
            local lgn = view(f, :, ind)
            interpolate!(buf, lgn, decint, precal, interp, tabmod)
            lgn .= buf
        end
        #        end
    elseif timeopt == SplitThreadsOpt
        #        @inbounds begin
        Threads.@threads for indth = 1:Threads.nthreads()
            local buf = view(tabbuf, :, Threads.threadid())
            #           @inbounds for ind in getitr(self)
            for ind in getitr(self)
                local decint, precal = getprecal(self, getalpha(extdata, self, ind))
                local lgn = view(f, :, ind)
                interpolate!(buf, lgn, decint, precal, interp, tabmod)
                lgn .= buf
            end
        end
        #        end
    end
    copydata!(self, f)
    return nextstate!(self)
end
