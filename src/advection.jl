# export Advection


include("mesh.jl")


abstract type InterpolationType{T, iscirc} end

# TODO ne plus avoir cette fonction
function get_kl_ku(order)
    ku = div(order,2)
    kl = order-1-ku
    return kl, ku
end

function splititr(nb, lgtot)
    lg, r = divrem(lgtot,nb)
    return vcat( map(x -> ((x-1)*(lg+1)+1):x*(lg+1), 1:r), map(x -> ((x-1)*lg+r+1):(x*lg+r), (r+1):nb) )
end

function splitvec(nb, v)
    return map(x -> v[x], splititr(nb, length(v)))
end


function transperm(a,b,n)
    p = collect(1:n)
    p[a],p[b] = b, a
    return p
end

oneitr(t)=map(x->x:x,t)

addcolon(ind,tup)=(tup[1:(ind-1)]...,:,tup[ind:end]...)
function _getitr(ind, sizeitr, nb)
    indtup = vcat(1:(ind-1),(ind+1):nb)
    return addcolon.(ind, Iterators.product(sizeitr[indtup]...))
end
function _getitr(ind, sizeitr, Nsum, nbth)
    itr = _getitr(ind, sizeitr, Nsum)
    return splitvec(nbth, itr)
end
function _getitr0(ind, sizeitr, Nsum)
    indtup = transperm(1,ind,Nsum)[2:Nsum]
    return addcolon.(1, Iterators.product(sizeitr[indtup]...))
end
function _getitr0(ind, sizeitr, Nsum, nbth)
    itr = _getitr0(ind, sizeitr, Nsum)
    return splitvec(nbth, itr)
end


@enum TimeOptimization NoTimeOpt=1 SimpleThreadsOpt=2 SplitThreadsOpt=3 MPIOpt=4


"""
    Advection{T, Nsp, Nv, Nsum}
    Advection(
    t_mesh_sp::NTuple{Nsp, UniformMesh{T}},
    t_mesh_v::NTuple{Nv, UniformMesh{T}},
    t_interp_sp::NTuple{Nsp, InterpolationType{T}},
    t_interp_v::NTuple{Nv, InterpolationType{T}};
    dt_base::T;
    tab_coef=[1//2, 1//1, 1//2],
)

Immutable structure that contains constant parameters for multidimensional advection

# Type parameters
- `T::DataType` : type of data
- `Nsp` : number of space dimensions
- `Nv` : number of velocity dimensions
- `Nsum` : the total number of dimensions (Nsum = Nsp + Nv)

# Arguments
- `t_mesh_sp::NTuple{Nsp, UniformMesh{T}}` : tuple of space meshes (one per space dimension)
- `t_mesh_v::NTuple{Nv, UniformMesh{T}}` : tuple of velocity meshes (one per velocity dimension)
- `t_interp_sp::NTuple{Nsp, InterpolationType{T}}` : tuple of space interpolations (one per space dimension)
- `t_interp_v::NTuple{Nv, InterpolationType{T}}` : tuple of velocity interpolations(one per velocity dimension)
- `dt_base::T` : time delta for one advection series

# Keywords
- `tab_coef=[1//2, 1//1, 1//2]` : coefficient table for one advection series, the
    coefficients at odd indexes is for space advection series, the coefficients at even indexes is for velocity advection series

# Implementation
- sizeall : tuple of the sizes of all dimensions (space before velocity)
- sizeitr : tuple of iterators of indexes of each dimension
- t_mesh_sp : tuple of space meshes
- t_mesh_v : tuple of velocity meshes
- t_interp_sp : tuple of space interpolation types
- t_interp_v : tuple of velocity interpolation types
- dt_base::T : time unit of an advection series
- tab_coef : coefficient table

# Throws
- `ArgumentError` : `Nsp` must less or equal to `Nv`.
"""

struct Advection{T, Nsp, Nv, Nsum, timeopt}
    sizeall
    sizeitr
    t_mesh_sp::NTuple{Nsp, UniformMesh{T}}
    t_mesh_v::NTuple{Nv, UniformMesh{T}}
    t_interp_sp::NTuple{Nsp, InterpolationType{T, true}}
    t_interp_v::NTuple{Nv, InterpolationType{T, true}}
    dt_base::T
    tab_coef
    tab_fct
    v_square
    nbsplit
    mpid
    function Advection(
    t_mesh_sp::NTuple{Nsp, UniformMesh{T}},
    t_mesh_v::NTuple{Nv, UniformMesh{T}},
    t_interp_sp::NTuple{Nsp, InterpolationType{T}},
    t_interp_v::NTuple{Nv, InterpolationType{T}},
    dt_base::T;
    tab_coef=[1//2, 1//1, 1//2],
    tab_fct=[identity,identity,identity],
    timeopt::TimeOptimization=NoTimeOpt
) where{T, Nsp, Nv}
        Nsp <= Nv || thrown(ArgumentError("Nsp=$Nsp must less or equal to Nv=$Nv"))
        sizeall=length.((t_mesh_sp..., t_mesh_v...))
        Nsum = Nsp + Nv
        sizeitr = ntuple(x -> 1:sizeall[x], Nsum)
        v_square = dotprod(t_mesh_v) .^ 2 # precompute for ke
        mpid = timeopt == MPIOpt ? MPIData() : missing        
        nbsplit = if timeopt == MPIOpt
            mpid.nb
        elseif timeopt == SplitThreadsOpt
            Threads.nthreads()
        else
            1
        end
        return new{T, Nsp, Nv, Nsum, timeopt}(
    sizeall,
    sizeitr,
    t_mesh_sp, t_mesh_v, 
    t_interp_sp, t_interp_v,
    dt_base, tab_coef, tab_fct,
    v_square,
    nbsplit,
    mpid
)
    end
end
"""
    sizeall(adv::Advection)

Return a tuple of the sizes of each dimensions

# Argument
- `adv::Advection` : Advection structure.
"""
sizeall(adv)=adv.sizeall
"""
    sizeitr(adv::Advection)

Return a tuple of iterators from one to the sizes of each dimensions

# Argument
- `adv::Advection` : Advection structure.
"""
sizeitr(adv)=adv.sizeitr


# 
"""
    AdvectionData{T,Nsp,Nv,Nsum,timeopt}
    AdvectionData(
    adv::Advection{T,Nsp,Nv,Nsum,timeopt}, 
    data::Array{T,Nsum},
    parext; 
    isthread=false)

Mutable structure that contains variable parameters of advection series

# Type parameters
- `T::DataType` : type of data
- `Nsp` : number of space dimensions
- `Nv` : number of velocity dimensions
- `Nsum` : the total number of dimensions (Nsum = Nsp + Nv)

# Arguments
- `adv::Advection{T,Nsp,Nv,Nsum}` : link to the constant data of this advection
- `data::Array{T,Nsum}` : Initial data of this advection
- `parext` : external data of this advection to compute alpha of each interpolations

# Keywords
- `isthread::Bool=false` : if false only one thread is using else all possible threads are using

# Implementation
- adv : link to the constant data of this advection
- state_coef : state that is the index of tab_coef, it is from one to lenth(tab_coef)
- state_dim : the dimension index, from 1 to Nsp in space states, from one to Nv in velocity state
- data : it is the working buffer
- t_buf : tuple of buffer that is used to get the linear data for interpolation, one buffer per thread
- parext : external data of this advection to compute alpha of each interpolations

# Methods to define
- `init!(self::AdvectionData)` : this method called at the beginning of each advection to initialize parext data. The `self.parext` mutable structure is the only data that init! can modify otherwise it leads to unpredictable behaviour.
- `getalpha(self::AdvectionData, ind)` : return the alpha number that is used for interpolation. 

"""
mutable struct AdvectionData{T,Nsp,Nv,Nsum,timeopt}
    adv::Advection{T,Nsp,Nv,Nsum,timeopt}
    state_coef # from 1 to length(adv.tab_coef)
    state_dim # from 1 to N
    data::Array{T,Nsum}
    bufdata::Vector{T}
    t_buf::NTuple{Nsum, Array{T,2}}
    parext
#    itrdataind
    function AdvectionData(
    adv::Advection{T,Nsp,Nv,Nsum,timeopt}, 
    data::Array{T,Nsum},
    parext, 
) where{T,Nsp,Nv,Nsum,timeopt}
        s = size(data)
        s == sizeall(adv) || thrown(ArgumentError("size(data)=$s it must be $(sizeall(adv))"))
        nbthr = timeopt == SimpleThreadsOpt || timeopt == SplitThreadsOpt ? Threads.nthreads() : 1
#        t_buf = ntuple(x -> Array{T,2}(undef, s[x], nbthr), Nsum)
        t_buf = ntuple(x -> zeros(T, s[x], nbthr), Nsum)
        datanew = Array{T,Nsum}(undef,s)
        bufdata = Vector{T}(undef, length(data))
        copyto!(datanew, data)
        # println("trace1")

        # getview(d,ind)=(view(d,ind...), ind)
        # itrdataind = map( x-> getview.((datanew,), adv.itr[x]), 1:adv.nbth)
        # println("trace2")

        return new{T, Nsp, Nv, Nsum, timeopt}(
    adv, 1, 1,  
    datanew, bufdata, t_buf, 
    parext
)
    end
end
getext(self)=self.parext
getdata(self)=self.data
getcur_t(adv::Advection, state_coef::Int)=adv.tab_fct[state_coef](adv.tab_coef[state_coef] * adv.dt_base)
getcur_t(self::AdvectionData) = getcur_t(self.adv, self.state_coef)
getstate_dim(self)=self.state_dim
isvelocity(adv::Advection{T, Nsp, Nv, Nsum, timeopt}, curid) where {T, Nsp, Nv, Nsum, timeopt} = (curid-1)%Nsum+1 > Nsp
isvelocitystate(state_coef::Int)=state_coef%2 == 0
isvelocitystate(self::AdvectionData)=isvelocitystate(self.state_coef)
function trans1(ind, n) 
    1 < ind <= n || thrown(ArgumentException("ind=$ind n=$n we must have 1 < ind <= n"))
    return ntuple(x-> (x == ind) ? 1 : ((x == 1) ? ind : x) , n)
end
function _getcurrentindice(self::AdvectionData{T,Nsp,Nv,Nsum}) where{T,Nsp,Nv,Nsum}
    return isvelocitystate(self)*Nsp+self.state_dim
end
getbufslgn(self::AdvectionData)=self.t_buf[_getcurrentindice(self)]
function getinterp(self::AdvectionData)
    t = isvelocitystate(self) ? self.adv.t_interp_v : self.adv.t_interp_sp
    return t[self.state_dim]
end
function getindsplit(self::AdvectionData{T,Nsp, Nv, Nsum}) where{T,Nsp,Nv,Nsum}
    return isvelocitystate(self) ? Nsp : Nv
end
       
# TODO precalculer dans Avection
addcolon(ind,tup)=(tup[1:(ind-1)]...,:,tup[ind:end]...)
# function getitr(self::AdvectionData{T, Nsp, Nv, Nsum}) where {T, Nsp, Nv, Nsum}
#     ind = _getcurrentindice(self)
#     indtup = vcat(1:(ind-1),(ind+1):Nsum)
#     return addcolon.(ind, Iterators.product(sizeitr(self.adv)[indtup]...))
# end
getitr(self)=self.adv.itr[_getcurrentindice(self)]
getitr(self, indth)=self.adv.itrth[_getcurrentindice(self)][indth]


"""
    nextstate!(self::AdvectionData{T, Nsp, Nv, Nsum})

Function called at the end of advection function to update internal state of AdvectionData structure

# Argument
- `self::AdvectionData{T, Nsp, Nv, Nsum}` : object to update

# return value
- `ret::Bool` : `true`is the series must continue
                `false`at the end of the series.
"""
function nextstate!(self::AdvectionData{T, Nsp, Nv, Nsum}) where{T, Nsp, Nv, Nsum}
    ret = true
    if self.state_dim == [Nv,Nsp][self.state_coef%2+1]
        self.state_dim = 1
        if self.state_coef == size(self.adv.tab_coef,1)
            self.state_coef = 1
            ret = false
        else
            self.state_coef += 1
        end
    else
        self.state_dim += 1
    end
    return ret
end
getperm(_::Any,advd::AdvectionData{T, Nsp, Nv, Nsum}) where{T, Nsp, Nv, Nsum} = (1:Nsum)
function getdata(advd::AdvectionData{T, Nsp, Nv, Nsum}) where{T, Nsp, Nv, Nsum}
    p = getperm(getext(advd),advd)
    if p == 1:Nsum
        # the case of identity permutation no copy is needed
        f = advd.data
    else
        ptr = pointer(advd.bufdata)
        f = unsafe_wrap(Array, ptr, sizeall(advd.adv)[p])
        permutedims!(f, advd.data, p)
    end
    return f
end
function copydata!(advd::AdvectionData{T, Nsp, Nv, Nsum, timeopt}, f) where{T, Nsp, Nv, Nsum, timeopt}
    if timeopt == MPIOpt && advd.adv.nbsplit != 1
        t_split = gett_split(getext(advd), advd)
        mpibroadcast(advd.adv.mpid, t_split, f)
    end
    p = getperm(getext(advd), advd)
    pinv = invperm(p)

    if f != advd.data
        permutedims!(advd.data, f, pinv)
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
function advection!(self::AdvectionData{T,Nsp, Nv, Nsum, timeopt}) where{T,Nsp, Nv, Nsum, timeopt}
    f = self.data
    tabbuf = getbufslgn(self)

    interp = getinterp(self)
    extdata = getext(self)
    initcoef!(extdata, self)
    curind =  _getcurrentindice(self)
    f = getdata(self)
    if timeopt == NoTimeOpt || timeopt == MPIOpt
        local buf=view(tabbuf, :, 1)
        for indfirst in getitrfirst(extdata, self)
            local decint, precal = getprecal(extdata, self, indfirst)
            for ind in getitrsecond(extdata, self, indfirst)
                local lgn = view(f,ind...)
                interpolate!(buf, lgn, decint, precal, interp)
                lgn .= buf
            end
        end
    elseif timeopt == SimpleThreadsOpt
        Threads.@threads for indfirst in collect(getitrfirst(extdata, self))
            local buf=view(tabbuf, :, Threads.threadid())
            local decint, precal = getprecal(extdata, self, indfirst)
            for ind in getitrsecond(extdata, self,indfirst)
                local lgn = view(f,ind...)
                interpolate!(buf, lgn, decint, precal, interp)
                lgn .= buf
            end
        end
    elseif timeopt == SplitThreadsOpt
        Threads.@threads for indth=1:Threads.nthreads()
            local buf=view(tabbuf, :, Threads.threadid())
            for indfirst in getitrfirst(extdata, self)
                local decint, precal = getprecal(extdata, self, indfirst)
                for ind in getitrsecond(extdata, self,indfirst)
                    local lgn = view(f,ind...)
                    interpolate!(buf, lgn, decint, precal, interp)
                    lgn .= buf
                end
            end
        end
    end
    copydata!(self, f)

    return nextstate!(self)
end
"""
    compute_ke(t_mesh_sp, t_mesh_v, f)

kinetic Energie

1/2∫∫ v^2 f(x,v,t) dv dx

# Arguments
- `t_mesh_sp::NTuple{Nsp, UniformMesh{T}}` : tuple of space meshes
- `t_mesh_v::NTuple{Nv, UniformMesh{T}}` : tuple of velocity meshes
- `f::Array{T,Nsum}` : function data.
"""
function compute_ke( 
    t_mesh_sp::NTuple{Nsp, UniformMesh{T}}, 
    t_mesh_v::NTuple{Nv, UniformMesh{T}}, 
    f::Array{T,Nsum}
) where {T, Nsp, Nv, Nsum}
    Nsum == Nsp+Nv || "Nsp=$Nsp, Nv=$Nv, Nsum=$Nsum, we must have Nsum==Nsp+Nv"
    szv=length.(t_mesh_v)
    dsp = prod(step, t_mesh_sp)
    dv = prod(step, t_mesh_v)
    sum_sp = Array{T,Nv}(undef,szv)
    sum_sp .= reshape(sum(f, dims = ntuple(x->x, Nsp)), szv )
    return  (dsp * dv ) * sum( dotprod(t_mesh_v) .^ 2 .* sum_sp)
end
"""
    compute_ke(t_mesh_sp, t_mesh_v, f)

kinetic Energie

1/2∫∫ v^2 f(x,v,t) dv dx

# Arguments
- `self::AdvectionData` : mutable structure of variables data.
"""
function compute_ke(self::AdvectionData{T, Nsp, Nv, Nsum}) where {T, Nsp, Nv, Nsum}
    Nsum == Nsp+Nv || "Nsp=$Nsp, Nv=$Nv, Nsum=$Nsum, we must have Nsum==Nsp+Nv"
    adv=self.adv
    szv=length.(adv.t_mesh_v)
    dsp = prod(step, adv.t_mesh_sp)
    dv = prod(step, adv.t_mesh_v)
    sum_sp = reshape(sum(getdata(self), dims = ntuple(x->x, Nsp)), szv )
    return (dsp * dv ) * sum(adv.v_square .* sum_sp)
end  

