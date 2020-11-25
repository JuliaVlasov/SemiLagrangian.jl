# export Advection


include("mesh.jl")

abstract type InterpolationType{T, iscirc} end

# TODO ne plus avoir cette fonction
function get_kl_ku(order)
    ku = div(order,2)
    kl = order-1-ku
    return kl, ku
end

addcolon(ind,tup)=(tup[1:(ind-1)]...,:,tup[ind:end]...)
function _getitr(ind, sizeitr, Nsum)
    indtup = vcat(1:(ind-1),(ind+1):Nsum)
    return addcolon.(ind, Iterators.product(sizeitr[indtup]...))
end

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

struct Advection{T,Nsp,Nv,Nsum}
    sizeall
    sizeitr
    t_mesh_sp::NTuple{Nsp, UniformMesh{T}}
    t_mesh_v::NTuple{Nv, UniformMesh{T}}
    t_interp_sp::NTuple{Nsp, InterpolationType{T, true}}
    t_interp_v::NTuple{Nv, InterpolationType{T, true}}
    dt_base::T
    tab_coef
    v_square
    itr
    function Advection(
    t_mesh_sp::NTuple{Nsp, UniformMesh{T}},
    t_mesh_v::NTuple{Nv, UniformMesh{T}},
    t_interp_sp::NTuple{Nsp, InterpolationType{T}},
    t_interp_v::NTuple{Nv, InterpolationType{T}},
    dt_base::T;
    tab_coef=[1//2, 1//1, 1//2],
) where{T, Nsp, Nv}
        Nsp <= Nv || thrown(ArgumentError("Nsp=$Nsp must less or equal to Nv=$Nv"))
        sizeall=length.((t_mesh_sp..., t_mesh_v...))
        Nsum = Nsp + Nv
        sizeitr = ntuple(x -> 1:sizeall[x], Nsum)
        v_square = dotprod(t_mesh_v) .^ 2 # precompute for ke
        itr = ntuple(x->_getitr(x, sizeitr, Nsum),Nsum)
        return new{T, Nsp, Nv, Nsum}(
    sizeall,
    sizeitr,
    t_mesh_sp, t_mesh_v, 
    t_interp_sp, t_interp_v,
    dt_base, tab_coef,
    v_square,
    itr
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
    AdvectionData{T,Nsp,Nv,Nsum}
    AdvectionData(
    adv::Advection{T,Nsp,Nv,Nsum}, 
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
mutable struct AdvectionData{T,Nsp,Nv,Nsum,isthr}
    adv::Advection{T,Nsp,Nv,Nsum}
    state_coef # from 1 to length(adv.tab_coef)
    state_dim # from 1 to N
    data::Array{T,Nsum}
    t_buf::NTuple{Nsum, Array{T,2}}
    parext
    function AdvectionData(
    adv::Advection{T,Nsp,Nv,Nsum}, 
    data::Array{T,Nsum},
    parext; 
    isthread::Bool=false
) where{T,Nsp,Nv,Nsum,isthr}
        s = size(data)
        s == sizeall(adv) || thrown(ArgumentError("size(data)=$s it must be $(sizeall(adv))"))
        nbthr = isthread ? Threads.nthreads() : 1
#        t_buf = ntuple(x -> Array{T,2}(undef, s[x], nbthr), Nsum)
        t_buf = ntuple(x -> zeros(T, s[x], nbthr), Nsum)
        datanew = Array{T,Nsum}(undef,s)
        copyto!(datanew, data)
        return new{T,Nsp, Nv, Nsum, isthread}(
    adv, 1, 1,  
    datanew, t_buf, 
    parext
)
    end
end
getext(self)=self.parext
getdata(self)=self.data
getcur_t(adv::Advection, state_coef::Int)=adv.tab_coef[state_coef] * adv.dt_base
getcur_t(self::AdvectionData) = getcur_t(self.adv, self.state_coef)
getstate_dim(self)=self.state_dim
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
       
# TODO precalculer dans Avection
addcolon(ind,tup)=(tup[1:(ind-1)]...,:,tup[ind:end]...)
# function getitr(self::AdvectionData{T, Nsp, Nv, Nsum}) where {T, Nsp, Nv, Nsum}
#     ind = _getcurrentindice(self)
#     indtup = vcat(1:(ind-1),(ind+1):Nsum)
#     return addcolon.(ind, Iterators.product(sizeitr(self.adv)[indtup]...))
# end
getitr(self)=self.adv.itr[_getcurrentindice(self)]

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
"""
    advection!(self::AdvectionData)

Advection function of a multidimensional function `f` discretized on `mesh`

# Argument
- `self::AdvectionData` : mutable structure of variables data

# Return value
- `true` : means that the advection series must continue
- `false` : means that the advection series is ended.
"""
function advection!(self::AdvectionData{T,Nsp, Nv, Nsum, isthr}) where{T,Nsp, Nv, Nsum, isthr}
    f = self.data
    tabbuf = getbufslgn(self)
#    tabbuf = ntuple(x->zeros(T, size(self.data)[_getcurrentindice(self)]), Threads.nthreads())
#    tabbuf = zeros(T,size(self.data)[_getcurrentindice(self)], Threads.nthreads())
    interp = getinterp(self)
    init!(self)
    global cl_obs
    clockbegin(cl_obs,_getcurrentindice(self))
#     if !self.isthread
# #        fl=true
#         buf = tabbuf[:,1]
#         for ind in getitr(self)
#             println("ind=$ind")
# #            buf=tabbuf[:,1]
#             alpha = getalpha(self, ind)
#             interpolate!(buf, f[ind...], alpha, interp)
#             f[ind...] .= buf
#         end
#     end
#    else
#        Threads.@threads :static for ind in getitr(self)
    if isthr
        Threads.@threads for ind in getitr(self)
            alpha = getalpha(self, ind)
            bufth = tabbuf[:,Threads.threadid()]
            # println("av ind=$ind alpha=$alpha norm(lgn)=$(norm(f333[ind...])) thrid=$(Threads.threadid())")
            # @assert size(buf) == size(f333[ind...]) "PB TAILLE !!!!"
            interpolate!(bufth, f[ind...], alpha, interp)
            f[ind...] .= bufth
            # println("ap ind=$ind alpha=$alpha norm(lgn)=$(norm(f333[ind...])) thrid=$(Threads.threadid()) norm(buf)=$(norm(buf))")
        end
    else
        buf=tabbuf[:,1]
        f=self.data
        println("trace")
        for ind in getitrbis(self)
            alpha = getalpha(self, ind)
            interpolate!(buf, f[ind...], alpha, interp)
            f[ind...] .= buf
        end
        @assert false "ON est passé par là"
    end
#    end
    clockend(cl_obs,_getcurrentindice(self))

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

