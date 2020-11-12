# export Advection


include("mesh.jl")

abstract type InterpolationType{T, iscirc} end


# getpermspace(N:Int)=vcat(circshif(1:N, -1), circshif((N+1):2N, -1))
# getpermspeed(N:Int)=vcat(circshif(1:N, -1), (N+1):2N)
# getpermswap(N:Int)=vcat((N+1):2N, 1:N)


"""
    Advection(interpolation_type, mesh, LBC, RBC)

Creates a 1d backward semi-lagrangian advection.

- `interp`   : Interpolation type (BsplineOld(degree), LagrangeOld(degree))
- `mesh`     : UniformMesh along advection direction

"""

struct Advection{T,Nsp,Nv,Nsum}
    t_mesh_sp::NTuple{Nsp, UniformMesh{T}}
    t_mesh_v::NTuple{Nv, UniformMesh{T}}
    t_mesh_all::NTuple{Nsum, UniformMesh{T}}
    t_interp_sp::NTuple{Nsp, InterpolationType{T, true}}
    t_interp_v::NTuple{Nv, InterpolationType{T, true}}
    t_interpall::NTuple{Nsum, InterpolationType{T, true}}
    dt_base::T
    tab_coef
    function Advection(
    t_mesh_sp::NTuple{Nsp, UniformMesh{T}},
    t_mesh_v::NTuple{Nv, UniformMesh{T}},
    t_interp_sp::NTuple{Nsp, InterpolationType{T}},
    t_interp_v::NTuple{Nv, InterpolationType{T}};
    dt_base::T,
    tab_coef=[1//2, 1//1, 1//2],
) where{T, Nsp, Nv}
        # Nsp+Nv == Nsum || thrown(ArgumentError("Nsp=$Nsp Nv=$Nv Nsum=$Nsum Nsp+Nv must be equal Nsum"))
        Nsp <= Nv || thrown(ArgumentError("Nsp=$Nsp must less or equal to Nv=$Nv"))

        t_meshall=(t_mesh_sp..., t_mesh_v...)
        t_interpall=(t_mesh_sp..., t_mesh_v...)
        # pfftbig = if isfftbig
        #     PrepareFftBig(length.(t_mesh_sp), T, numdims=N, ntuple(x->x,N))
        # else
        #     missing
        # end
        # dimsperm = vcat(circshift(1:N,-1),(N+1):2N)
        # swapperm = vcat((N+1):2N,1:N)
        # bascperm = swapperm[dimsperm]
        Nsum = Nsp+nv
        return new{T, Nsp, Nv, Nsum}(
    t_mesh_sp, t_mesh_v, t_meshall, 
    t_interp_sp, t_interp_v, t_interpall,
    dt_base, tab_coef,
 )
    end
end
function mkworkbufs(s, nbthreads, T::DataType)
    res = Vector{Vector{T}}(undef,nbthreads)
    for i=1:nbthreads
        res[i] = Vector{T}(undef, s)
    end
    res
end

mutable struct AdvectionData{T,Nsp,Nv,Nsum}
    adv::Advection{T,Nsp,Nv,Nsum}
    state_coef # from 1 to length(adv.tab_coef)
    state_dim # from 1 to N
    data::Array{T,Nsum}
    t_buf::NTuple{Nsum, Vector{Vector{T}}}
    parext
    isthread   
    function AdvectionData(
    adv::Advection{T,Nsp,Nv,Nsum}, 
    data::Array{T,Nsum},
    parext; 
    isthread=false
) where{T,Nsp,Nv,Nsum}
        s = size(data)
        s == length.(adv.t_meshall) || thrown(ArgumentError("size(data)=$s it must be $(length.(adv.t_meshall))"))
        nbthr = isthread ? Threads.nthread() : 1
        t_buf = ntuple(x ->mkworkbufs(length(t_mesall[x]), nbthr, T), Nsum)
        # s_sp = s[1:N]
        # tabsize =Vector{Vector{Int}}(undef, N2)
        # svar = s
        # for i=1:N2
        #     tabsize[i] = svar
        #     svar = svar[adv.dimsperm]
        #     if i == N
        #         svar = svar[adv.swapperm]
        #     end
        # end
        # t_data = ntuple( x-> Array{T,N2}(undef, tabsize[x]), N2)
        # rho = Array{T,N}(undef, s_sp)
        # t_elffield = ntuple( x-> Array{T,N}(undef, s_sp), N)
        # t_buf = ntuple(x->getworkbuf(size(t_data[x],1), nbthreads, T), N2)

        datanew = Array{T,Nsum}(undef,s)
        copyto!(datanew, data)
        return new{T,Nsp, Nv, Nsum}(
    adv, 1, 1,  
    datanew, t_buf, 
    parext, isthread
)
    end
end
# function getalpha(self::AdvectionData{T, Nsp, Nv, Nsum}, ind) where{T, Nsp, Nv, Nsum}
#     thrown(MethodError("Method getalpha must be defined"))
# end
getext(self)=self.parext
getcur_t(self) = self.adv.tab_coef[self.state_coef] * self.adv.dt_base
getstate_dim(self)=self.state_dim
isvelocitystate(self::AdvectionData)=self.state_coef%2 == 0
function _getcurrentindice(self::AdvectionData{T,Nsp,Nv,Nsum}) where{T,Nsp,Nv,Nsum}
    return isvelocitystate(self)*Nsp+self.state_dim
end
# function _getnextindice(self::AdvectionData{T,N2,N}) where{T,N2,N}
#     if self.state_coef == size(self.adv.tab_coef,1) && self.state_dim == N
#         return 1
#     else
#         return _getcurrentindice(self)%N2+1
#     end
# end
# getdata(self::AdvectionData)=self.t_data[_getcurrentindice(self)]
getbufslgn(self::AdvectionData)=self.t_buf[_getcurrentindice(self)]
function getinterp(self::AdvectionData)
    if isvelocitystate(self)
        return self.adv.t_interp_v[self.state_dim]
    else
        return self.adv.t_interp_sp[self.state_dim]
    end
end
       
# function getalpha(self::AdvectionData{T, N2, N}, ind) where{T,N2, N}
#     coef = - self.adv.dt_base * self.adv.tab_coef[self.state_coef]
#     if isspeedstate(self)
#         return coef * self.t_elffield[self.state_dim][ind...]
#     else
#         mesh = self.adv.t_mesh_v[self.state_dim]
#         return  coef/mesh.step*mesh.points[ind]
#     end
# end
addcolon(ind,tup)=(tup[1:(ind-1)]...,:,tup[ind:end]...)
function getitr(self)
    ind = _getcurrentindice(self)
    indtup = vcat(1:(ind-1),(ind+1):N2)
    return addcolon.(Iterator.product(lentgh.(t_meshall)[indtup]))
end

# function getitrfirst(self::AdvectionData)
#     if isspeedstate(self)
#         s = size(self.t_elffield[1])
#         return Iterators.product((1:sz for sz in s)...)
#     else
#         return self.adv.t_mesh[self.state_dim].points
#     end
# end
# function getitrsecond(self::AdvectionData, indfirst)
#     s = size(getdata(self))
#     tupbegin=(1:s[i] for i=2:N)
#     tupend = if isspeedstate(self)
#         (ind:ind for ind in indfirst)
#     else
#         (((i == self.state_dim) ? (indfirst:indfirst) : 1:s[N+i]) for i=1:N)
#     end
#     return Iterators.product((tupdeb...,tupend...)...)
# end
function nextstate!(self::AdvectionData{T, N2, N}) where{T,N2, N}
    ret = true
    if self.state_dims == [Nv,Nsp][self.state_coef%2+1]
        self.state_dims = 1
        if self.state_coef == size(self.adv.tab_coef,1)
            self.state_coef = 1
            ret = false
        else
            self.state_coef += 1
         end
    else
        self.state_dims += 1
    end
    return ret
 end

    
# function changestate!(self::AdvectionData)
# #    curbuf = getbuf(self)
# #    newbuf = getbuf(self)
# #    permutedims!(newbuf, curbuf, perm)
#     if nextstate!(self)
#         self.advextinit!(self.parext, self.data)
#         # compute_charge!(self.rho, self.adv.t_mesh_v, newbuf)
#         # compute_elfield!(self.t_elf, self.adv.t_mesh_sp, self.rho, self.adv.pfft)
#     end
#     return self.state_dim != 1 || self.state_coef != 1
# end

    




"""
    advection!(f, v, dt)

Advection of a 2d function `f` discretized on a 2d `mesh`
along the input axis at velocity `v`. This function is
created from `Advector` callable type.

```julia
mesh = UniformMesh( -π, π, 64 )
advection! = Advection( mesh, BsplineOld(3), :periodic )

f = exp.( - mesh.points .^ 2 )

dt = 0.5
v  = ones( Float64, mesh.length)

advection!( f, v, dt )

"""
addlgn(x)=(:,x...)
function advection!(self::AdvectionData{T, N2, N}) where{T,N2, N}
    f = self.data
    tabbuf = getbufslgn(self)
    interp = getinterp(self)
    init!(self)
    if isthread
        for ind in getitr(self)
            alpha = getalpha(self, ind)
            buf=tabbuf[1]
            interpolate!(buf, f[ind...], alpha, interp)
            f[ind...] .= buf
        end
    else
        Threads.@threads for ind in getitr(self)
            alpha = getalpha(self, ind)
            buf=tabbuf[Threads.getthreadid()]
            interpolate!(buf, f[ind...], alpha, interp)
            f[ind...] .= buf
        end
    end
    return nextstate!(self)
end


# function advection!(self::Advection{T, N},
#      f::Array{T}, 
#      v::Union{Tuple{Vector{T},N},Array{T,N}}, 
#      dt::T,
#      dims::Ntuple{N,Int}
#      dimsother::Ntuple{N,Int}
# ) where {T,N,N2}
#     maxalpha = 0
#     minalpha = 100000
#     for (inddim, dim) in enumerate(dims) 
#         tabbuf = Vector{Vector{T}}(undef, Threads.nthreads())
#         for i=1:size(tabbuf,1)
#             tabbuf[i] = Vector{T}(undef,size(f,dim));
#         end
#         Threads.@threads for j in CartesianIndices()=1:size(v,1)
#             value = v[j]
#         #for (j, value) in enumerate(v) # jchunk
#             buf = tabbuf[Threads.threadid()]
#             #                println("value=$value dt=$dt step =$(self.mesh.step)")
#             alpha = - value * dt / self.mesh.step
#     #                println("j=$j alpha=$alpha")
#             maxalpha = max(maxalpha,abs(alpha))
#             minalpha = min(minalpha,abs(alpha))
#             interpolate!( buf, f[:, j], alpha, self.interp)
#             f[:,j] .= buf
#         end
#     end

#     # nbthr = Threads.nthreads()

#     # @sync for numtr=1:nbthr
#     #     @spawn begin
#     #         buf = Vector{T}(undef,size(f,1));
#     #         for j=numtr:nbthr:size(v,1)
#     #             value = v[j]
#     #             alpha = - value * dt / self.mesh.step
#     #             interpolate!( buf, f[:, j], alpha, self.interp)
#     #             f[:,j] .= buf
#     #         end
#     #     end
#     # end

#     #    @sync for jchunk in Iterators.partition(1:nj, nj÷nthreads())
# #        @spawn begin

# #        end
# #    end
# #    println("maxaplha = $maxalpha minaplha = $minalpha")
# end

# old one dimension version

# function advection!(self::Advection{T}, f::Array{T,2}, v::Vector{T}, dt::T) where {T}
#     maxalpha = 0
#     minalpha = 100000
#     tabbuf = Vector{Vector{T}}(undef, Threads.nthreads())
#     for i=1:size(tabbuf,1)
#         tabbuf[i] = Vector{T}(undef,size(f,1));
#     end

#     Threads.@threads for j=1:size(v,1)
#         value = v[j]
#     #for (j, value) in enumerate(v) # jchunk
#         buf = tabbuf[Threads.threadid()]
#         #                println("value=$value dt=$dt step =$(self.mesh.step)")
#         alpha = - value * dt / self.mesh.step
# #                println("j=$j alpha=$alpha")
#         maxalpha = max(maxalpha,abs(alpha))
#         minalpha = min(minalpha,abs(alpha))
#         interpolate!( buf, f[:, j], alpha, self.interp)
#         f[:,j] .= buf
#     end

#     # nbthr = Threads.nthreads()

#     # @sync for numtr=1:nbthr
#     #     @spawn begin
#     #         buf = Vector{T}(undef,size(f,1));
#     #         for j=numtr:nbthr:size(v,1)
#     #             value = v[j]
#     #             alpha = - value * dt / self.mesh.step
#     #             interpolate!( buf, f[:, j], alpha, self.interp)
#     #             f[:,j] .= buf
#     #         end
#     #     end
#     # end

#     #    @sync for jchunk in Iterators.partition(1:nj, nj÷nthreads())
# #        @spawn begin

# #        end
# #    end
# #    println("maxaplha = $maxalpha minaplha = $minalpha")
# end
