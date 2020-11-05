# export Advection


include("mesh.jl")

abstract type InterpolationType{T, iscirc} end

function get_kl_ku(order)
    ku = div(order,2)
    kl = order-1-ku
    return kl, ku
end


# getpermspace(N:Int)=vcat(circshif(1:N, -1), circshif((N+1):2N, -1))
# getpermspeed(N:Int)=vcat(circshif(1:N, -1), (N+1):2N)
# getpermswap(N:Int)=vcat((N+1):2N, 1:N)


"""
    Advection(interpolation_type, mesh, LBC, RBC)

Creates a 1d backward semi-lagrangian advection.

- `interp`   : Interpolation type (BsplineOld(degree), LagrangeOld(degree))
- `mesh`     : UniformMesh along advection direction

"""

struct Advections{T,N}
    t_mesh_x::NTuple{N, UniformMesh{T}}
    t_mesh_v::NTuple{N, UniformMesh{T}}
    t_interp_x::NTuple{N, InterpolationType{T, true}}
    t_interp_v::NTuple{N, InterpolationType{T, true}}
    dt_base::T
    tab_coef
    pfftbig
    dimsperm
    swapperm
    bascperm
    function Advections(
    t_mesh_x::NTuple{N, UniformMesh{T}},
    t_mesh_v::NTuple{N, UniformMesh{T}},
    t_interp_x::NTuple{N, InterpolationType{T}},
    t_interp_v::NTuple{N, InterpolationType{T}};
    dt_base::T,
    tab_coef=[1//2, 1//1, 1//2],
    isfftbig=true ) where{T, N}
        pfftbig = if isfftbig
            PrepareFftBig(length.(t_mesh_x), T, numdims=N, ntuple(x->x,N))
        else
            missing
        end
        dimsperm = vcat(circshift(1:N,-1),(N+1):2N)
        swapperm = vcat((N+1):2N,1:N)
        bascperm = swapperm[dimsperm]
        return new{T,N}(
    t_mesh_x, t_mesh_v, 
    t_interp_x, t_interp_v, 
    dt_base, tab_coef,
    pfftbig,
    dimsperm, swapperm, bascperm
)
    end
end
function getworkbufs(s, nbthreads, T::DataType)
    res = Vector{Vector{T}}(undef,nbthreads)
    for i=1:nbthreads
        res[i] = Vector{T}(undef, s)
    end
    res
end

mutable struct AdvectionsData{T,N2,N}
    adv::Advections{T,N}
    state_coef # from 1 to length(adv.tab_coef)
    state_dim # from 1 to N
    t_data::NTuple{N2, Array{T,N2}}
    rho::Array{T,N}
    t_elfield::NTuple{N, Array{T,N}}
    t_buf::Ntuple{N2,Vector{Vector{T}}}
    function AdvectionsData(
    adv::Advections{T,N}, 
    data::Array{T,N2}; 
    nbthreads=1
) where{T,N2,N}
        @assert N2 != 2*N "N2=$N2 must the double of N=$N"
        s = size(data)
        s_x = s[1:N]
        tabsize =Vector{Vector{Int}}(undef, N2)
        svar = s
        for i=1:N2
            tabsize[i] = svar
            svar = svar[adv.dimsperm]
            if i == N
                svar = svar[adv.swapperm]
            end
        end
        t_data = ntuple( x-> Array{T,N2}(undef, tabsize[x]), N2)
        rho = Array{T,N}(undef, s_x)
        t_elffield = ntuple( x-> Array{T,N}(undef, s_x), N)
        t_buf = ntuple(x->getworkbuf(size(t_data[x],1), nbthreads, T), N2)


        copyto!(t_data[1], data)
        return new{T,N2,N}(adv, 1, 1, t_data, rho, t_elfield, t_buf)
    end
end
function _getcurrentindice(self::AdvectionsData{T,N2,N}) where{T,N2,N}
    ((self.state_coef-1)*N+self.state_dim-1)%N2+1
end
function _getnextindice(self::AdvectionsData{T,N2,N}) where{T,N2,N}
    if self.state_coef == size(self.adv.tab_coef,1) && self.state_dim == N
        return 1
    else
        return _getcurrentindice(self)%N2+1
    end
end
getdata(self::AdvectionsData)=self.t_data[_getcurrentindice(self)]
getbufslgn(self::AdvectionsData)=self.t_buf[_getcurrentindice(self)]
isspeedstate(self::AdvectionsData)=self.state_coef%2 == 0
function getinterp(self::AdvectionsData)
    if isspeedstate(self)
        return self.adv.t_interp_v[self.state_dim]
    else
        return self.adv.t_interp_x[self.state_dim]
    end
end
       
function getalpha(self::AdvectionsData{T, N2, N}, ind) where{T,N2, N}
    coef = - self.adv.dt_base * self.adv.tab_coef[self.state_coef]
    if isspeedstate(self)
        return coef * self.t_elffield[self.state_dim][ind...]
    else
        mesh = self.adv.t_mesh_v[self.state_dim]
        return  coef/mesh.step*mesh.points[ind]
    end
end
function getitrfirst(self::AdvectionsData)
    if isspeedstate(self)
        s = size(self.t_elffield[1])
        return Iterators.product((1:sz for sz in s)...)
    else
        return self.adv.t_mesh[self.state_dim].points
    end
end
function getitrsecond(self::AdvectionsData, indfirst)
    s = size(getdata(self))
    tupbegin=(1:s[i] for i=2:N)
    tupend = if isspeedstate(self)
        (ind:ind for ind in indfirst)
    else
        (((i == self.state_dim) ? (indfirst:indfirst) : 1:s[N+i]) for i=1:N)
    end
    return Iterators.product((tupdeb...,tupend...)...)
end
function nextstate!(self::AdvectionsData{T, N2, N}) where{T,N2, N}
    perm = self.adv.dimsperm
    if self.state_dims == N
        self.state_dims = 1
        if self.state_coef == size(self.adv.tab_coef,1)
            self.state_coef = 1
        else
            perm = self.adv.bascperm
            self.state_coef += 1
        end
    else
        self.state_dims += 1
    end
    return perm
end

    
function changestate!(self::AdvectionsData)
    curbuf = getbuf(self)
    perm = nextstate!(self)
    newbuf = getbuf(self)
    permutedims!(newbuf, curbuf, perm)
    if self.state_dim == 1 && ispeedstate(self)
        compute_charge!(self.rho, self.adv.t_mesh_v, newbuf)
        compute_elfield!(self.t_elf, self.adv.t_mesh_x, self.rho, self.adv.pfft)
    end
    return self.state_dim != 1 || self.state_coef != 1
end

    




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
function advection!(self::AdvectionsData{T, N2, N}) where{T,N2, N}
    f = getdata(self)
    tabbuf = getbufslgn(self)
    interp = getinterp(self)
    for indfirst in getitrfirst(self)
        alpha = getalpha(self,indfirst)
        buf=tabbuf[1]
        for indsecond in getitrsecond(self, indfirst)
            indglob=(:,indsecond...)
            interpolate!(buf,f[indglob...], alpha, interp)
            f[indglob...] .= buf
        end
    end
    return changestate!(self)
end


function advection!(self::Advection{T, N},
     f::Array{T}, 
     v::Union{Tuple{Vector{T},N},Array{T,N}}, 
     dt::T,
     dims::Ntuple{N,Int}
     dimsother::Ntuple{N,Int}
) where {T,N,N2}
    maxalpha = 0
    minalpha = 100000
    for (inddim, dim) in enumerate(dims) 
        tabbuf = Vector{Vector{T}}(undef, Threads.nthreads())
        for i=1:size(tabbuf,1)
            tabbuf[i] = Vector{T}(undef,size(f,dim));
        end
        Threads.@threads for j in CartesianIndices()=1:size(v,1)
            value = v[j]
        #for (j, value) in enumerate(v) # jchunk
            buf = tabbuf[Threads.threadid()]
            #                println("value=$value dt=$dt step =$(self.mesh.step)")
            alpha = - value * dt / self.mesh.step
    #                println("j=$j alpha=$alpha")
            maxalpha = max(maxalpha,abs(alpha))
            minalpha = min(minalpha,abs(alpha))
            interpolate!( buf, f[:, j], alpha, self.interp)
            f[:,j] .= buf
        end
    end

    # nbthr = Threads.nthreads()

    # @sync for numtr=1:nbthr
    #     @spawn begin
    #         buf = Vector{T}(undef,size(f,1));
    #         for j=numtr:nbthr:size(v,1)
    #             value = v[j]
    #             alpha = - value * dt / self.mesh.step
    #             interpolate!( buf, f[:, j], alpha, self.interp)
    #             f[:,j] .= buf
    #         end
    #     end
    # end

    #    @sync for jchunk in Iterators.partition(1:nj, nj÷nthreads())
#        @spawn begin

#        end
#    end
#    println("maxaplha = $maxalpha minaplha = $minalpha")
end

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
