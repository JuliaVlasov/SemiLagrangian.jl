# export Advection


include("mesh.jl")

abstract type InterpolationType{T, iscirc} end

function get_kl_ku(order)
    ku = div(order,2)
    kl = order-1-ku
    return kl, ku
end


getpermspace(N:Int)=vcat(circshif(1:N, -1), circshif((N+1):2N, -1))
getpermspeed(N:Int)=vcat(circshif(1:N, -1), (N+1):2N)
getpermswap(N:Int)=vcat((N+1):2N, 1:N)


"""
    Advection(interpolation_type, mesh, LBC, RBC)

Creates a 1d backward semi-lagrangian advection.

- `interp`   : Interpolation type (BsplineOld(degree), LagrangeOld(degree))
- `mesh`     : UniformMesh along advection direction

"""
struct Advection{T,N}
    t_mesh_x::NTuple{N, UniformMesh{T}}
    t_mesh_v::NTuple{N, UniformMesh{T}}
    t_interp_x::NTuple{N, InterpolationType{T, true}}
    t_interp_v::NTuple{N, InterpolationType{T, true}}
    pfft

    function Advection(
    t_mesh_x::NTuple{N, UniformMesh{T}},
    t_mesh_v::NTuple{N, UniformMesh{T}},
    t_interp_x::NTuple{N, InterpolationType{T}},
    t_interp_v::NTuple{N, InterpolationType{T}};
    isfft=true
    ) where{T, N}
        pfft = if isfft
            PrepareFftBig(length.(t_mesh_x), T, numdims=N, ntuple(x->x,N))
        else
            missing
        end
        return new{T,N}(t_mesh_x, t_mesh_v, t_interp_x, t_interp_v, pfft)
    end
end

mutable struct AdvectionData{T,N2,N}
    state_isspeed::Bool
    indstate::Int  # from 1 to N
    # t_space::Vector{Array{T,N2}}
    # t_speed::Vector{Array{T,N2}}
    t_space::NTuple{N, Array{T,N2}}
    t_speed::NTuple{N, Array{T,N2}}
    function AdvectionData(adv::Advection{T,N}, data::Array{T,N2}) where{T,N2,N}
        @assert N2 != 2*N "N2=$N2 must the double of N=$N"
        s = size(data)
        t_space = Vector{Array{T,N2}}(undef,N)
        t_speed = Vector{Array{T,N2}}(undef,N)
        permspace = getpermspace(N)
        permspeed = getpermspeed(N)
        permswap = getpermswap(N)
        for i in 1:N
            t_space[i] = Array{T,N2}(undef, s)
            s = s[permspace]
        end
        s = s[vcat((N+1):2N, 1:N)]
        perm = vcat(permswap)
        for i in 1:N
            t_speed =  Array{T,N2}(undef, s)
            s = s[permspeed]
        end
        copyto!(t_space[1], data)
        return new{T,N2,N}(false,1,t_space, t_speed)
    end
end
function getnewstates(
    advdata::AdvectionData{T,N2,N}; flend=false
) where{T,N2,N}
    if advdata.indstate == N
        newboolstate = !advdata.state_isspeed
        if flend 
            newboolstate = !newboolstate
        end
        return newboolstate, 1
    else
        return advdata.state_isspeed, advdata.indstate +1
    end
end
function newstate(advdata::AdvectionData{T,N2,N}; flend=false)
    (advdata.state_isspeed, advdata.indstate) = getnewstates(advdata)
end
function getindfrom(advdata::AdvectionData{T,N2,N}; flend=false) where{T,N2,N}
    return advdata.indstate + (advdata.isspeed ? N : 0)
end
function getindto(advdata::AdvectionData{T,N2,N}; flend=false) where{T,N2,N}
    ind, fl = getnewstates(advdata) 
    return advdata.indstate + (advdata.isspeed ? N : 0)
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
