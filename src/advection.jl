# export Advection


include("mesh.jl")

abstract type InterpolationType{T, iscirc} end

function get_kl_ku(order)
    ku = div(order,2)
    kl = order-1-ku
    return kl, ku
end





"""
    Advection(interpolation_type, mesh, LBC, RBC)

Creates a 1d backward semi-lagrangian advection.

- `interp`   : Interpolation type (BsplineOld(degree), LagrangeOld(degree))
- `mesh`     : UniformMesh along advection direction

"""
struct Advection{T}

    mesh::UniformMesh{T}
    interp::InterpolationType{T}
    f1d::Vector{T}
    function Advection(mesh::UniformMesh{T}, interp::InterpolationType{T}) where{T}
        return new{T}(mesh, interp, zeros(mesh.length))
    end
    # function Advection(mesh::UniformMesh{T}, interp::BsplineOld) where{T}
    #     parfft = if T == BigFloat 
    #         PrepareFftBig(mesh.length, T, ndims=1)
    #     else
    #         missing
    #     end
    #     return new{T}(mesh, interp, zeros(mesh.length), parfft)
    # end
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
function advection!(self::Advection{T}, f::Array{T,2}, v::Vector{T}, dt::T) where {T}
    buf = zeros(T,size(f,1));
#    @sync for jchunk in Iterators.partition(1:nj, nj÷nthreads())
#        @spawn begin
    maxalpha = 0
    minalpha = 100000
            for (j, value) in enumerate(v) # jchunk
                alpha = - value * dt / self.mesh.step
#               println("j=$j alpha=$alpha")
                maxalpha = max(maxalpha,abs(alpha))
                minalpha = min(minalpha,abs(alpha))
                interpolate!( self, buf, f[:, j], alpha, self.interp)
                f[:,j] .= buf
            end
#        end
#    end
#    println("maxaplha = $maxalpha minaplha = $minalpha")
end

