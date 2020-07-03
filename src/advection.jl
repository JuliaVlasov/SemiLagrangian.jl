export Advection

"""
    Advection(interpolation_type, mesh, LBC, RBC)

Creates a 1d backward semi-lagrangian advection.

- `interp`   : Interpolation type (Bspline(degree), Lagrange(degree))
- `mesh`     : UniformMesh along advection direction

"""
struct Advection{T}

    mesh::UniformMesh{T}
    interp::InterpolationType
    f1d::Vector{T}
    parfft
    function Advection(mesh::UniformMesh{T}, interp::InterpolationType) where{T}
        return new{T}(mesh, interp, zeros(mesh.length), missing)
    end
    function Advection(mesh::UniformMesh{T}, interp::Bspline) where{T}
        parfft = if T == BigFloat 
            PrepareFftBig(mesh.length, T, ndims=1)
        else
            missing
        end
        return new{T}(mesh, interp, zeros(mesh.length), parfft)
    end
end

"""
    advection!(f, v, dt)

Advection of a 2d function `f` discretized on a 2d `mesh`
along the input axis at velocity `v`. This function is
created from `Advector` callable type.

```julia
mesh = UniformMesh( -π, π, 64 )
advection! = Advection( mesh, Bspline(3), :periodic )

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
            for j in eachindex(v) # jchunk
                alpha = - v[j] * dt / self.mesh.step
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
