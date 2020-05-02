export Advection

"""
    Advection(interpolation_type, mesh, LBC, RBC)

Creates a 1d backward semi-lagrangian advection.

- `interp`   : Interpolation type (Bspline(degree), Lagrange(degree))
- `mesh`     : UniformMesh along advection direction

"""
struct Advection

    mesh::UniformMesh
    interp::InterpolationType
    f1d::Vector{Float64}

    function Advection(mesh::UniformMesh, interp::InterpolationType)
        new(mesh, interp, zeros(mesh.length))
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
function (self::Advection)(f::Array{Float64,2}, v::Vector{Float64}, dt::Float64)

#    @sync for jchunk in Iterators.partition(1:nj, nj÷nthreads())
#        @spawn begin
            for j in eachindex(v) # jchunk
                alpha = - v[j] * dt / self.mesh.step
                interpolate!( self.f1d, view(f,:, j), alpha, self.interp)
                f[:,j] .= self.f1d
            end
#        end
#    end

end
