export Advection

"""
    Advection(interpolation_type, mesh, LBC, RBC)

Creates a 2d backward semi-lagrangian advection.

- `interp`   : Interpolation type (Bspline(degree), Lagrange(degree))
- `mesh`     : UniformMesh along advection direction

"""
struct Advection2d{T}

    mesh_x::UniformMesh{T}
    mesh_y::UniformMesh{T}
    interp::InterpolationType
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
function advection!(self::Advection2d{T}, f::Array{T,2}, v::Array{Tuple{T,T}, 2}, dt::T) where {T}
    buf = zeros(T,size(f));
#    @sync for jchunk in Iterators.partition(1:nj, nj÷nthreads())
#        @spawn begin
    tabv = fill((zero(T),zero(T)), size(v))
    for j in eachindex(v) # jchunk
        tabv[j] = (- v[j][1] * dt / self.mesh_y.step, - v[j][2] * dt / self.mesh_x.step)
    end
    interpolate!( self, buf, f, tabv, self.interp)
    f .= buf
#        end
#    end
end
