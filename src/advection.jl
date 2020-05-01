export Advection

"""
    Advection(interpolation_type, mesh, LBC, RBC)

Creates a 1d backward semi-lagrangian advection.

- `interp`   : Interpolation type (Bspline(degree), Lagrange(degree))
- `mesh`     : UniformMesh along advection direction
- `LBC, RBC` : Boundary conditions type (:periodic, :Hermite)

"""
struct Advection

    mesh::UniformMesh
    interp::InterpolationType
    adv::AbstractAdvection

    function Advection(mesh::UniformMesh, interp::InterpolationType, bc::Symbol)

        new(interp, adv, dims)

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

    p = self.interp.p
    dims = self.dims
    nx, nv = size(f)

    if dims == 1
        @sync for jv in Iterators.partition(1:nv, nv÷nthreads())
            @spawn begin
                for j in jv
                    alpha = v[j] * dt
                    f[:, j] .= interpolate(f[:, j], self.adv, alpha)
                end
            end
        end
    else
        @sync for ix in Iterators.partition(1:nx, nx÷nthreads())
            @spawn begin
                for i in eachindex(ix)
                    alpha = v[i] * dt
                    f[i, :] .= interpolate(f[i, :], self.adv, alpha)
                end
            end
        end
    end

end
