import VlasovBase: UniformMesh

export Advection

"""
    Advection(interpolation_type, mesh, LBC, RBC)

Creates a 1d backward semi-lagrangian advection.

- `interp`   : Interpolation type (Bspline(degree), Lagrange(degree))
- `mesh`     : UniformMesh along advection direction
- `LBC, RBC` : Boundary conditions type (:periodic, :Hermite)

"""
mutable struct Advection 
    
    interp :: InterpolationType 
    adv    :: AbstractAdvection
    dims   :: Int
    
    function Advection( interp :: InterpolationType, 
                        mesh   :: UniformMesh,
                        dims   :: Int,
                        LBC    :: Symbol,
                        RBC    :: Symbol)

        if ( (LBC,RBC) == (:periodic, :periodic))

            @assert isodd(interp.p)
            adv = BsplinePeriodic(interp, mesh)

        elseif ( (LBC, RBC) == (:Hermite, :Hermite))

            adv = BsplineHermite(interp, mesh)

        else

            throw( " This kind of advection is not yet implemented ")

        end 

        new( interp, adv, dims )

    end

    
end

"""
    advection!(f, v, dt)

Advection of a 2d function `f` discretized on a 2d `mesh`
along the input axis at velocity `v`. This function is
created from `Advector` callable type.

```julia
mesh = UniformMesh( -π, π, 64 )
advection! = Advection( Bspline(3), mesh, 1, :periodic, :periodic )

f = exp.( - mesh.points .^ 2 )

dt = 0.5
v  = ones( Float64, mesh.length)

advection!( f, v, dt )




"""
function (self :: Advection)(f  :: Array{Float64,2}, 
                             v  :: Vector{Float64}, 
                             dt :: Float64)

    p    = self.interp.p
    dims = self.dims

    if (dims == 1)
        for j in eachindex(v)
            alpha = v[j] * dt
            f[:,j] .= interpolate(f[:,j], self.adv, alpha)
        end
    else
        for i in eachindex(v)
            alpha = v[i] * dt
            f[i,:] .= interpolate(f[i,:], self.adv, alpha)
        end
    end

end
