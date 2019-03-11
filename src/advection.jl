import VlasovBase: UniformMesh

export Advection

"""
    Advection(interpolation_type, mesh, LBC, RBC)

Backward semi-lagrangian advection using spline interpolation.
Domain is periodic and `p` is the spline degree

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

            adv = BsplinePeriodic(interp, mesh)

        elseif ( (LBC, RBC) == (:Hermite, :Hermite))

            adv = BsplineHermite(interp, mesh)

        else

            throw( " This advection is not yet implemented ")

        end 


        new( interp, adv, dims )

    end

    
end

"""
    advection! = Advection(interpolation_type, mesh, LBC, RBC)

Advection of a 2d function `f` discretized on a 2d `mesh`
along the input axis at velocity `v`

"""
function (self :: Advection)(f    :: Array{Float64,2}, 
                             v    :: Vector{Float64}, 
                             dt   :: Float64)

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
