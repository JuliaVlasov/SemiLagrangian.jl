import VlasovBase: UniformMesh

export Advection

"""
    PeriodicAdvection(p, mesh)

Backward semi-lagrangian advection using spline interpolation.
Domain is periodic and `p` is the spline degree

"""
mutable struct Advection{T} 
    
    p        :: Int64 
    adv      :: AbstractAdvection
    
    function Advection( p    :: Int, 
                        mesh :: UniformMesh,
                        LBC  :: Symbol,
                        RBC  :: Symbol)

        if ( (LBC,RBC) == (:periodic, :periodic))

            adv = BsplinePeriodic(p, mesh)

        elseif ( (LBC, RBC) == (:Hermite, :Hermite))

            adv = BsplineHermite(p, mesh)

        else

            throw( " This advection is not yet implemented ")

        end 


        new(adv, LBC, RBC)

    end

    
end

function (adv :: Advection)(f    :: Array{Float64,2}, 
                            v    :: Vector{Float64}, 
                            dt   :: Float64)
    
end
