import VlasovBase: UniformMesh

export Advection

"""
    PeriodicAdvection(p, mesh)

Backward semi-lagrangian advection using spline interpolation.
Domain is periodic and `p` is the spline degree

"""
mutable struct Advection{T} 
    
    p        :: Int64 
    mesh     :: UniformMesh
    
    function Advection{T}( p         :: Int, 
                           mesh      :: UniformMesh,
                           dimension :: Int, 
                           LBC       :: Symbol,
                           RBC       :: Symbol) where {T<: Union{Float64, ComplexF64}}

        new(p, mesh, dimension, LBC, RBC)

    end

    
end

function (adv :: Advection)(f    :: Array{ComplexF64,2}, 
                            v    :: Vector{Float64}, 
                            dt   :: Float64)
    
end

function (adv :: Advection)(f    :: Array{Float64,2}, 
                            v    :: Vector{Float64}, 
                            dt   :: Float64)
    
end
