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
import VlasovBase: UniformMesh

"""
    PeriodicAdvection(p, mesh)

Backward semi-lagrangian advection using spline interpolation.
Domain is periodic and `p` is the spline degree

"""
mutable struct PeriodicAdvection <: AbstractAdvection
    
    p        :: Int64 
    mesh     :: UniformMesh
    modes    :: Vector{Float64}
    eig_bspl :: Vector{Float64}
    eigalpha :: Vector{Complex{Float64}}
    
    function PeriodicAdvection( p, mesh )
        nx        = mesh.nx
        modes     = zeros(Float64, nx)
        modes    .= [2π * i / nx for i in 0:nx-1]
        eig_bspl  = zeros(Float64, nx)
        eig_bspl  = zeros(Float64, nx)
        eig_bspl .= bspline(p, -div(p+1,2), 0.0)
        for i in 1:div(p+1,2)-1
            eig_bspl .+= bspline(p, i-(p+1)÷2, 0.0) * 2 .* cos.(i * modes)
        end
        eigalpha  = zeros(Complex{Float64}, nx)
        new( p, mesh, modes, eig_bspl, eigalpha )
    end
    
end

function (adv :: PeriodicAdvection)(f    :: Array{Complex{Float64},2}, 
                                    v    :: Vector{Float64}, 
                                    dt   :: Float64)
    
   nx = adv.mesh.nx
   nv = length(v)
   dx = adv.mesh.dx
    
   fft!(f,1)
    
   @inbounds for j in 1:nv
      alpha = dt * v[j] / dx
      # compute eigenvalues of cubic splines evaluated at displaced points
      ishift = floor(-alpha)
      beta   = -ishift - alpha
      fill!(adv.eigalpha,0.0im)
      for i in -div(adv.p-1,2):div(adv.p+1,2)
         adv.eigalpha .+= (bspline(adv.p, i-div(adv.p+1,2), beta) 
                        .* exp.((ishift+i) * 1im .* adv.modes))
      end
          
      # compute interpolating spline using fft and properties of circulant matrices
      
      @inbounds for i in eachindex(adv.eigalpha)
         f[i,j] *= adv.eigalpha[i] / adv.eig_bspl[i]
      end
        
   end
        
   ifft!(f,1)
    
end
