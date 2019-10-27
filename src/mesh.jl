import Statistics: mean

export UniformMesh

"""

    UniformMesh(start, stop, length)

1D uniform mesh data.

length   : number of points
length-1 : number of cells

If you want remove the last point for periodic domain, set endpoint=false

    - `step`  : size step
    - `points`: Array with node positions
    - `width` : Distance between left and right edges.

"""
struct UniformMesh

   start    :: Float64
   stop     :: Float64
   length   :: Int
   step     :: Float64
   points   :: Vector{Float64}
   width    :: Float64

   function UniformMesh(start, stop, length::Int; endpoint=true)

       if (endpoint)
           points = range(start, stop=stop, length=length)
       else
           points = range(start, stop=stop, length=length+1)[1:end-1]
       end

       step = points.step

       width = stop - start

       new( start, stop, length, step, points, width)

   end

end

export compute_charge!

"""
    compute_rho!( rho, mesh_v, fvx)
 
 Compute charge density

 ρ(x,t) = ∫ f(x,v,t) dv

"""
function compute_charge!(rho   :: Vector{Float64}, 
                         meshv :: UniformMesh,
                         fvx   :: Array{Float64,2})

    dv   =  meshv.step
    rho .=  dv .* vec(sum(fvx, dims=1))
    rho .= rho .- mean(rho)

end
 

export compute_e!

"""

    compute_e!( e, mesh, ρ)

    ∇.e = - ρ

Inplace computation of electric field. Fields e and rho are
already allocated.

"""
function compute_e!(e     :: Vector{Float64},
		    meshx :: UniformMesh,
		    rho   :: Vector{Float64})

   nx = meshx.length
   k = 2π / (meshx.stop - meshx.start)
   modes  = zeros(Float64, nx)
   modes .= k .* vcat(0:nx÷2-1,-nx÷2:-1)
   modes[1] = 1.0
   e .= real(ifft(-1im .* fft(rho) ./ modes))

end
