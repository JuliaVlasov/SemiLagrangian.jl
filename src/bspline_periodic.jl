using FFTW, LinearAlgebra

export Bspline

struct Bspline <: InterpolationType

    p :: Int

    function Bspline( p :: Int )
	
	    @assert (p & 1 == 1)
	    new( p )

    end
end


"""
    bspline(p, j, x)

Return the value at x in [0,1[ of the B-spline with
integer nodes of degree p with support starting at j.
Implemented recursively using the de Boor's recursion formula

Original python program written by 
Eric Sonnendrucker (Max-Planck-Institut fur Plasmaphysik - Garching (Germany))

"""
function bspline(p::Int, j::Int, x::Float64)

   if p == 0
       if j == 0
           return 1.0
       else
           return 0.0
       end
   else
       w = (x - j) / p
       w1 = (x - j - 1) / p
   end
   ( w       * bspline(p - 1, j    , x) +
    (1 - w1) * bspline(p - 1, j + 1, x))
end


mutable struct BsplinePeriodic <: AbstractAdvection

    p        :: Int
    mesh     :: UniformMesh
    modes    :: Vector{ComplexF64}
    eig_bspl :: Vector{ComplexF64}
    eigalpha :: Vector{ComplexF64}

    function BsplinePeriodic( bspl :: Bspline, 
                              mesh :: UniformMesh )

        n = mesh.length
        modes = 2π .* (0:n-1) ./ n
        eig_bspl = zeros(Complex{Float64},n)
        eigalpha = zeros(Complex{Float64},n)
        eig_bspl .= bspline(bspl.p, -div(bspl.p+1,2), 0.0)
        for j in 1:div(bspl.p+1,2)-1
           eig_bspl .+= (bspline(bspl.p, j-div(bspl.p+1,2), 0.0)
              .* 2 .* cos.(j * modes))
        end
    
        new( bspl.p, mesh, modes, eig_bspl, eigalpha)

    end

end

"""
    interpolate( f, bspl, alpha)

Compute the interpolating spline of degree p of odd
degree of a 1D function f on a periodic uniform mesh, at
all points x-alpha. f type is Vector{Float64}.

Derived from a Python program written by 
Eric Sonnendrucker (Max-Planck-Institut fur Plasmaphysik - Garching)

"""
function interpolate( f     :: Vector{Float64}, 
                      bspl  :: BsplinePeriodic, 
                      alpha :: Float64)

   p     = bspl.p
   delta = bspl.mesh.step
   ishift = floor(- alpha / delta)
   beta = - ishift - alpha / delta
   fill!(bspl.eigalpha, 0.0)
   for j in -div(p-1,2):div(p+1,2)
      bspl.eigalpha .+= (bspline(p, j-div(p+1,2), beta)
         .* exp.((ishift + j) * 1im .* bspl.modes))
   end

   real(ifft(fft(f) .* bspl.eigalpha ./ bspl.eig_bspl))

end

"""
    interpolate!( f, bspl, alpha)

Compute the interpolating spline of degree p of odd
degree of a 1D function f on a periodic uniform mesh, at
all points x-alpha

Array representing f is complex and modified inplace.

"""
function interpolate!( f     :: Vector{Complex{Float64}}, 
                       bspl  :: BsplinePeriodic, 
                       alpha :: Float64)

   p = bspl.p
   n = bspl.mesh.length
   delta = bspl.mesh.step
   ishift = floor(- alpha / delta)
   beta = - ishift - alpha / delta
   fill!(bspl.eigalpha,0.0)
   for j in -div(p-1,2):div(p+1,2)
      bspl.eigalpha .+= (bspline(p, j-div(p+1,2), beta)
         .* exp.((ishift + j) * 1im .* bspl.modes))
   end

   fft!(f)
   f .* bspl.eigalpha ./ bspl.eig_bspl
   ifft!(f)

end


