using FFTW, LinearAlgebra
import VlasovBase:UniformMesh

export BSpline

struct BSpline

    p :: Int64

    function BSpline( p )
	
	@assert (p & 1 == 1)
	new( p )

    end
end


"""
    bspline(p, j, x)

Return the value at x in [0,1[ of the B-spline with
integer nodes of degree p with support starting at j.
Implemented recursively using the de Boor's recursion formula

Derived from a Python program written by 
Eric Sonnendrucker (Max-Planck-Institut fur Plasmaphysik - Garching))

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


"""
    interpolate( p, f, delta, alpha)

Compute the interpolating spline of degree p of odd
degree of a 1D function f on a periodic uniform mesh, at
all points x-alpha. f type is Vector{Float64}.

Derived from a Python program written by 
Eric Sonnendrucker (Max-Planck-Institut fur Plasmaphysik - Garching)

"""
function interpolate(p     :: Int, 
                     f     :: Vector{Float64}, 
                     delta :: Float64, 
                     alpha :: Float64)

   n = size(f)[1]
   modes = 2π * (0:n-1) / n
   eig_bspl = zeros(Complex{Float64},n)
   eig_bspl .= bspline(p, -div(p+1,2), 0.0)
   for j in 1:div(p+1,2)-1
      eig_bspl .+= (bspline(p, j-div(p+1,2), 0.0)
         * 2 * cos.(j * modes))
   end
   ishift = floor(- alpha / delta)
   beta = - ishift - alpha / delta
   eigalpha = zeros(Complex{Float64},n)
   for j in -div(p-1,2):div(p+1,2)
      eigalpha .+= (bspline(p, j-div(p+1,2), beta)
         .* exp.((ishift + j) * 1im .* modes))
   end

   real(ifft(fft(f) .* eigalpha ./ eig_bspl))

end

"""
    interpolate( p, f, delta, alpha)

Compute the interpolating spline of degree p of odd
degree of a 1D function f on a periodic uniform mesh, at
all points x-alpha
input f is the Fourier transform and complex
you have to in2erse transform after return

"""
function interpolate(p     :: Int, 
                     f     :: Vector{Complex{Float64}}, 
                     delta :: Float64, 
                     alpha :: Float64)

   n = size(f)[1]
   modes = 2 * pi * (0:n-1) / n
   eig_bspl = zeros(Complex{Float64},n)
   eig_bspl .= bspline(p, -div(p+1,2), 0.0)
   for j in 1:div(p+1,2)-1
      eig_bspl .+= (bspline(p, j-div(p+1,2), 0.0)
         * 2 * cos.(j * modes))
   end
   ishift = floor(- alpha / delta)
   beta = - ishift - alpha / delta
   eigalpha = zeros(Complex{Float64},n)
   for j in -div(p-1,2):div(p+1,2)
      eigalpha .+= (bspline(p, j-div(p+1,2), beta)
         .* exp.((ishift + j) * 1im .* modes))
   end

   f .* eigalpha ./ eig_bspl

end


"""
    advection!( mesh, f, v, dt, interp, axis)

Advection of a 2d function `f` discretized on a 2d `mesh`
along the input axis at velocity `v`

"""
function advection!(f      :: Array{Float64,2}, 
                    mesh1  :: UniformMesh,
                    mesh2  :: UniformMesh,
                    v      :: Vector{Float64}, 
                    dt     :: Float64,
                    interp :: BSpline,
                    axis   :: Int64 )

    @assert ( axis == 1 || axis == 2 )

    p = interp.p

    if (axis == 1)
        @simd for j in eachindex(v)
            alpha = v[j] * dt
            @inbounds f[:,j] .= interpolate(p, f[:,j], mesh1.step, alpha)
        end
    else
        @simd for i in eachindex(v)
            alpha = v[i] * dt
            @inbounds f[i,:] .= interpolate(p, f[i,:], mesh2.step, alpha)
        end
    end

end

"""
    advection!(f, mesh, v, n2, dt, interp)

Advection of a 2d function `f` along its first dimension with
velocity `v`. Since the fft are computed inplace, the function 
must be represented by a Array{Complex{Float64},2}.

"""
function advection!(f::Array{Complex{Float64},2}, 
                    mesh::UniformMesh, 
                    v::Vector{Float64}, 
                    n2::Int, 
                    dt::Float64,
                    interp::BSpline)

    p  = interp.p
    n1 = mesh.length
    delta1 = mesh.step
    modes = [2π * i / n1 for i in 0:n1-1]
    # compute eigen2alues of degree p b-spline matrix
    eig_bspl = zeros(Float64, n1)
    eig_bspl .= bspline(p, -div(p+1,2), 0.0)
    for i in 1:div(p+1,2)-1
       eig_bspl .+= bspline(p, i - div(p+1,2), 0.0) * 2 .* cos.(i * modes)
    end
    eigalpha = zeros(Complex{Float64}, n1)

    fft!(f,1)

    @simd for j in 1:n2
       alpha = dt * v[j] / delta1

       # compute eigen2alues of cubic splines evaluated at displaced points
       ishift = floor(-alpha)
       beta   = -ishift - alpha
       fill!(eigalpha,0.0im)
       for i in -div(p-1,2):div(p+1,2)
          eigalpha .+= (bspline(p, i-div(p+1,2), beta)
                         .* exp.((ishift+i) * 1im .* modes))
       end

       # compute interpolating spline using fft and properties of circulant matrices
       @inbounds f[:,j] .*= eigalpha ./ eig_bspl

    end

    ifft!(f,1)

end
