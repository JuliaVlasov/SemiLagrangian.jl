import VlasovBase:UniformMesh

"""

   exact(tf, mesh1, mesh2)

   Julia function to compute exact solution "

```math
\frac{d f}{dt} +  (y \frac{df}{delta1} - x \frac{df}{delta2}) = 0 
```

"""
function exact(tf, mesh1::UniformMesh, mesh2::UniformMesh)

    f = zeros(Float64,(mesh1.length,mesh2.length))
    for (i, x) in enumerate(mesh1.points), (j, y) in enumerate(mesh2.points)
        xn = cos(tf) * x - sin(tf) * y
        yn = sin(tf) * x + cos(tf) * y
        f[i,j] = exp(-(xn-1)*(xn-1)/0.1)*exp(-(yn-1)*(yn-1)/0.1)
    end

    f

end

" Function to compute error "
function error1(f, f_exact)
    maximum(abs.(f .- f_exact))
end

import FFTW: FFTWPlan, plan_fft, fft!, ifft!

struct Fourier

    px :: FFTWPlan
    py :: FFTWPlan
    kx :: Vector{Float64}
    ky :: Vector{Float64}

    function Fourier( mesh1 :: UniformMesh,
		      mesh2 :: UniformMesh,
		      f  :: Array{Complex{Float64},2},
                      fᵗ :: Array{Complex{Float64},2})

        px = plan_fft(f,  1)
        py = plan_fft(fᵗ, 1)

        n1   = mesh1.length
        x1min = mesh1.start
        x1max = mesh1.stop
        kx = 2π/(x1max-x1min)*[0:n1÷2-1;n1÷2-n1:-1]

        n2   = mesh2.length
        ymin = mesh2.start
        ymax = mesh2.stop
        ky = 2π/(ymax-ymin)*[0:n2÷2-1;n2÷2-n2:-1]

	new( px, py, kx, ky)

    end
end

function advection_y!( f ::Array{Complex{Float64},2}, 
		       fᵗ::Array{Complex{Float64},2}, 
		       f̂ᵗ::Array{Complex{Float64},2}, 
		       mesh1::UniformMesh, 
		       dt::Float64, 
		       fourier::Fourier)

    exky = exp.( 1im * dt * fourier.ky .* transpose(mesh1.points))

    transpose!(fᵗ,f)
    mul!(f̂ᵗ,  fourier.py, fᵗ)
    f̂ᵗ .= f̂ᵗ .* exky
    ldiv!(fᵗ, fourier.py, f̂ᵗ)
    transpose!(f,fᵗ)

end

function advection_x!( f ::Array{Complex{Float64},2}, 
		       f̂ ::Array{Complex{Float64},2}, 
		       mesh2::UniformMesh, 
		       dt::Float64, 
		       fourier::Fourier)

    ekxy = exp.(-1im * dt * fourier.kx .* transpose(mesh2.points))

    mul!(f̂,  fourier.px, f)
    f̂ .= f̂ .* ekxy 
    ldiv!(f, fourier.px, f̂)

end

function rotation_2d_fft(tf, nt, mesh1::UniformMesh, mesh2::UniformMesh)

    dt = tf/nt

    n1 = mesh1.length
    x1min, x1max = mesh1.start, mesh1.stop
    delta1 = mesh1.step

    n2 = mesh2.length
    ymin, ymax = mesh2.start, mesh2.stop
    delta2 = mesh2.step

    f  = zeros(Complex{Float64},(n1,n2))
    f̂  = similar(f)
    fᵗ = zeros(Complex{Float64},(n2,n1))
    f̂ᵗ = similar(fᵗ)
    
    fourier = Fourier( mesh1, mesh2, f, fᵗ)
    
    f .= exact(0.0, mesh1, mesh2)
    
    for n=1:nt

        advection_y!( f, fᵗ, f̂ᵗ, mesh1, tan(dt/2), fourier)
	advection_x!( f, f̂, mesh2, sin(dt), fourier)
        advection_y!( f, fᵗ, f̂ᵗ, mesh1, tan(dt/2), fourier)

    end
    real(f)
end

tf, nt = 200π, 1000

mesh1 = UniformMesh(-π, π, 128; endpoint=false)
mesh2 = UniformMesh(-π, π, 256; endpoint=false)

fc = rotation_2d_fft(tf, nt, mesh1, mesh2)
fe = exact(tf, mesh1, mesh2)

@testset "Rotation test with Fourier advections " begin

println(error1(fc, fe))
@test rotation_2d_fft(tf, nt, mesh1, mesh2) ≈ fe atol = 1e-10

end
