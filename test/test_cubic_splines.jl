import Splittings:compute_interpolants, interpolate, CubicSpline

@testset "CubicSpline interpolation" begin

    function interpolation_test(n1 = 128)
    
        x1min, x1max = 0.0, 1.0
        x = collect(range(x1min,stop=x1max,length=n1))
        y = sin.(2π*x)
        x_new = zeros(Float64,n1)
        y_new = zeros(Float64,n1)
        coeffs = compute_interpolants(n1, y) 
        for (i, xi) in enumerate(x)
            x_new[i] = xi - 0.1
            x_new[i] = x1min + mod(x_new[i]-x1min,x1max-x1min) 
            y_new[i] = interpolate(coeffs, n1, x1min, x1max, x_new[i])
        end
        maximum(abs.(sin.(2π*(x.-0.1)) - y_new))
    
    end

    @test ≈(interpolation_test(), 0.0, atol=1e-7)

end


import LinearAlgebra: transpose
import Splittings: UniformMesh, advection!

@testset " BSL advections with cubic splines " begin
  
  n1, n2 = 128, 128
  x1min, x1max = -5, 10
  x2min, x2max = -5, 10
  mesh1 = UniformMesh(x1min, x1max, n1)
  mesh2 = UniformMesh(x2min, x2max, n2)

  f  = zeros(Float64,(n1,n2))
  f .= exp.(-mesh1.points.^2) * transpose(exp.(-mesh2.points.^2))
  fᵗ = zeros(Float64,(n2,n1))

  dt =  0.5

  e = ones(Float64, n1)
  v = ones(Float64, n2)

  advection!(f, mesh1,  v, dt, CubicSpline(), 1)
  advection!(f, mesh2,  e, dt, CubicSpline(), 2)
  advection!(f, mesh1, -v, dt, CubicSpline(), 1)
  advection!(f, mesh2, -e, dt, CubicSpline(), 2)

  f0 =  exp.(-mesh1.points.^2) * transpose(exp.(-mesh2.points.^2))
  println( maximum( abs.(f .- f0)))

  @test f ≈ f0 atol=1e-3

end
