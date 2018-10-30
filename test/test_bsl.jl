import Splittings: UniformMesh, advection!, BSpline

@testset "BSL advections" begin
  
  p = 5
  n1, n2 = 128, 128
  x1min, x1max = -10, 10
  x2min, x2max = -10, 10
  mesh1 = UniformMesh(x1min, x1max, n1; endpoint=false)
  mesh2 = UniformMesh(x2min, x2max, n2; endpoint=false)

  f  = zeros(Complex{Float64},(n1,n2))
  f .= exp.(-mesh1.points.^2) * transpose(exp.(-mesh2.points.^2))
  fᵗ = zeros(Complex{Float64},(n2,n1))

  dt =  0.5

  e = ones(Float64, n1)
  v = ones(Float64, n2)

  advection!(f,  mesh1,  v, n2, dt, BSpline(p))
  transpose!(fᵗ, f)
  advection!(fᵗ, mesh2,  e, n1, dt, BSpline(p))
  transpose!(f,  fᵗ)
  advection!(f,  mesh1, -v, n2, dt, BSpline(p))
  transpose!(fᵗ, f)
  advection!(fᵗ, mesh2, -e, n1, dt, BSpline(p))
  transpose!(f,  fᵗ)

  f0 =  exp.(-mesh1.points.^2) * transpose(exp.(-mesh2.points.^2))
  println( maximum( abs.(real(f) -f0)))

  @test real(f) ≈ f0 atol=1e-6

end
