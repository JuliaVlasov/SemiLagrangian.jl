import SemiLagrangian: LagrangePeriodicAdvection, Lagrange

@testset "Lagrange periodic advections" begin
  
  p = 5

  n1, n2 = 200, 200

  x1min, x1max = -10, 10
  x2min, x2max = -10, 10

  mesh1 = UniformMesh(x1min, x1max, n1; endpoint=false)
  mesh2 = UniformMesh(x2min, x2max, n2; endpoint=false)

  advection1! = LagrangePeriodicAdvection( mesh1, Lagrange(p) )
  advection2! = LagrangePeriodicAdvection( mesh2, Lagrange(p) )

  f  = zeros(Float64,(n1,n2))
  f .= exp.(-mesh1.points.^2) .* transpose(exp.(-mesh2.points.^2))
  f0 = copy(f)
  fᵗ = zeros(Float64,(n2,n1))

  dt = 0.1

  v1 = ones(Float64, n1)
  v2 = ones(Float64, n2)

  advection1!( f, v2, dt)

  transpose!(fᵗ, f)

  advection2!(fᵗ, v1, dt)

  transpose!(f,  fᵗ)

  advection1!(f,  -v2,  dt)

  transpose!(fᵗ, f)

  advection2!(fᵗ, -v1,  dt)

  transpose!(f,  fᵗ)

  @show maximum( abs.(f .- f0))

  @test f ≈ f0 atol=1e-6

end
