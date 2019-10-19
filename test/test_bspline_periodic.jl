import SemiLagrangian: BsplinePeriodic, interpolate!

@testset "Bspline periodic advections" begin
  
  p = 5

  n1, n2 = 128, 128

  x1min, x1max = -10, 10
  x2min, x2max = -10, 10

  mesh1 = UniformMesh(x1min, x1max, n1; endpoint=false)
  mesh2 = UniformMesh(x2min, x2max, n2; endpoint=false)

  bspl1 = BsplinePeriodic( Bspline(p), mesh1 )
  bspl2 = BsplinePeriodic( Bspline(p), mesh2 )

  f  = zeros(ComplexF64,(n1,n2))
  f .= exp.(-mesh1.points.^2) .* transpose(exp.(-mesh2.points.^2))
  f0 = copy(f)
  fᵗ = zeros(ComplexF64,(n2,n1))

  dt = 0.5

  v1 = ones(Float64, n1)
  v2 = ones(Float64, n2)

  for j = 1:n2
      alpha = dt * v2[j]
      interpolate!(f[:,j],  bspl1,  alpha)
  end

  transpose!(fᵗ, f)

  for j = 1:n1
      alpha = dt * v1[j]
      interpolate!(fᵗ[:,j],  bspl2,  alpha)
  end

  transpose!(f,  fᵗ)

  for j = 1:n2
      alpha = - dt * v2[j]
      interpolate!(f[:,j],  bspl1,  alpha)
  end

  transpose!(fᵗ, f)

  for j = 1:n1
      alpha = - dt * v1[j]
      interpolate!(fᵗ[:,j],  bspl2,  alpha)
  end

  transpose!(f,  fᵗ)

  @show maximum( abs.(real(f) .- f0))

  @test real(f) ≈ f0 atol=1e-6

end
