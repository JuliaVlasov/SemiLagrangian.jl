

@testset "Lagrange periodic advections" begin
  
  for p in 3:2:9

     n1, n2 = 200, 200

     x1min, x1max = -big"10.0", big"10.0"
     x2min, x2max = -big"10.0", big"10.0"

     mesh1 = UniformMesh(x1min, x1max, n1; endpoint=false)
     mesh2 = UniformMesh(x2min, x2max, n2; endpoint=false)

     adv1 = Advection( mesh1, LagrangeNew(p) )
     adv2 = Advection( mesh2, LagrangeNew(p) )

     f  = zeros(typeof(x1min),(n1,n2))
     f .= exp.(-mesh1.points.^2) .* transpose(exp.(-mesh2.points.^2))
     f0 = copy(f)
     fᵗ = zeros(typeof(x1min),(n2,n1))

     dt = big"0.1"

     v1 = ones(typeof(x1min), n1)
     v2 = ones(typeof(x1min), n2)

     advection!( adv1, f, v2, dt)

     transpose!(fᵗ, f)

     advection!(adv2, fᵗ, v1, dt)

     transpose!(f,  fᵗ)

     advection!(adv1, f,  -v2,  dt)

     transpose!(fᵗ, f)

     advection!( adv2, fᵗ, -v1,  dt)

     transpose!(f,  fᵗ)

     @show maximum( abs.(f .- f0))

     @test f ≈ f0 atol=1e-6

  end

end
