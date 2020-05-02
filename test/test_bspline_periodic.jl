using SemiLagrangian

@testset "Bspline periodic advections" begin
  
    for p in (3, 5, 7)

        n1, n2 = 128, 128

        x1min, x1max = -10, 10
        x2min, x2max = -10, 10

        mesh1 = UniformMesh(x1min, x1max, n1; endpoint=false)
        mesh2 = UniformMesh(x2min, x2max, n2; endpoint=false)

        advection1! = BsplinePeriodicAdvection( mesh1, Bspline(p) )
        advection2! = BsplinePeriodicAdvection( mesh2, Bspline(p) )

        f  = zeros(ComplexF64,(n1,n2))
        f .= exp.(-mesh1.points.^2) .* transpose(exp.(-mesh2.points.^2))
        f0 = copy(f)
        fᵗ = zeros(ComplexF64,(n2,n1))

        dt = 0.5

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

        @show maximum( abs.(real(f) .- f0))

        @test real(f) ≈ real(f0) atol=1/10^p

    end

end
