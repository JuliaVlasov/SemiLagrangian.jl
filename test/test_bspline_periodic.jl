# using SemiLagrangian

include("../src/bspline_periodic.jl")
using Test

@testset "Bspline periodic advections" begin
  
    for p in (3, 5, 7, 9, 11, 13)

        n1, n2 = 128, 128

        x1min, x1max = -10.0, 10.0
        x2min, x2max = -10.0, 10.0

        mesh1 = UniformMesh(x1min, x1max, n1; endpoint=false)
        mesh2 = UniformMesh(x2min, x2max, n2; endpoint=false)

        adv1 = BsplinePeriodicAdvection( mesh1, Bspline(p) )
        adv2 = BsplinePeriodicAdvection( mesh2, Bspline(p) )

        f  = zeros(ComplexF64,(n1,n2))
        f .= exp.(-mesh1.points.^2) .* transpose(exp.(-mesh2.points.^2))
        f0 = copy(f)
        fᵗ = zeros(ComplexF64,(n2,n1))

        dt = 0.5

        v1 = ones(Float64, n1)
        v2 = ones(Float64, n2)

        advection!(adv1, f, v2, dt)

        transpose!(fᵗ, f)

        advection!(adv2, fᵗ, v1, dt)

        transpose!(f,  fᵗ)

        advection!(adv1, f,  -v2,  dt)

        transpose!(fᵗ, f)

        advection!(adv2, fᵗ, -v1,  dt)

        transpose!(f,  fᵗ)

        @show maximum( abs.(real(f) .- f0))

        @test real(f) ≈ real(f0) atol=1/10^p

    end

end

# @testset "Bspline periodic advections with big correctbug" begin
  
# # for p in (3, 5, 7, 9, 11, 13)
#     for p in (3)
#         println("trace test big p=$p")
#         n1, n2 = 128, 128

#         x1min, x1max = -big"10.0", big"10.0"
#         x2min, x2max = -big"10.0", big"10.0"
#         x1minf, x1maxf = -10.0, 10.0
#         x2minf, x2maxf = -10.0, 10.0

#         mesh1 = UniformMesh(x1min, x1max, n1; endpoint=false)
#         mesh2 = UniformMesh(x2min, x2max, n2; endpoint=false)
#         mesh1f = UniformMesh(x1minf, x1maxf, n1; endpoint=false)
#         mesh2f = UniformMesh(x2minf, x2maxf, n2; endpoint=false)

#         adv1 = BsplinePeriodicAdvection( mesh1, Bspline(p) )
#         adv2 = BsplinePeriodicAdvection( mesh2, Bspline(p) )
#         adv1f = BsplinePeriodicAdvection( mesh1f, Bspline(p) )
#         adv2f = BsplinePeriodicAdvection( mesh2f, Bspline(p) )

#         f  = zeros(Complex{BigFloat},(n1,n2))
#         f .= exp.(-mesh1.points.^2) .* transpose(exp.(-mesh2.points.^2))
#         f0 = copy(f)
#         fᵗ = zeros(Complex{BigFloat},(n2,n1))
#         ff  = zeros(ComplexF64,(n1,n2))
#         ff .= exp.(-mesh1f.points.^2) .* transpose(exp.(-mesh2f.points.^2))
#         ff0 = copy(ff)
#         ffᵗ = zeros(ComplexF64,(n2,n1))

#         println("normdiff f0 ff0=$(norm(f0-ff0))")
#         println("1normdiff f ff=$(norm(f-ff))")

#         dt = big"0.5"
#         dtf = 0.5

#         v1 = ones(BigFloat, n1)
#         v2 = ones(BigFloat, n2)
#         v1f = ones(Float64, n1)
#         v2f = ones(Float64, n2)

#         advection!(adv1, f, v2, dt)
#         advection!(adv1f, ff, v2f, dtf)
#         println("2normdiff eigalpha=$(norm(adv1.eigalpha-adv1f.eigalpha))")
#         println("2normdiff eig_bspl=$(norm(adv1.eig_bspl-adv1f.eig_bspl))")
#         println("2normdiff f ff=$(norm(f-ff))")

#         transpose!(fᵗ, f)
#         transpose!(ffᵗ, ff)

#         advection!(adv2, fᵗ, v1, dt)
#         advection!(adv2f, ffᵗ, v1f, dtf)

#         transpose!(f,  fᵗ)
#         transpose!(ff,  ffᵗ)
#         println("3normdiff f ff=$(norm(f-ff))")

#         advection!(adv1, f,  -v2,  dt)
#         advection!(adv1f, ff,  -v2f,  dtf)
#         println("4normdiff f ff=$(norm(f-ff))")

#         transpose!(fᵗ, f)
#         transpose!(ffᵗ, ff)

#         advection!(adv2, fᵗ, -v1,  dt)
#         advection!(adv2f, ffᵗ, -v1f,  dtf)

#         transpose!(f,  fᵗ)
#         transpose!(ff,  ffᵗ)
#         println("5normdiff f ff=$(norm(f-ff))")

#         @show maximum( abs.(real(f) .- f0))
#         @show maximum( abs.(real(ff) .- ff0))

#         @test real(f) ≈ real(f0) atol=1/10^p

#     end

# end

@testset "Bspline periodic advections with big" begin
  
    for p in (3, 5, 7, 9, 11, 13)
        println("trace test big p=$p")
        n1, n2 = 128, 128

        x1min, x1max = -big"10.0", big"10.0"
        x2min, x2max = -big"10.0", big"10.0"

        mesh1 = UniformMesh(x1min, x1max, n1; endpoint=false)
        mesh2 = UniformMesh(x2min, x2max, n2; endpoint=false)
 
        adv1 = BsplinePeriodicAdvection( mesh1, Bspline(p) )
        adv2 = BsplinePeriodicAdvection( mesh2, Bspline(p) )

        f  = zeros(Complex{BigFloat},(n1,n2))
        f .= exp.(-mesh1.points.^2) .* transpose(exp.(-mesh2.points.^2))
        f0 = copy(f)
        fᵗ = zeros(Complex{BigFloat},(n2,n1))

        dt = big"0.5"

        v1 = ones(BigFloat, n1)
        v2 = ones(BigFloat, n2)

        advection!(adv1, f, v2, dt)

        transpose!(fᵗ, f)

        advection!(adv2, fᵗ, v1, dt)

        transpose!(f,  fᵗ)

        advection!(adv1, f,  -v2,  dt)

        transpose!(fᵗ, f)

        advection!(adv2, fᵗ, -v1,  dt)

        transpose!(f,  fᵗ)

        @show maximum( abs.(real(f) .- f0))

        @test real(f) ≈ real(f0) atol=1/10^p

    end

end

