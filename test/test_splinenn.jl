import LinearAlgebra: transpose
import Splittings: Geometry
import Splittings: SplineNN, interpolate!

@testset "Natural Spline" begin

    geom = Geometry( 101, 101, -5.0, -5.0, 0.1, 0.1 )
    x1, x2 = geom.x1grid, geom.x2grid
    dt = 0.5
    fc = exp.(-x1.^2) .* transpose(exp.(-x2.^2))
    spline = SplinePP( geom )
    interpolate!( spline, fc, dt, dt )
    fe = exp.(-(x1.-dt).^2) .* transpose(exp.(-(x2.-dt).^2))
    error = maximum( abs.( fe .- fc ))
    println( " error  = $error ")
    @test fc ≈ fe 

    geom = Geometry( 51, 101, -5.0, -5.0, 0.2, 0.1 )
    x1, x2 = geom.x1grid, geom.x2grid
    dt = 0.5
    fc = exp.(-x1.^2) .* transpose(exp.(-x2.^2))
    spline = SplineNN( geom )
    interpolate!( spline, fc, dt, dt )
    fe = exp.(-(x1.-dt).^2) .* transpose(exp.(-(x2.-dt).^2))
    error = maximum( abs.( fe .- fc ))
    println( " error  = $error ")
    @test error < 1e-4

end
