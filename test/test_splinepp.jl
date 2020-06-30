import LinearAlgebra: transpose
# import SemiLagrangian: SplinePP, interpolate!


function test_geometry(geom::Geometry{T}) where{T<:AbstractFloat}
    x1, x2 = geom.x1grid, geom.x2grid
    dt = one(T)/2
    fc = exp.(-x1.^2) .* transpose(exp.(-x2.^2))
    spline = SplinePP( geom )
    interpolate!( spline, fc, dt, dt )
    fe = exp.(-(x1.-dt).^2) .* transpose(exp.(-(x2.-dt).^2))
    error = maximum( abs.( fe .- fc ))
    return error, fc, fe
end
@testset "Spline periodic periodic" begin

    error, _1, _2 = test_geometry(Geometry( 51, 101, -5.0, -5.0, 0.2, 0.1 ))
    println( " error  = $error ")
    @test error < 1e-4
    
    error, _1, _2 = test_geometry(
    Geometry( 51, 101, big"-5.0", big"-5.0", big"0.2", big"0.1" )
)
    println( " error  = $error ")
    @test error < 1e-4

    error, fc, fe = test_geometry(Geometry( 101, 101, -5.0, -5.0, 0.1, 0.1 ))
    println( " error  = $error ")
    @test fc ≈ fe 

    error, _1, _2 = test_geometry(Geometry( 64, 128, (-2π,2π), (-2π,2π),:perxy ))
    println( " error  = $error ")
    @test error < 1e-4

end
