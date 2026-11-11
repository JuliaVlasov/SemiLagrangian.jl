using Test

@testitem "CubicSpline interpolation" begin
    nx = 100
    alpha = 0.2
    u = Float64[cos(2π * (i - 1) / nx) for i in 1:nx]
    u_out = zeros(nx)
    expected = Float64[cos(2π * ((i - 1) + alpha) / nx) for i in 1:nx]

    interpolant = CubicSpline(nx)
    interpolate!(u_out, u, 0.0, interpolant)
    @test maximum(abs.(u_out - u)) < 1.0e-14

    interpolate!(u_out, u, alpha, interpolant)
    @test maximum(abs.(u_out - expected)) < 1.0e-7

end
