using Test

@testitem "CubicSpline interpolation" begin
    nx = 100
    alpha = 0.2
    u = Float64[cos(2π * (i - 1) / nx) for i = 1:nx]
    u_out = zeros(nx)
    expected = Float64[cos(2π * ((i - 1) + alpha) / nx) for i = 1:nx]

    interpolant = CubicSpline(nx)
    interpolate!(u_out, interpolant, u, 0.0)
    @test maximum(abs.(u_out - u)) < 1e-14

    interpolate!(u_out, interpolant, u, alpha)
    @test maximum(abs.(u_out - expected)) < 1e-7

end
