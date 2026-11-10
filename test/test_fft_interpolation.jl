using Test

@testitem "Spectral interpolation" begin
    N = 100
    alpha = 0.2
    u = Float64[cos(2π * (i - 1) / N) for i in 1:N]
    u_out = zeros(N)
    expected = Float64[cos(2π * ((i - 1) + alpha) / N) for i in 1:N]

    interpolant = Spectral(N)

    interpolate!(u_out, u, 0.0, interpolant)
    @test maximum(abs.(u_out - u)) < 1.0e-14

    interpolate!(u_out, u, alpha, interpolant)
    @test maximum(abs.(u_out - expected)) < 1.0e-14

end
