using Test

@testitem "Spectral interpolation" begin

    N = 100
    alpha = 0.2
    u = Float64[cos(2π * (i - 1) / N) for i = 1:N]
    u_out = zeros(N)
    expected = Float64[cos(2π * ((i - 1) + alpha) / N) for i = 1:N]

    interpolant = Spectral(N)

    interpolate!(u_out, interpolant, u, 0.0)
    @test maximum(abs.(u_out - u)) < 1e-14

    interpolate!(u_out, interpolant, u, alpha)
    @test maximum(abs.(u_out - expected)) < 1e-14

end
