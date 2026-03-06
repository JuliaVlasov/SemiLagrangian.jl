using Test

@testitem "Lagrange interpolation" begin
    nx = 100
    alpha = 0.2
    u = Float64[cos(2π * (i - 1) / nx) for i = 1:nx]
    u_out = zeros(nx)
    expected = Float64[cos(2π * ((i - 1) + alpha) / nx) for i = 1:nx]

    for order in [3, 5, 7, 9]

        interpolant = Lagrange(nx, order)
        interpolate!(u_out, interpolant, u, 0.0)
        @test maximum(abs.(u_out - u)) < 1e-14

        interpolate!(u_out, interpolant, u, alpha)
        # FFT-based Lagrange has excellent accuracy for smooth functions
        # Tolerance decreases with order
        tol = 10.0^(-(2 + order/2))
        @test maximum(abs.(u_out - expected)) < tol

    end
end
