"""
# Test Program for Lagrange Interpolation 1D

This test program validates the Lagrange interpolation functions by:
1. Using a known smooth test function: f(x) = cos(2π*x/num_points)
2. Testing various interpolation orders (3, 4, 5, 6 points)
3. Testing different boundary condition types
4. Verifying that interpolation errors are within expected tolerances

## Test Cases
The program runs 5 test cases:
1. Order 5, no BC: Tests interior interpolation only
2. Order 3, periodic BC: Tests basic periodic boundary conditions
3. Order 5, periodic with extended output: Tests periodic with fp[n+1] = fp[1]
4. Order 4, centered periodic: Tests even-order stencil with auto interval shift
5. Order 6, centered periodic: Tests higher-order even stencil

## Expected Output
If all tests pass, the program prints "PASSED."
If any test fails, the program prints "FAILED."
"""

using PeriodicInterpolation1D

@testmodule CommonHelpers begin


    using PeriodicInterpolation1D


    """
        test_interpolation(num_points, fi, alpha, xp, order, type, tolerance) -> Bool

    Run a single interpolation test case.

    # Arguments
    - `num_points`: Number of interpolation points
    - `fi`: Known function values at grid points
    - `alpha`: Displacement parameter (offset from grid points)
    - `xp`: x values with offset (where to interpolate)
    - `order`: Interpolation order (stencil size: 3, 4, 5, 6, etc.)
    - `type`: Boundary condition type
      - 0: no BC
      - 1: periodic
      - 2: periodic with last value (fp[n+1] = fp[1])
      - 3: centered periodic (even stencils)
    - `tolerance`: Maximum allowed interpolation error

    # Returns
    - Updated test status (false if this test failed, unchanged otherwise)

    # Description
    Performs interpolation, computes maximum error compared to exact function values,
    and checks if error is within tolerance.
    """
    function test_interpolation(
        f::Function,
        num_points::Int,
        fi::Vector{Float64},
        alpha::Float64,
        xp::Vector{Float64},
        order::Int,
        tolerance::Float64,
    )


        fp = zeros(Float64, num_points)
        diff = 0.0

        pmessage = string(order)

        println("Test fixed_periodic with order ", pmessage, " .")
        interpolant = FastLagrange(order)
        interpolate!(fp, interpolant, fi, alpha)

        for i = 1:num_points
            diff = max(diff, abs(f(xp[i], num_points) - fp[i]))
        end

        println("error = ", diff)

        return diff < tolerance

    end

end


@testsnippet SharedData begin

    f(x::Float64, num_points::Int) = cos(2 * π * x / num_points)
    g(x::Float64, num_points::Int) = sin(2 * π * x / num_points)

    num_points = 100
    alpha = 0.2

    # Allocate arrays
    xi = zeros(Float64, num_points)
    fi = zeros(Float64, num_points)
    gi = zeros(Float64, num_points)
    xp = zeros(Float64, num_points)

    # Data initialization
    xmin = 0.0
    xmax = Float64(num_points - 1)
    l = xmax - xmin

    for i = 1:(num_points)
        xi[i] = Float64(i - 1)
        fi[i] = f(xi[i], num_points)
        gi[i] = g(xi[i], num_points)
        xp[i] = xi[i] + alpha
    end

end

@testitem "periodic order 3" tags = [:Lagrange] setup = [CommonHelpers, SharedData] begin
    @test CommonHelpers.test_interpolation(f, num_points, fi, alpha, xp, 3, 8.0e-6)
    @test CommonHelpers.test_interpolation(g, num_points, gi, alpha, xp, 3, 1.0e-2)
end

@testitem "periodic order 5" tags = [:Lagrange] setup = [CommonHelpers, SharedData] begin
    @test CommonHelpers.test_interpolation(f, num_points, fi, alpha, xp, 5, 1.0e-8)
    @test CommonHelpers.test_interpolation(g, num_points, gi, alpha, xp, 5, 1.0e-8)
end

@testitem "periodic order 7" tags = [:Lagrange] setup = [CommonHelpers, SharedData] begin
    @test CommonHelpers.test_interpolation(f, num_points, fi, alpha, xp, 7, 1.0e-8)
    @test CommonHelpers.test_interpolation(g, num_points, gi, alpha, xp, 7, 1.0e-8)
end

@testitem "periodic order 9" tags = [:Lagrange] setup = [CommonHelpers, SharedData] begin
    @test CommonHelpers.test_interpolation(f, num_points, fi, alpha, xp, 9, 1.0e-14)
    @test CommonHelpers.test_interpolation(g, num_points, gi, alpha, xp, 9, 1.0e-14)
end

@testitem "periodic order 11" tags = [:Lagrange] setup = [CommonHelpers, SharedData] begin
    @test CommonHelpers.test_interpolation(f, num_points, fi, alpha, xp, 11, 1.0e-14)
    @test CommonHelpers.test_interpolation(g, num_points, gi, alpha, xp, 11, 1.0e-14)
end
