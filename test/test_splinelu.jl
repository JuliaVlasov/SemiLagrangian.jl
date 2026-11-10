using LinearAlgebra
using Random
using SemiLagrangian:
    PeriodicBC,
    get_kl_ku,
    getbspline,
    LuSpline,
    BSplineLU,
    solve,
    preprocess,
    decal,
    get_order,
    split_shift,
    stencil_weights,
    interpolate!
import SemiLagrangian: solve

Random.seed!(5431221)
# create a band or circular matrix from a vector of non-zero data
function topl(n, t, iscirc = true)
    res = zeros(Rational{BigInt}, n, n)
    kl, ku = get_kl_ku(size(t, 1))
    for i in 1:n
        for (j, v) in enumerate(t)
            ind = i + j - kl - 1
            if 1 <= ind <= n
                res[i, ind] = v
            elseif iscirc
                if ind < 1
                    ind += n
                else
                    ind -= n
                end
                res[i, ind] = v
            end
        end
    end
    return res
end
# inplace lu decomposition
function decLU(A::Matrix{T}) where {T}
    n = size(A, 1)
    for k in 1:n
        pivot = A[k, k]
        for i in (k + 1):n
            A[i, k] /= pivot
        end
        for i in (k + 1):n, j in (k + 1):n
            # if A[i,k] != 0 && A[k,j] != 0
            #     println("trace avant k=$k i=$i j=$j A[i,j]=$(A[i,j]) A[i,k]=$(A[i,k]) A[k,j]=$(A[k,j])")
            # end
            A[i, j] -= A[i, k] * A[k, j]
            # if A[i,k] != 0 && A[k,j] != 0
            #     println("trace après k=$k i=$i j=$j A[i,j]=$(A[i,j])")
            # end
        end
    end
    return A
end
function getLU(A::Matrix{T}) where {T}
    n = size(A, 1)
    L = zeros(T, n, n)
    U = zeros(T, n, n)
    for i in 1:n
        L[i, i] = 1
        for j in 1:(i - 1)
            L[i, j] = A[i, j]
        end
        for j in i:n
            U[i, j] = A[i, j]
        end
    end
    return L, U
end

# function testTools()

function solve(A::Matrix{T}, Y::Vector{T}) where {T}
    L, U = getLU(A)
    n = size(A, 1)
    Y1 = zeros(T, n)
    for i in 1:n
        Y1[i] = Y[i] - sum(Y1[1:(i - 1)] .* A[i, 1:(i - 1)])
    end
    #   @assert Y1 == (L^(-1))*Y "Erreur 1"
    #   @assert isapprox(Y1, (L^(-1))*Y, atol=1e-60) "Erreur 1"
    X = zeros(T, n)
    for i in n:-1:1
        X[i] = (Y1[i] - sum(X[(i + 1):n] .* A[i, (i + 1):n])) / A[i, i]
    end
    #    @assert X == (U^(-1))*Y1 "Erreur 2"
    #    @assert isapprox(X,(U^(-1))*Y1,atol=1e-60) "Erreur 2"
    #    @assert X == ((L*U)^(-1))*Y "Erreur3"
    #    @assert isapprox(X, ((L*U)^(-1))*Y, atol=1e-60) "Erreur3"
    return X, Y1
end

function test_decLU(n)
    dividende = 10000
    A = big.((rand(Int, n, n) .% dividende) .// dividende)
    B = copy(A)
    decLU(A)
    L, U = getLU(A)
    @test L * U == B
    b = big.((rand(Int, n) .% dividende) .// dividende)
    x, _ = solve(A, b)
    return @test B * x == b
end

function test_solve(T, n, order)
    A = topl(n, getbspline(order, 0).(1:order))
    bsp = BSplineLU(order, n, Float64)
    b = rand(n)
    #    y = zeros(n)
    y = preprocess(bsp, b)
    diff = norm(A * y - b)
    @show diff
    return diff
end

function test_splu(n, order, iscirc, isLU; type = Rational{BigInt}, tol = NaN, perf = false)
    t = convert.(type, getbspline(big(order), 0).(1:order))
    A = convert.(type, topl(n, t, iscirc))
    Acp = copy(A)
    ku = div(order, 2)
    kl = order - 1 - ku
    if isLU
        perf && @time decLU(A)
        !perf && decLU(A)
    end
    spA = LuSpline(A, ku, kl; iscirc = iscirc, isLU = isLU)
    perf && @time spB = LuSpline(n, t; iscirc = iscirc, isLU = isLU)
    !perf && (spB = LuSpline(n, t; iscirc = iscirc, isLU = isLU))
    if isnan(tol)
        @test spA == spB
    end

    dividende = 10000
    b = big.((rand(Int, n) .% dividende) .// dividende)
    b = convert.(type, b)
    return if isLU
        x, y = solve(A, b)
        x2, y2 = solve(spB, b)

        if isnan(tol)
            #            @show norm(x-x2,Inf)
            @test x == x2
        else
            println("norm=$(norm(x - x2))")
            res = Acp * x
            println("normres=$(norm(res - b))")
            @test isapprox(res, b, atol = tol)
            @test isapprox(x, x2, atol = tol)
        end
    end
end

function test_perf(n, order, iscirc)
    t = getbspline(big(order), 0).(1:order)
    A = topl(n, t, iscirc)
    ku = div(order, 2)
    kl = order - 1 - ku
    @time decLU(A)

    spA = LuSpline(A, ku, kl; iscirc = iscirc, isLU = true)
    @time spB = LuSpline(n, t; iscirc = iscirc, isLU = true)

    dividende = 10000
    b = big.((rand(Int, n) .% dividende) .// dividende)
    x, y = solve(A, b)
    x2, y2 = solve(spB, b)
    @test y == y2
    return @test x == x2
end

function test_interface()
    getpar(bsp::BSplineLU{T, tedge, order}) where {T, tedge, order} = (T, tedge, order)
    for i in 3:2:29
        bsp = BSplineLU(i, 104, big"0.")
        (_1, _2, order) = getpar(bsp)
        @test i == get_order(bsp)
        @test i == order
        @test getpar(bsp) == (BigFloat, PeriodicBC, order)
    end
    return
end

@testset "test decLU that is a tool for test" begin
    test_decLU(30)
end

@testset "test BSplineLU preprocessing" begin
    @test test_solve(Float64, 128, 9) < 1.0e-10
end

@testset "even-order BSplineLU" begin
    T = Rational{BigInt}
    for order in (2, 4, 6), n in (16, 17)
        spline = getbspline(order, 0)
        interpolant = BSplineLU(order, n, T)
        coefficients = T[mod(7i + 3, 11) - 5 for i in 1:n]
        half = T(1 // 2)
        radius = order ÷ 2
        evaluate(i, shift) = sum(
            coefficients[mod1(i + j, n)] * spline(radius + half + shift - j)
                for j in (-radius - 2):(radius + 2)
        )
        input = [evaluate(i, zero(T)) for i in 1:n]
        @test preprocess(interpolant, input) == coefficients

        for shift in T.((-5 // 4, -1 // 2, -1 // 4, 0, 1 // 4, 1 // 2, 3 // 4, 5 // 4))
            expected = [evaluate(i, shift) for i in 1:n]
            output = similar(input)
            interpolate!(output, input, shift, interpolant)
            @test output == expected

            integer, fractional = split_shift(interpolant, shift)
            weights = stencil_weights(interpolant, fractional)
            interpolate!(output, input, integer, weights, interpolant)
            @test output == expected
        end

        shifts = [T(mod(i, 7) - 3) / 4 for i in 1:n]
        expected = [evaluate(i, shifts[i]) for i in 1:n]
        output = similar(input)
        interpolate!(output, copy(input), ind -> (shifts[ind[1]],), [interpolant])
        @test output == expected
    end

    for order in (2, 4, 6), n in (64, 65)
        input = cos.(2pi .* (0:n-1) ./ n)
        interpolant = BSplineLU(order, n)
        output = similar(input)
        for shift in (-1.25, -0.5, 0.0, 0.5, 0.75, 1.25)
            interpolate!(output, input, shift, interpolant)
            expected = cos.(2pi .* ((0:n-1) .+ shift) ./ n)
            @test maximum(abs.(output .- expected)) < 1.0e-4
        end
        alternating = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
        interpolate!(output, alternating, 0.0, interpolant)
        @test output ≈ alternating
    end

    setprecision(128) do
        n = 64
        input = [cos(2 * big(pi) * i / n) for i in 0:(n - 1)]
        interpolant = BSplineLU(4, n, BigFloat)
        output = similar(input)
        shift = big"0.75"
        interpolate!(output, input, shift, interpolant)
        expected = [cos(2 * big(pi) * (i + shift) / n) for i in 0:(n - 1)]
        @test maximum(abs.(output .- expected)) < big"1e-7"
    end

    n, m = 8, 9
    splines = (getbspline(2, 0), getbspline(4, 0))
    interpolants = [BSplineLU(2, n, T), BSplineLU(4, m, T)]
    coefficients = T[mod(3i + 5j, 13) - 6 for i in 1:n, j in 1:m]
    half = T(1 // 2)
    evaluate2d(i, j, sx, sy) = sum(
        coefficients[mod1(i + k, n), mod1(j + l, m)] *
            splines[1](1 + half + sx - k) * splines[2](2 + half + sy - l)
            for k in -3:3, l in -4:4
    )
    input = [evaluate2d(i, j, zero(T), zero(T)) for i in 1:n, j in 1:m]
    shifts(i, j) = (T(mod(i, 5) - 2) / 4, T(mod(j, 5) - 2) / 4)
    expected = [evaluate2d(i, j, shifts(i, j)...) for i in 1:n, j in 1:m]
    output = similar(input)
    interpolate!(output, copy(input), ind -> shifts(ind[1], ind[2]), interpolants)
    @test output == expected
end

@testset "test LuSpline" begin
    Random.seed!(1899221)
    test_splu(30, 9, true, false)
    test_splu(30, 10, true, false)
    test_splu(30, 9, false, false)
    test_splu(30, 10, false, false)
    test_splu(30, 9, false, true)
    test_splu(30, 10, false, true)
    test_splu(30, 9, true, true)
    test_splu(31, 10, true, true)
    test_splu(30, 9, true, true; type = BigFloat, tol = 1.0e-70)
    test_splu(31, 10, true, true; type = BigFloat, tol = 1.0e-70)
    test_splu(30, 11, false, true; type = BigFloat, tol = 1.0e-70)
    test_splu(31, 10, false, true; type = BigFloat, tol = 1.0e-70)
    test_splu(30, 3, true, false)
    test_splu(30, 4, true, false)
    test_splu(30, 3, false, false)
    test_splu(30, 4, false, false)
    test_splu(30, 3, false, true)
    test_splu(30, 4, false, true)
    test_splu(30, 3, true, true)
    test_splu(31, 4, true, true)
    test_splu(30, 3, true, true; type = BigFloat, tol = 1.0e-70)
    test_splu(31, 4, true, true; type = BigFloat, tol = 1.0e-70)
    test_splu(30, 5, false, true; type = BigFloat, tol = 1.0e-70)
    test_splu(31, 4, false, true; type = BigFloat, tol = 1.0e-70)
    test_interface()
end

@testset "test perf" begin
    test_splu(1000, 9, true, true; type = Float64, tol = 1.0e-12, perf = true)
end
