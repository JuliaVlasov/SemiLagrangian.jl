include("../src/spline.jl")
include("../src/matspline.jl")
using LinearAlgebra
using Random
using Test
Random.seed!(54321)
function test_decLU(n)
    dividende=10000
    A = big.((rand(Int, n, n) .% dividende) .// dividende)
    B = copy(A)
    decLU(A)
    L, U = getLU(A)
    @test L*U == B
    b = big.((rand(Int, n) .% dividende) .// dividende)
    x, _ = sol(A, b)
    @test B*x == b
end

function test_splu(n, order, iscirc, isLU; type=Rational{BigInt}, tol=NaN, perf=false)
   t = convert.(type, getbspline(big(order),0).(1:order))
    A = convert.(type, topl(n, t, iscirc))
    Acp = copy(A)
    ku = div(order,2)
    kl = order-1-ku
    if isLU
        perf && @time decLU(A)
        !perf && decLU(A)
    end
    spA = LuSpline(A, ku, kl; iscirc=iscirc, isLU=isLU)
    perf && @time spB = LuSpline(n, t; iscirc=iscirc, isLU=isLU)
    !perf && (spB = LuSpline(n, t; iscirc=iscirc, isLU=isLU))
    if isnan(tol)
        @test spA == spB
    end

    dividende = 10000
    b = big.((rand(Int, n) .% dividende) .// dividende)
    b = convert.(type, b)
    if isLU
        x, y = sol(A, b)
        x2, y2 = sol(spB, b)
        
        if isnan(tol)
            @test x == x2
        else
            println("norm=$(norm(x-x2))")
            res=Acp*x
            println("normres=$(norm(res-b))")
            @test isapprox(res, b, atol=tol)
            @test isapprox(x, x2, atol=tol)
        end
    end
end

function test_perf(n ,order, iscirc)
    t = getbspline(big(order),0).(1:order)
    A = topl(n, t, iscirc)
    ku = div(order,2)
    kl = order-1-ku
    @time  decLU(A)
    
    spA = LuSpline(A, ku, kl; iscirc=iscirc, isLU=true)
    @time spB = LuSpline(n, t; iscirc=iscirc, isLU=true)

    dividende = 10000
    b = big.((rand(Int, n) .% dividende) .// dividende)
        x, y = sol(A, b)
        x2, y2 = sol(spB, b)
        @test y == y2
        @test x == x2
    
end


@testset "test decLU that is a tool for test" begin
    test_decLU(30)
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
test_splu(30, 9, true, true, type=BigFloat, tol=1e-70 )
test_splu(31, 10, true, true, type=BigFloat, tol=1e-70 )
test_splu(30, 11, false, true, type=BigFloat, tol=1e-70 )
test_splu(31, 10, false, true, type=BigFloat, tol=1e-70 )
end

@testset "test perf" begin
    test_splu(1000,9,true, true, type=Float64, tol=1e-12, perf=true)    
end
