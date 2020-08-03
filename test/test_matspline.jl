include("../src/bspline.jl")
include("../src/matspline.jl")
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

function test_splu(n, order, iscirc, isLU)
   t = getbspline(big(order),0).(1:order)
    A = topl(n, t, iscirc)
    ku = div(order,2)
    kl = order-1-ku
    if isLU
        decLU(A)
    end
    spA = LuSpline(A, ku, kl; iscirc=iscirc, isLU=isLU)
    spB = LuSpline(n, t; iscirc=iscirc, isLU=isLU)

    @test spA == spB

    dividende = 10000
    b = big.((rand(Int, n) .% dividende) .// dividende)
    if isLU
        x, y = sol(A, b)
        x2, y2 = sol(spB, b)
        @test y == y2
    end
end

function test_perf(n ,order, iscirc)
    t = getbspline(big(order),0).(1:order)
    A = topl(n, t, iscirc)
    ku = div(order,2)
    kl = order-1-ku
    @time  decLU(A)
    
    spA = LuSpline(A, ku, kl; iscirc=iscirc, isLU=isLU)
    @time spB = LuSpline(n, t; iscirc=iscirc, isLU=isLU)

    dividende = 10000
    b = big.((rand(Int, n) .% dividende) .// dividende)
    if isLU
        x, y = sol(A, b)
        x2, y2 = sol(spB, b)
        @test y == y2
        @test x == x2
    end
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
end
