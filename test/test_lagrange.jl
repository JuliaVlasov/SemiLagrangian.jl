using LinearAlgebra
using Polynomials
using SemiLagrangian: Lagrange


function test_base_lagrange(order)
    lag = Lagrange(order, Rational{BigInt})
    tab = rationalize.(BigInt, rand(order), tol=1/1000000)
    fct = Polynomial(tab)
    dec = div(order, 2)
    for i = 1:10
        value = rationalize.(BigInt, rand(), tol=1/1000000)
        res = 0//1
        for j=0:order
            res += lag.tabfct[j+1](value)*fct(j-dec)
        end
        @test res == fct(value)
    end
end

@time @testset "test base interpolation" begin
    for ord=3:27
        test_base_lagrange(ord)
    end
end




