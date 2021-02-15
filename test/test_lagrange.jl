using LinearAlgebra
using Polynomials
using SemiLagrangian: Lagrange, get_precal, get_allprecal,  interpolate!, get_fact_order


function test_base_lagrange(order)
    lag = Lagrange(order)
    tab = rationalize.(BigInt, rand(order), tol=1/1000000)
    fct = Polynomial(tab)
    dec = div(order, 2)
    for i = 1:10
        value = rationalize.(BigInt, rand(), tol=1/1000000)
        res = 0//1
        for j=0:order
            res += lag.lagpol[j+1](value)*fct(j-dec)
        end
        res /= get_fact_order(lag)
        @test res == fct(value)
    end
end

@time @testset "test base interpolation" begin
    for ord=3:27
        test_base_lagrange(ord)
    end
end




