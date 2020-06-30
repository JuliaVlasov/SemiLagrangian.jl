include("../src/bspline.jl")
using Test

function test_add()
    @time @testset "test spline object" begin
        p0 = getbspline(20,0)
        p3 = getbspline(20,3)
        @test isapprox(sum(p0.(range(big"1.0",length=21))),1, atol=1e-60)
        @test sum(p0.(range(one(Rational{BigInt}),length=21))) == 1//1
        for i = 0:20
            x = rand(BigFloat)
            res0=p0(x+i)
            res3=p3(x+i+3)
            @test isapprox( res0, res3,atol=1e-60)
        end
    end
end

test_add()

