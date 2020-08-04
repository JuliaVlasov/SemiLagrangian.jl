include("../src/bspline.jl")
include("../src/matspline.jl")

using LinearAlgebra
using Test

function test_spline(order, prec)
    setprecision(prec) do
        @time @testset "test spline object of order $order" begin
           p0 = getbspline(order,0)
            p3 = getbspline(order,3)
            @test p3 == decal(p0,3)
            @test isapprox(sum(p0.(range(big"1.0",length=order+1))),1, atol=1e-60)
            @test sum(p0.(range(one(Rational{BigInt}),length=order+1))) == 1             
            for i = 0:order
                x = rand(BigFloat)
                res0=p0(x+i)
                res3=p3(x+i+3)
                @test isapprox( res0, res3,atol=1e-60)
                @test order == Polynomials.degree(p0[i])
                pol1 = p0[i]
                pol2 = p0[i+1]
                for j=1:order
                    @test pol1(i+1) == pol2(i+1)
                    pol1 = derivative(pol1)
                    pol2 = derivative(pol2)
                end
            end
        end
    end
end
for order = 1:15
    test_spline( order, order < 20 ? 256 : 512 )
end

function test_interpolation(type::DataType, order, iscirc::Bool, n,  tol)
    
    sp = BSplineNew(order, n, zero(type); iscirc=iscirc)
    coef = convert(type, iscirc ? 1 : 1.111)
    # fct(v,n) = exp( -cos(2pi*coef*v/n)^2)
    fct(v,n) = cos(2pi*coef*v/n)
    fp = fct.(convert.(type,(collect(1:n))),n)
    fi = zeros(type, n)
    value = convert(type, big"0.38571390114441619187615524132001118762519")
    for i=1:n
        fi .= fp
        fpref = fct.(convert.(type,(collect(1:n))) .+ i*value,n)
        interpolate!(missing, fp, fi, value, sp)
        # for i=1:number
        #     println("i=$i norm=$(norm(fpref[i]-fp[i]))")
        # end
        println("i=$i norm=$(norm(fpref-fp))")
        @test isapprox(fpref, fp, atol=tol)
    end
end

@testset "test interpolation bspline" begin
    test_interpolation(BigFloat, 11, true, 100, 1e-8)
    test_interpolation(BigFloat, 21, true, 100, 1e-14)
end
