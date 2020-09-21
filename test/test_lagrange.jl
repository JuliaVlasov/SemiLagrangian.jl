using Test
using LinearAlgebra
include("../src/advection.jl")
include("../src/lagrange.jl")
function test_lagrange(type::DataType, order, iscirc::Bool, number,  tol)
    @testset "Lagrange interpolation function type=$type order=$order iscirc=$iscirc" begin
        lag = Lagrange(type, order; iscirc=iscirc)
#        println("lag=$lag")
        n = number
        # whennot circuar coef != 1 set the function unperiodic
        coef = iscirc ? 1. : 1.111
        fct(v,n) = cos(2big(pi)*coef*v/n)
        tf = fct.(BigFloat.(collect(1:n)),n)
        value = big"0.456231986"
        normax=0.0
        for i=1:number
            pol = polinterpol(lag, tf, i)
            ref=fct(i+value,n)
            res=pol(value)
            normax = max(normax, norm(ref-res))
#            println("i=$i norm=$(norm(res-ref,Inf))")
#            println("i=$i norm=$(norm(res-ref,Inf)) type(pol)=$(typeof(pol)) pol=$pol")
            @test isapprox(ref, res, atol=tol)
        end
        println("normax=$normax")
    end
end

function test_interpolation(type::DataType, order, iscirc::Bool, number,  tol, nb=1)
    
    lag = Lagrange(type, order; iscirc=iscirc)
    coef = iscirc ? 1. : 1.111
    n = number
 #   fct(v,n) = exp( -cos(2big(pi)*coef*v/n)^2)
 #    fct(v,n) = exp( -(50*(v-n/2)/n)^2)
    fct(v,n) = exp( -(0.02*cos(2big(pi)*coef*v/n))^2)
 fp = fct.(BigFloat.(collect(1:n)),n)
    fi = zeros(type, number)
    value = big"0.38571390114441619187615524132001118762519"
    for i=1:nb
        fi .= fp
        fpref = fct.(BigFloat.(collect(1:n)).+i*value,n)
        interpolate!(fp, fi, value, lag)
        println("i=$i norm = $(norm(fpref-fp,Inf))")
        # for i=1:number
        #     println("i=$i norm=$(norm(fpref[i]-fp[i]))")
        # end
        println("i=$i norm=$(norm(fpref-fp))")
        @test isapprox(fpref, fp, atol=tol)
    end
end

# test_lagrange(BigFloat, 17,true, 1000, 1e-45)
# test_lagrange(BigFloat, 23,true, 1000, 1e-50)
# test_lagrange(BigFloat, 17,false, 1000, 1e-38)
# test_lagrange(BigFloat, 23,false, 1000, 1e-50)
# @testset "Lagrange interpolation order=17 iscirc=true" begin
#     @time test_interpolation(BigFloat, 17,true, 1000, 1, 1e-20, 10000)
#     for i= 2:13
#         @time test_interpolation(BigFloat, 17,true, 1000, i, 1e-20, 100)
#     end
# end

@time test_interpolation(BigFloat, 43,true, 200, 1e-1, 200 )
# test_interpolation(BigFloat, 23,true, 1000, 12, 1e-40)
# test_interpolation(BigFloat, 17,false, 1000, 10, 1e-38)
# test_interpolation(BigFloat, 23,false, 1000, 12, 1e-40)






