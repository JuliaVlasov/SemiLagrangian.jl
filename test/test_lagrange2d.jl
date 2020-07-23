using Test
using LinearAlgebra
include("../src/lagrange2d.jl")
function test_lagrange(type::DataType, order, iscirc::Bool, number,  tol)
    @testset "Lagrange interpolation 2d function type=$type order=$order iscirc=$iscirc" begin
        lag = Lagrange2d(type, order; iscirc=iscirc)
#        println("lag=$lag")
        n = number
        # whennot circuar coef != 1 set the function unperiodic
        coef = iscirc ? 1. : 1.111
        fct(x,y,n) = cos(2big(pi)*coef*x/n+0.43)*sin(2big(pi)*coef*y/n+0.2)
        tf = [fct(i,j,n) for i=1:n, j=1:n]
        value_x = big"0.456231986"
        value_y = big"0.3676956102"
        normax=0.0
        for i=1:number, j=1:number
            pol = polinterpol(lag, tf, (i,j))
            ref=fct(i+value_x,j+value_y,n)
            res=pol(value_x, value_y)
            normax = max(normax, norm(ref-res))
#            println("i=$i norm=$(norm(res-ref,Inf))")
#            println("i=$i norm=$(norm(res-ref,Inf)) type(pol)=$(typeof(pol)) pol=$pol")
            @test isapprox(ref, res, atol=tol)
        end
        println("normax=$normax")
    end
end

function test_interpolation(type::DataType, order, iscirc::Bool, number, granularity,  tol, nb=1)
    
    lag = LagrangeNew(type, order; iscirc=iscirc, granularity=granularity)
    coef = iscirc ? 1. : 1.111
    n = number
    fct(v,n) = exp( -cos(2big(pi)*coef*v/n)^2)
    fp = fct.(BigFloat.(collect(1:n)),n)
    fi = zeros(type, number)
    value = big"0.485713901"
    for i=1:nb
        fi .= fp
        fpref = fct.(BigFloat.(collect(1:n)).+i*value,n)
        interpolate!(missing, fp, fi, value, lag)
        println("i=$i granularity=$granularity norm = $(norm(fpref-fp,Inf))")
        # for i=1:number
        #     println("i=$i norm=$(norm(fpref[i]-fp[i]))")
        # end
        @test isapprox(fpref, fp, atol=tol)
    end
end

test_lagrange(BigFloat, 7,true, 100, 1e-45)
test_lagrange(BigFloat, 5,true, 100, 1e-50)
# test_lagrange(BigFloat, 17,false, 1000, 1e-38)
# test_lagrange(BigFloat, 23,false, 1000, 1e-50)
# @testset "Lagrange interpolation order=17 iscirc=true" begin
#     for i= 1:25
#         @time test_interpolation(BigFloat, 17,true, 1000, i, 1e-20)
#     end
# end
# @time test_interpolation(BigFloat, 17,true, 200, 1, 1e-10, 200 )
# test_interpolation(BigFloat, 23,true, 1000, 12, 1e-40)
# test_interpolation(BigFloat, 17,false, 1000, 10, 1e-38)
# test_interpolation(BigFloat, 23,false, 1000, 12, 1e-40)






