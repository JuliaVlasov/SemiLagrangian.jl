using LinearAlgebra
using Polynomials
using SemiLagrangian: Lagrange, get_precal, interpolate!


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
        res /= lag.fact_order
        @test res == fct(value)
    end
end

function test_interpolation(T::DataType, order, iscirc::Bool, number,  tol, nb=1)
    
    lag = Lagrange(order, T; iscirc=iscirc)
    n = number
 #   fct(v,n) = exp( -cos(2big(pi)*coef*v/n)^2)
 #    fct(v,n) = exp( -(50*(v-n/2)/n)^2)
    fct(v,n) = exp( -(2*cos(2T(big(pi)*v/n)))^2)

    tabv = T.([-big"1.28561390114441619187615524132001118762519",
            -big"0.885901390114441619187615524132001118762519",
            -big"0.3859416191876155241320011187619",
            big"0.186666659416191876155241320011187619",
            big"0.590999232323232323232365566787878898898",
            big"1.231098015934444444444444788888888878878"
        ])
    fi = zeros(T, number)
    fp = zeros(T, number)
    for valuebig in tabv
        decint = convert(Int,floor(valuebig))
        value = valuebig-decint
        if order%2 == 0 && value > 0.5
            value -= 1
            decint += 1
        end
        precal = get_precal(lag, value)
        nmax=0
        fp .= fct.(T.(collect(1:n)),n)
        for i=1:nb
            fi .= fp
            fpref = fct.(T.(collect(1:n)) .+ i*valuebig, n)
            interpolate!(fp, fi, decint, precal, lag)

            nmax = max(nmax,norm(fpref-fp))
           @test isapprox(fpref, fp, atol=tol)
        end
#        println("order = $order value=$valuebig,nmax=$nmax")
    end
end

@time @testset "test base interpolation" begin
    for ord=3:27
        test_base_lagrange(ord)
    end
end


@time @testset "Interpolation Lagrange" begin
    test_interpolation(Float64, 3,true, 200, 1e-3, 10 )
    test_interpolation(Float64, 4,true, 200, 1e-4, 10 )
    test_interpolation(Float64, 5,true, 200, 1e-5, 10 )
    test_interpolation(BigFloat, 23,true, 1000, 1e-34, 10)
    test_interpolation(BigFloat, 47,true, 1000, 1e-59, 10)
end

# TODO cas non circulire
# test_interpolation(BigFloat, 23,false, 1000, 12, 1e-40)

#




