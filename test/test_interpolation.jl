
using LinearAlgebra
using SemiLagrangian: InterpolationType, Lagrange, B_SplineLU, B_SplineFFT, get_type, 
    interpolate!, isbspline, getbspline, get_precal, get_allprecal, get_order, EdgeType, InsideEdge, CircEdge

function test_interp(interp::InterpolationType{Rational{BigInt}, edge}, dec,  sz) where {edge}

    @time @testset "test interpolation  $interp dec=$dec" begin    
        fct = if edge == CircEdge
            order=3
            getbspline(big(order),0)
        else
            order=0
            x -> x^3-x^2-x//6+1//4
        end
        mesh=collect((order+1)*big.(0:(sz-1)))//sz
        deb = fct.(mesh)
 #       println("deb=$deb")
        ref=if edge == CircEdge
            decf = mod(dec,1//1)
            decint = div(dec,1//1)
            circshift(fct.(mesh .+ (order+1)*decf//sz), (-decint,))
        else
            fct.(mesh .+ (order+1)*dec//sz)
        end
        res = zeros(Rational{BigInt}, sz)
        interpolate!(res, deb, dec, interp)
        if edge == CircEdge
            res2 = zeros(Rational{BigInt}, sz)
            decint = convert(Int, floor(dec))
            decfloat = dec - decint
            precal = get_precal(interp,decfloat)
            # println("outside : precal=$precal decint=$decint decfloat=$decfloat")
            # println("order=$(get_order(interp)) order from precal=$(size(precal,1)-1)")   
            interpolate!(res2, deb, decint, precal, interp)
            println("norm=$(float(norm(res2-ref,Inf)))")
            if isbspline(interp)
                @test res2 == ref
            else
                @test size(filter( x->x==0,res2-ref),1)==120
            end
        end
        println("norm=$(float(norm(res-ref,Inf)))")
        diff = res-ref
        # for i=1:sz
        #     println("i=$i diff=$(convert(Float64,diff[i])), $(diff[i])")
        # end
        fl = edge != CircEdge || isbspline(interp)
        println("fl=$fl isbspline=$(isbspline(interp))")
        if fl
            @test res == ref
        else
            @test size(filter( x->x==0,res-ref),1)==120
        end

    end
end
function test_interp(interp, sz)
    for i=0:5
        println("test_interp i=$i")
        test_interp(interp, big"1235"//10240 +i, sz)
    end
end
# function test_interpolation2(T::DataType, order, edge::EdgeType, number,  nb, tol, islu=true)
    
#     n = number
#     sp = if (islu)
#         B_SplineLU(order, n, zero(T); edge=edge)
#     else
#         B_SplineFFT(order, n, zero(T))
#     end
#  #   fct(v,n) = exp( -cos(2big(pi)*coef*v/n)^2)
#  #    fct(v,n) = exp( -(50*(v-n/2)/n)^2)
#     fct1(v,n) = exp( -(2*cos(2T(big(pi)*v/n)))^2)
#     fct2(v,n)=cos(2T(big(pi)*v/n))
#     tabfct = [fct1, fct2]

#     tabv = T.([            big"0.186666659416191876155241320011187619",
#                     -big"1.58561390114441619187615524132001118762519",
#     -big"1.28561390114441619187615524132001118762519",
#     -big"0.885901390114441619187615524132001118762519",
#             -big"0.3859416191876155241320011187619",
#            big"0.590999232323232323232365566787878898898",
#             big"1.231098015934444444444444788888888878878"
#         ])
#     ifct=0
#     for fct in tabfct
#         ifct += 1
#         fi = zeros(T, number)
#         fp = zeros(T, number)
#         @show typeof(fp), ifct
#         ival=0
#         for valuebig in tabv
#             ival += 1
#             decint = convert(Int,floor(valuebig))
#             value = valuebig-decint
#             if order%2 == 0
#                 if value < 0.5
#                     value -= 1
#                     decint += 2
#                 else
#                     value -= 0
#                     decint += 1
#                 end
#             end
#             precal = get_precal(sp, value)
#             nmax=0
#             fp .= fct.(T.(collect(1:n)),n)
# #            @show typeof(fp), ifct, ival
#             for i=1:nb
#                 fi .= fp
#                 fpref = fct.(T.(collect(1:n)) .+ i*valuebig, n)
# #                @show typeof(fp), ifct, ival, i
#                 interpolate!(fp, fi, decint, precal, sp)

#                 nmax = max(nmax,norm(fpref-fp))
#                 if order%2 == 0
#                     @show i, nmax, ival
#                 end
# #                @test isapprox(fpref, fp, atol=tol)
#             end
#             println("order = $order value=$valuebig,nmax=$nmax ifct=$ifct")
#         end
#     end
# end

function test_interpfloat(interp::InterpolationType{T, edge}, sz, tol, nb=100) where {T,edge}

    tabdec = T.([big"0.345141526199181716726626262655544",
    -big"0.3859416191876155241320011187619",
    -big"1.28561390114441619187615524132001118762519",
    -big"0.885901390114441619187615524132001118762519",
    -big"5.678513256790098898776656565545454544544545", # only for edge==CircEdge
    big"4.9876651456677809099887665655556565565656565", # only for edge==CircEdge
    big"0.186666659416191876155241320011187619",
    big"0.590999232323232323232365566787878898898",
    big"1.231098015934444444444444788888888878878"
])
    mesh=T.(collect(big.(0:(sz-1)))/sz)

    fct1(x) = T(cos(2big(pi)*x+big"0.25"))
    fct2(x) = T(exp(-(cos(2big(pi)*x+big"0.25")-1)^2))
    tabfct=[fct1, fct2]
    @time @testset "test interpolation  $interp" begin    
        nmax=0

        for (i_fct, fct) in enumerate(tabfct), dec in tabdec
            deb = fct.(mesh)
    #       println("deb=$deb")
            fp = deb
            fi = zeros(T, sz)

            decint = convert(Int,floor(dec))
            value = dec-decint
            if get_order(interp)%2 == 0 && value > 0.5 
                value -= 1
                decint += 1
            end
            if edge != CircEdge && abs(decint) > 2
                continue
            end
            precal = edge == CircEdge ? get_precal(interp, value) : get_allprecal(interp, decint, value)
            for i=1:nb
                fi .= fp
                ref=fct.(mesh .+ i*dec/sz)   
                interpolate!(fp, fi, decint, precal, interp)
                nmax = max(nmax, float(norm(fp-ref,Inf)))
                if isbspline(interp) && edge != CircEdge
                    @show typeof(interp), i_fct, dec, sz, nb, nmax 
                else
                    @test isapprox(fp, ref, atol=tol)
                end
            end
        end
        @show typeof(interp), sz, nb, nmax 
    end
end


test_interp(Lagrange(3, Rational{BigInt}; edge=InsideEdge),big"3"//1024,128)

test_interp(Lagrange(3, Rational{BigInt}; edge=CircEdge), 128)

test_interp(B_SplineLU(3, 128, Rational{BigInt}), 128)

test_interpfloat(Lagrange(3, BigFloat,edge=CircEdge), 128, 1e-3, 100)
test_interpfloat(Lagrange(3, Float64,edge=CircEdge), 128, 1e-3, 100)

test_interpfloat(Lagrange(7, BigFloat,edge=InsideEdge), 128, 1e-3, 3)
test_interpfloat(Lagrange(7, Float64,edge=InsideEdge), 128, 1e-3, 3)

test_interpfloat(Lagrange(21, BigFloat, edge=CircEdge), 256, 1e-20)
test_interpfloat(Lagrange(9, Float64, edge=CircEdge), 256, 1e-10)

test_interpfloat(Lagrange(4, BigFloat,edge=CircEdge), 256, 1e-5)
test_interpfloat(Lagrange(4, Float64,edge=CircEdge), 256, 1e-5)

test_interpfloat(Lagrange(22, BigFloat,edge=CircEdge),256,1e-20)
test_interpfloat(Lagrange(12, Float64,edge=CircEdge),256,1e-10)

test_interpfloat(B_SplineLU(3,256,BigFloat), 256, 1e-5)
test_interpfloat(B_SplineLU(3,256,Float64), 256, 1e-5)

test_interpfloat(B_SplineLU(21,256,BigFloat), 256, 1e-30)
test_interpfloat(B_SplineLU(11,256,Float64), 256, 1e-12)

test_interpfloat(B_SplineFFT(3, 256, BigFloat), 256, 1e-5)
test_interpfloat(B_SplineFFT(3, 256, Float64), 256, 1e-5)

test_interpfloat(B_SplineFFT(21, 256, BigFloat), 256, 1e-30)
test_interpfloat(B_SplineFFT(11, 256, Float64), 256, 1e-12)


# test_interpfloat(B_SplineLU(7,1024,BigFloat; edge=InsideEdge),1024, 1e-4, 1)

# test_interpfloat(B_SplineLU(21,1024,BigFloat; edge=InsideEdge),1024, 1e-18, 5)
