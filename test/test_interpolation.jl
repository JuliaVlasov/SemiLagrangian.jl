
using LinearAlgebra
using DoubleFloats
using SemiLagrangian:
    AbstractInterpolation,
    Lagrange,
    B_SplineLU,
    B_SplineFFT,
    interpolate!,
    isbspline,
    getbspline,
    getprecal,
    get_allprecal,
    get_order,
    EdgeType,
    InsideEdge,
    CircEdge,
    interpolatemod!,
    getinverse,
    OpTuple

function test_interp(
    interp::AbstractInterpolation{Rational{BigInt},edge},
    dec,
    sz,
) where {edge}
    @time @testset "test interpolation  $interp dec=$dec" begin
        fct = if edge == CircEdge
            order = 3
            getbspline(big(order), 0)
        else
            order = 0
            x -> x^3 - x^2 - x // 6 + 1 // 4
        end
        mesh = collect((order + 1) * big.(0:(sz-1))) // sz
        deb = fct.(mesh)
        #       println("deb=$deb")
        ref = if edge == CircEdge
            decf = mod(dec, 1 // 1)
            decint = div(dec, 1 // 1)
            circshift(fct.(mesh .+ (order + 1) * decf // sz), (-decint,))
        else
            fct.(mesh .+ (order + 1) * dec // sz)
        end
        res = zeros(Rational{BigInt}, sz)
        interpolate!(res, deb, dec, interp)
        if edge == CircEdge
            res2 = zeros(Rational{BigInt}, sz)
            decint = convert(Int, floor(dec))
            decfloat = dec - decint
            precal = getprecal(interp, decfloat)
            # println("outside : precal=$precal decint=$decint decfloat=$decfloat")
            # println("order=$(get_order(interp)) order from precal=$(size(precal,1)-1)")   
            interpolate!(res2, deb, decint, precal, interp)
            println("norm=$(float(norm(res2-ref,Inf)))")
            if isbspline(interp)
                @test res2 == ref
            else
                @test size(filter(x -> x == 0, res2 - ref), 1) == 120
            end
        end
        println("norm=$(float(norm(res-ref,Inf)))")
        diff = res - ref
        # for i=1:sz
        #     println("i=$i diff=$(convert(Float64,diff[i])), $(diff[i])")
        # end
        fl = edge != CircEdge || isbspline(interp)
        println("fl=$fl isbspline=$(isbspline(interp))")
        if fl
            @test res == ref
        else
            @test size(filter(x -> x == 0, res - ref), 1) == 120
        end
    end
end
function test_interp(interp, sz)
    for i = 0:5
        println("test_interp i=$i")
        test_interp(interp, big"1235" // 10240 + i, sz)
    end
end

function test_inv0(
    t_interp::Vector{I},
    dec::NTuple{N,T},
    sz::NTuple{N,Int},
) where {T,N,I<:AbstractInterpolation{T,CircEdge}}
    res = ntuple(x -> zeros(T, sz), N)
    res2 = ntuple(x -> zeros(T, sz), N)
    for ind in CartesianIndices(sz)
        for i = 1:N
            res[i][ind] = T(ind.I[i] - 1)
        end
    end
    decbegin = ntuple(x -> dec[x] < 0 ? -Int(floor(dec[x])) : 0, N)
    decend = ntuple(x -> dec[x] > 0 ? Int(ceil(dec[x])) : 0, N)
    for i = 1:N
        interpolatemod!(res2[i], res[i], ind -> dec, t_interp, T(sz[i]), decbegin, decend)
    end

    resinv = ntuple(x -> zeros(T, sz), N)
    resinv2 = ntuple(x -> zeros(T, sz), N)

    calinverse!(resinv, res2)
    calinverse!(resinv2, resinv)

    note = resinv .- res

    note[1] .= mod.(note[1] .+ (sz[1] / 2 + dec[1]), sz[1]) .- sz[1] / 2
    note[2] .= mod.(note[2] .+ (sz[2] / 2 + dec[2]), sz[2]) .- sz[2] / 2

    @show norm.(note)

    for i = 1:N
        @test norm(note[i]) < 1e-12
        @test norm(res2[i] - resinv2[i]) < 1e-12
    end
end

function test_interp2d(
    t_interp::Vector{I},
    coeff::T,
    sz::Tuple{Int,Int},
) where {T,N,I<:AbstractInterpolation{T,CircEdge}}
    prec = 1000 * eps(T)
    @show prec
    ref = zeros(T, sz)
    ref2 = zeros(T, sz)
    cmplxref = zeros(Complex{T}, sz)
    opref = zeros(OpTuple{2,T}, sz)
    bufdec = zeros(Complex{T}, sz)
    opbufdec = zeros(OpTuple{2,T}, sz)
    res1 = zeros(T, sz)
    res2 = zeros(T, sz)
    res3 = zeros(T, sz)
    res4 = zeros(T, sz)
    opres1 = zeros(OpTuple{2,T}, sz)
    cmplxres2 = zeros(Complex{T}, sz)
    cmplxres3 = zeros(Complex{T}, sz)
    cmplxres4 = zeros(Complex{T}, sz)

    for ind in CartesianIndices(sz)
        x = 2T(pi) * sum(ind.I ./ sz)
        bufdec[ind] = coeff * (cos(x + 1) + im * cos(x + 2))
        opbufdec[ind] = coeff * OpTuple((cos(x + 1), cos(x + 2)))
        ref[ind] = sin(x + 3) + cos(x - 1)
        ref2[ind] = sin(x + 1) + 3cos(-x + 2) / 5
        cmplxref[ind] = coeff * (sin(x + 4) + im * cos(x + 5))
        opref[ind] = coeff * OpTuple((sin(x + 4), cos(x + 5)))
    end

    @time interpolate!(res1, ref, opbufdec, t_interp)
    @time interpolate!(res2, ref, bufdec, t_interp)
    @time interpolate!(res3, ref, ind -> reim(bufdec[ind]), t_interp)
    @time interpolate!(res4, ref2, bufdec, t_interp)

    norm2 = norm(res2 - res3)
    norm1 = norm(res1 - res3)

    @test norm1 < prec
    @test norm2 < prec
    @show norm1, norm2

    @time interpolate!(opres1, opref, opbufdec, t_interp)
    @time interpolate!(cmplxres2, cmplxref, bufdec, t_interp)
    @time interpolate!(cmplxres3, cmplxref, ind -> reim(bufdec[ind]), t_interp)
    @time interpolate!(cmplxres4, ref + im * ref2, bufdec, t_interp)

    optocmplx(v::OpTuple{2,T}) where {T} = v.v[1] + im * v.v[2]

    norm1 = norm(optocmplx.(opres1) - cmplxres3)
    norm2 = norm(cmplxres2 - cmplxres3)

    norm4 = norm(res2 + im * res4 - cmplxres4)

    @test norm1 < prec
    @test norm2 < prec
    @test norm4 < prec
    @show norm1, norm2, norm4
end

function test_inv(
    t_interp::Vector{I},
    coeff::T,
    sz::NTuple{N,Int},
) where {T,N,I<:AbstractInterpolation{T,CircEdge}}
    res = ntuple(x -> zeros(T, sz), N)
    res2 = ntuple(x -> zeros(T, sz), N)
    res3 = ntuple(x -> zeros(T, sz), N)
    dec = ntuple(x -> zeros(T, sz), N)
    resnew = ntuple(x -> zeros(T, sz), N)
    for ind in CartesianIndices(sz)
        for i = 1:N
            dec[i][ind] = coeff * cos(i + 2T(pi) * sum(ind.I ./ sz))
            res[i][ind] = T(ind.I[i] - 1)
        end
    end
    extr = extrema.(dec)
    decbegin = ntuple(x -> -Int(floor(extr[x][1])), N)
    decend = ntuple(x -> Int(ceil(extr[x][2])), N)
    for i = 1:N
        interpolatemod!(
            res2[i],
            res[i],
            ind -> ntuple(x -> dec[x][ind], N),
            t_interp,
            T(sz[i]),
            decbegin,
            decend,
        )
    end
    for i = 1:N
        interpolatemod!(
            res3[i],
            res2[i],
            ind -> ntuple(x -> dec[x][ind], N),
            t_interp,
            T(sz[i]),
            decbegin,
            decend,
        )
    end

    resinv = ntuple(x -> zeros(T, sz), N)
    resinv2 = ntuple(x -> zeros(T, sz), N)

    calinverse!(resinv, res2)

    calinverse!(resinv2, resinv)

    @show Float32.(res[1]) + im * Float32.(res[2])
    @show Float32.(res2[1]) + im * Float32.(res2[2])
    @show Float32.(resinv2[1]) + im * Float32.(resinv2[2])
    @show Float32.(dec[1]) + im * Float32.(dec[2])
    @show Float32.(resinv[1]) + im * Float32.(resinv[2])
    @show Float32.(resinv[1][CartesianIndex(10, 10):CartesianIndex(12, 12)]) +
          im * Float32.(resinv[2][CartesianIndex(10, 10):CartesianIndex(12, 12)])

    resfmr =
        ntuple(x -> mod.(res2[x] - res[x] - dec[x] .+ sz[x] / 2, sz[x]) .- sz[x] / 2, 2)
    resfmr2 =
        ntuple(x -> mod.(res3[x] - res2[x] - dec[x] .+ sz[x] / 2, sz[x]) .- sz[x] / 2, 2)
    @show norm(resfmr)
    @show norm(resfmr2)

    decnew = ntuple(x -> resinv[x] - res[x], 2)
    extr = extrema.(decnew)
    decbegin = ntuple(x -> -Int(floor(extr[x][1])), N)
    decend = ntuple(x -> Int(ceil(extr[x][2])), N)
    for i = 1:N
        interpolatemod!(
            resnew[i],
            res2[i],
            ind -> ntuple(x -> decnew[x][ind], N),
            t_interp,
            T(sz[i]),
            decbegin,
            decend,
        )
    end
    @show resnew

    for i = 1:N
        @test norm(res2[i] - resinv2[i]) < 1e-12
    end
end
function test_invfmr(
    t_interp::Vector{I},
    coeff::T,
    sz::NTuple{N,Int},
) where {T,N,I<:AbstractInterpolation{T,CircEdge}}
    res = ntuple(x -> zeros(T, sz), N)
    res2 = ntuple(x -> zeros(T, sz), N)
    res3 = ntuple(x -> zeros(T, sz), N)
    dec = ntuple(x -> zeros(T, sz), N)
    dec2 = ntuple(x -> zeros(T, sz), N)
    resnew = ntuple(x -> zeros(T, sz), N)
    for ind in CartesianIndices(sz)
        for i = 1:N
            dec[i][ind] = coeff * cos(i + 2T(pi) * sum(ind.I ./ sz))
            res[i][ind] = T(ind.I[i] - 1)
        end
    end
    extr = extrema.(dec)
    decbegin = ntuple(x -> -Int(floor(extr[x][1])), N)
    decend = ntuple(x -> Int(ceil(extr[x][2])), N)
    for i = 1:N
        interpolatemod!(
            res2[i],
            res[i],
            ind -> ntuple(x -> dec[x][ind], N),
            t_interp,
            T(sz[i]),
            decbegin,
            decend,
        )
    end
    @show res
    @show res2
    extr = extrema.(ntuple(x -> -dec[x], N))
    decbegin = ntuple(x -> -Int(floor(extr[x][1])), N)
    decend = ntuple(x -> Int(ceil(extr[x][2])), N)
    for i = 1:N
        interpolatemod!(
            dec2[i],
            -dec[i],
            ind -> ntuple(x -> dec[x][ind], N),
            t_interp,
            T(sz[i]),
            decbegin,
            decend,
        )
    end
    for i = 1:N
        dec2[i] .= mod.(dec2[i] .+ sz[i] / 2, sz[i]) .- sz[i] / 2
    end
    @show dec
    @show dec2
    extr = extrema.(dec2)
    decbegin = ntuple(x -> -Int(floor(extr[x][1])), N)
    decend = ntuple(x -> Int(ceil(extr[x][2])), N)

    for i = 1:N
        interpolatemod!(
            res3[i],
            res2[i],
            ind -> ntuple(x -> dec2[x][ind], N),
            t_interp,
            T(sz[i]),
            decbegin,
            decend,
        )
    end

    @show res2
    @show res3

    res4 = ntuple(x -> res[x] - res3[x], N)
    for i = 1:N
        res4[i] .= mod.(res4[i] .+ sz[i] / 2, sz[i]) .- sz[i] / 2
    end

    @show res4

    @show norm(res4)
    res5 = ntuple(x -> zeros(T, sz), N)
    ind = 1
    while norm(res4) > eps(Double64) * length(res4[1])
        if norm(res4) > 1
            res4 = ntuple(x -> res4[x] / 2, N)
        end
        extr = extrema.(res4)
        decbegin = ntuple(x -> -Int(floor(extr[x][1])), N)
        decend = ntuple(x -> Int(ceil(extr[x][2])), N)
        for i = 1:N
            interpolatemod!(
                res5[i],
                res3[i],
                ind -> ntuple(x -> res4[x][ind], N),
                t_interp,
                T(sz[i]),
                decbegin,
                decend,
            )
        end
        res4 = ntuple(x -> res[x] - res5[x], N)
        for i = 1:N
            res4[i] .= mod.(res4[i] .+ sz[i] / 2, sz[i]) .- sz[i] / 2
            res3[i] .= res5[i]
        end

        @show res4

        @show ind, norm(res4)

        ind += 1
    end
end
function test_getinv(
    t_interp::Vector{I},
    coeff::T,
    sz::NTuple{N,Int},
) where {T,N,I<:AbstractInterpolation{T,CircEdge}}
    res = ntuple(x -> zeros(T, sz), N)
    res2 = ntuple(x -> zeros(T, sz), N)
    res3 = ntuple(x -> zeros(T, sz), N)
    dec = ntuple(x -> zeros(T, sz), N)
    dec2 = ntuple(x -> zeros(T, sz), N)
    resnew = ntuple(x -> zeros(T, sz), N)
    for ind in CartesianIndices(sz)
        for i = 1:N
            dec[i][ind] = coeff * cos(i + 2T(pi) * sum(ind.I ./ sz))
            res[i][ind] = T(ind.I[i] - 1)
        end
    end

    sz_2 = div.(sz, 2)

    for i = 1:N
        interpolatemod!(res2[i], res[i], dec, t_interp, T(sz[i]))
    end

    dec2 = getinverse(dec, t_interp)

    for i = 1:N
        interpolatemod!(res3[i], res2[i], dec2, t_interp, T(sz[i]))

        res2[i] .= mod.(res[i] .- res3[i] .+ sz_2[i], sz[i]) .- sz_2[i]
    end

    note = norm(res2)
    @show note, sqrt(eps(T))
    @test note < sqrt(eps(T))
end

@testset "test inverse" begin
    T = Double64
    test_getinv([Lagrange(11, T), Lagrange(11, T)], T(0.00911), (128, 100))
    # test_inv0([Lagrange(11,T), Lagrange(11,T)], (zero(T),zero(T)), (20,30))
    # test_inv0([Lagrange(11,T), Lagrange(11,T)], (T(pi)/10,T(pi)/9), (20,30))
    # test_inv([Lagrange(11, T), Lagrange(11, T)], T(0.011), (20, 30))
end

function test_interpfloat(
    interp::AbstractInterpolation{T,edge},
    sz,
    tol,
    nb = 100,
) where {T,edge}
    tabdec =
        T.([
            big"0.345141526199181716726626262655544",
            -big"0.3859416191876155241320011187619",
            -big"1.28561390114441619187615524132001118762519",
            -big"0.885901390114441619187615524132001118762519",
            -big"5.678513256790098898776656565545454544544545", # only for edge==CircEdge
            big"4.9876651456677809099887665655556565565656565", # only for edge==CircEdge
            big"0.186666659416191876155241320011187619",
            big"0.590999232323232323232365566787878898898",
            big"1.231098015934444444444444788888888878878",
        ])
    mesh = T.(collect(big.(0:(sz-1))) / sz)

    fct1(x) = T(cos(2big(pi) * x + big"0.25"))
    fct2(x) = T(exp(-(cos(2big(pi) * x + big"0.25") - 1)^2))
    tabfct = [fct1, fct2]
    @time @testset "test interpolation  $interp" begin
        nmax = 0

        for (i_fct, fct) in enumerate(tabfct), dec in tabdec
            deb = fct.(mesh)
            #       println("deb=$deb")
            fp = deb
            fi = zeros(T, sz)

            decint = convert(Int, floor(dec))
            value = dec - decint
            if get_order(interp) % 2 == 0 && value > 0.5
                value -= 1
                decint += 1
            end
            if edge != CircEdge && abs(decint) > 2
                continue
            end
            precal = if edge == CircEdge
                getprecal(interp, value)
            else
                get_allprecal(interp, decint, value)
            end
            for i = 1:nb
                fi .= fp
                ref = fct.(mesh .+ i * dec / sz)
                interpolate!(fp, fi, decint, precal, interp)
                nmax = max(nmax, float(norm(fp - ref, Inf)))
                if isbspline(interp) && edge != CircEdge
                    @show typeof(interp), i_fct, dec, sz, nb, nmax
                else
                    @test isapprox(fp, ref, atol = tol)
                end
            end
        end
        @show typeof(interp), sz, nb, nmax
    end
end
T = Double64
test_interp2d([Lagrange(11, T), Lagrange(11, T)], T(1.25), (128, 100))

test_interp(Lagrange(3, Rational{BigInt}; edge = InsideEdge), big"3" // 1024, 128)

test_interp(Lagrange(3, Rational{BigInt}; edge = CircEdge), 128)

test_interp(B_SplineLU(3, 128, Rational{BigInt}), 128)

test_interpfloat(Lagrange(3, BigFloat; edge = CircEdge), 128, 1e-3, 100)
test_interpfloat(Lagrange(3, Float64; edge = CircEdge), 128, 1e-3, 100)

test_interpfloat(Lagrange(7, BigFloat; edge = InsideEdge), 128, 1e-3, 3)
test_interpfloat(Lagrange(7, Float64; edge = InsideEdge), 128, 1e-3, 3)

test_interpfloat(Lagrange(21, BigFloat; edge = CircEdge), 256, 1e-20)
test_interpfloat(Lagrange(9, Float64; edge = CircEdge), 256, 1e-10)

test_interpfloat(Lagrange(4, BigFloat; edge = CircEdge), 256, 1e-5)
test_interpfloat(Lagrange(4, Float64; edge = CircEdge), 256, 1e-5)

test_interpfloat(Lagrange(22, BigFloat; edge = CircEdge), 256, 1e-20)
test_interpfloat(Lagrange(12, Float64; edge = CircEdge), 256, 1e-10)

test_interpfloat(B_SplineLU(3, 256, BigFloat), 256, 1e-5)
test_interpfloat(B_SplineLU(3, 256, Float64), 256, 1e-5)

test_interpfloat(B_SplineLU(21, 256, BigFloat), 256, 1e-30)
test_interpfloat(B_SplineLU(11, 256, Float64), 256, 1e-12)

test_interpfloat(B_SplineFFT(3, 256, BigFloat), 256, 1e-5)
test_interpfloat(B_SplineFFT(3, 256, Float64), 256, 1e-5)

test_interpfloat(B_SplineFFT(21, 256, BigFloat), 256, 1e-30)
test_interpfloat(B_SplineFFT(11, 256, Float64), 256, 1e-12)
