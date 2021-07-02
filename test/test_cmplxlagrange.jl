
using Polynomials
using SemiLagrangian: getpoly
function test_cmplx(n, decref::CartesianIndex{2}=CartesianIndex(0,0))
#    fctref = Polynomial(big.(rationalize.(rand(n^2-1), tol=0.001)+ im*rationalize.(rand(n^2-1), tol=0.001)))
    tabref = [big(i+decref.I[1]+im*(j+decref.I[2])) + 0.1*(rand()+im*rand()-0.5-0.5*im) for i=1:n, j=1:n]
    moy = sum(tabref)/length(tabref)
    tabref = tabref .- moy

    p = getpoly(tabref)
    @show p
    for i=1:n, j=1:n
        res = p(tabref[i,j]) - i -im*j
#        @show n,decref,i,j,res
        @test isapprox(p( tabref[i,j] ), i +im*j, atol=1e-27)
    end
end

@testset "Lagrange complexe" begin
    for n=3:2:11
        println("")
        test_cmplx(n)
        test_cmplx(n, CartesianIndex(35,47))
        test_cmplx(n, CartesianIndex(-3557,4702))
    end
end

