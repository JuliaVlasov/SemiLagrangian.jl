include("../src/lapack.jl")
using Test
using LinearAlgebra
using LinearAlgebra.LAPACK: pttrf!, pttrs!, gbtrf!, gbtrs!
using Random

function test_ptt()
    @testset "functions pttrf and pttrs" begin
        D=[1.,2,3,4,5,6,7]
        E=[0.1,0.2,0.3,0.4,0.5,0.6]
        D1 = copy(D)
        D2 = copy(D)
        E1 = copy(E)
        E2 = copy(E)
        res = pttrf!(D1,E1)
        resgen = pttrfgen!(D2,E2)
        @test res == resgen
        @test D1 == D2
        @test E1 == E2

        B = [10.,11.,12,13,14,15,16]
        B1 = copy(B)
        B2 = copy(B)
        res = pttrs!(D1, E1, B1)
        resgen = pttrsgen!(D2, E2, B2)
        @test res == resgen
        @test B1 == B2
    end
end
function test_pttbig()
    @testset "functions pttrf and pttrs with BigFloat" begin
        D=[1.,2,3,4,5,6,7]
        E=[0.1,0.2,0.3,0.4,0.5,0.6]
        D1 = copy(D)
        D2 = big.(D)
        E1 = copy(E)
        E2 = big.(E)
        res = pttrf!(D1,E1)
        resgen = pttrfgen!(D2,E2)
        @test isapprox(res[1], resgen[1], atol=1e-15, rtol=1e-14)
        @test isapprox(res[2], resgen[2], atol=1e-15, rtol=1e-14)
        @test isapprox(D1, D2, atol=1e-15, rtol=1e-14)
        @test isapprox(E1, E2, atol=1e-15, rtol=1e-14)

        B = [10.,11.,12,13,14,15,16]
        B1 = copy(B)
        B2 = big.(B)
        res = pttrs!(D1, E1, B1)
        resgen = pttrsgen!(D2, E2, B2)
        @test isapprox(res, resgen, atol=1e-15, rtol=1e-14)
        @test isapprox(B1, B2, atol=1e-15, rtol=1e-14)
    end
end

function luPDecompose( A )
    n = size(A,1)
    p = collect(1:n)
    kl=3
    for i=1:n
        maxA=abs(A[i,i])
        imax=i
        borne = min(i+kl,n)
        for k=i+1:borne
            if ( (absA = abs(A[i,k])) > maxA )
                maxA = absA
                imax = k
            end
        end
        if imax != i
            pnew = p[imax]
            Anew = A[imax,:]
            for j=imax:-1:i+1
                p[j] = p[j-1]
                A[j,:] = A[j-1,:]
            end
            p[i] = pnew
            A[i,:] = Anew
        end
        for j=i+1:n
            A[j,i] /= A[i,i]
            for k=i+1:n
                A[j,k] -= A[j,i] * A[i,k]
            end
        end
    end
    return p
end

# INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
# OUTPUT: x - solution vector of A*x=b
#
function LUPSolve( A, p, b)
    n=size(A,1)
    bcp =copy(b)
    for i=1:n
        b[i] = bcp[p[i]]
        for k=1:i-1
            b[i] -= A[i,k] * b[k]
        end
    end

    for i = n:-1:1
        for k=i+1:n
            b[i] -= A[i,k] * b[k]
        end
        b[i] = b[i] / A[i,i]
    end
end

function testLUP()
    n=1000
   A = rand(n,n)
     ku=2
    kl=3
    for i=1:n, j=1:n
        if j > i+ku || j < i-kl
            A[i,j] = 0
        end
    end
 
    Acp = copy(A)
    P=luPDecompose(A)

    b = rand(n);
    bcp = copy(b);

    LUPSolve(A,P,b)

    res = norm(Acp*b-bcp)
    println("res=$res")

    kuMax = klMax = 0
    for i=1:n, j=1:n
        if A[i,j] != 0
            if j>i && j-i> klMax
                klMax = j-i
            end
            if  j<i && i-j> kuMax
                kuMax = i-j
            end
        end
    end
    println("kuMax=$kuMax klMax=$klMax")
end


function test_gbt(isbig, seedval)
    @testset "functions gbtrf isbig=$isbig" begin
        n = 10
        Random.seed!(seedval)
        A= rand(n,n)
        ku=3
        kl=4
        for i=1:n, j=1:n
            if j > i+ku || j < i-kl
                A[i,j] = 0
            end
        end
        AB = zeros(2kl+ku+1,n)
        for j=1:n
            for i = max(1,j-ku):min(n,j+kl)
                AB[kl+ku+1+i-j,j] = A[i,j]
            end
        end

        AB1 = copy(AB)
        AB2 = isbig ? big.(AB) : copy(AB)

        res1 = gbtrf!(kl,ku,n,AB1)
        res2 = gbtrfgen!(kl,ku,n,AB2)

        p1 = res1[2]
        p2 = res2[2]

        @test p1 == p2

        println("norm(AB1-AB2)=$(norm(AB1-AB2))")

        @test isapprox(AB1,AB2,atol=1e-15,rtol=1e-15) 

        B = rand(n,3)
        B1 = copy(B)
        B2 = isbig ? big.(B) : copy(B)

        res_1 = gbtrs!('N', kl, ku, n, AB1, p1, B1)
        res_2 = gbtrsgen!('N', kl, ku, n, AB2, p1, B2)

        println("norm(B1-B2)=$(norm(B1-B2,Inf))")
        @test isapprox(B1, B2, atol=1e-14, rtol=1e-14)
    end
end

test_ptt()
test_pttbig()
test_gbt(false,6789)
test_gbt(true,1234)
# testLUP()
