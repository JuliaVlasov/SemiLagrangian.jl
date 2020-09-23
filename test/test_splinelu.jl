include("../src/advection.jl")
include("../src/bspline.jl")
include("../src/bsplinelu.jl")

using LinearAlgebra
using Random
using Test
Random.seed!(5431221)
# create a band or circular matrix from a vector of non-zero data
function topl(n, t, iscirc=true)
    res=zeros(Rational{BigInt},n,n)
    kl, ku = get_kl_ku(size(t,1))
    for i=1:n
        for (j,v) in enumerate(t)
            ind = i+j-kl-1
            if 1 <= ind <= n
                res[i, ind] = v
            elseif iscirc
                if ind < 1
                    ind += n
                else
                    ind -= n
                end
                res[i,ind] = v
            end
        end
    end
    return res
end
# inplace lu decomposition
function decLU( A::Matrix{T} ) where{T}
    n = size(A,1)
    for k=1:n
        pivot = A[k,k]
        for i=k+1:n
            A[i,k] /= pivot
        end
        for i=k+1:n, j=k+1:n
            # if A[i,k] != 0 && A[k,j] != 0
            #     println("trace avant k=$k i=$i j=$j A[i,j]=$(A[i,j]) A[i,k]=$(A[i,k]) A[k,j]=$(A[k,j])")
            # end
            A[i,j] -= A[i,k]*A[k,j]
            # if A[i,k] != 0 && A[k,j] != 0
            #     println("trace après k=$k i=$i j=$j A[i,j]=$(A[i,j])")
            # end
        end
    end
    return A
end
function getLU(A::Matrix{T}) where{T}
    n = size(A,1)
    L = zeros(T,n,n)
    U = zeros(T,n,n)
    for i=1:n
        L[i,i] = 1
        for j=1:i-1
            L[i,j] = A[i,j]
        end
        for j=i:n
            U[i,j] = A[i,j]
        end
    end
    return L, U
end

# function testTools()

function sol( A::Matrix{T}, Y::Vector{T}) where{T}
    L, U = getLU(A)
    n = size(A,1)
    Y1 = zeros(T,n)
    for i=1:n
        Y1[i] = Y[i] - sum(Y1[1:i-1] .* A[i, 1:i-1])
    end
 #   @assert Y1 == (L^(-1))*Y "Erreur 1" 
 #   @assert isapprox(Y1, (L^(-1))*Y, atol=1e-60) "Erreur 1" 
    X = zeros(T,n)
    for i=n:-1:1
        X[i] = (Y1[i] - sum(X[i+1:n] .* A[i, i+1:n]))/A[i,i]
    end
#    @assert X == (U^(-1))*Y1 "Erreur 2" 
#    @assert isapprox(X,(U^(-1))*Y1,atol=1e-60) "Erreur 2" 
#    @assert X == ((L*U)^(-1))*Y "Erreur3"
#    @assert isapprox(X, ((L*U)^(-1))*Y, atol=1e-60) "Erreur3"
    return X, Y1
end



function test_decLU(n)
    dividende=10000
    A = big.((rand(Int, n, n) .% dividende) .// dividende)
    B = copy(A)
    decLU(A)
    L, U = getLU(A)
    @test L*U == B
    b = big.((rand(Int, n) .% dividende) .// dividende)
    x, _ = sol(A, b)
    @test B*x == b
end

function test_splu(n, order, iscirc, isLU; type=Rational{BigInt}, tol=NaN, perf=false)
   t = convert.(type, getbspline(big(order),0).(1:order))
    A = convert.(type, topl(n, t, iscirc))
    Acp = copy(A)
    ku = div(order,2)
    kl = order-1-ku
    if isLU
        perf && @time decLU(A)
        !perf && decLU(A)
    end
    spA = LuSpline(A, ku, kl; iscirc=iscirc, isLU=isLU)
    perf && @time spB = LuSpline(n, t; iscirc=iscirc, isLU=isLU)
    !perf && (spB = LuSpline(n, t; iscirc=iscirc, isLU=isLU))
    if isnan(tol)
        @test spA == spB
    end

    dividende = 10000
    b = big.((rand(Int, n) .% dividende) .// dividende)
    b = convert.(type, b)
    if isLU
        x, y = sol(A, b)
        x2, y2 = sol(spB, b)
        
        if isnan(tol)
            @test x == x2
        else
            println("norm=$(norm(x-x2))")
            res=Acp*x
            println("normres=$(norm(res-b))")
            @test isapprox(res, b, atol=tol)
            @test isapprox(x, x2, atol=tol)
        end
    end
end

function test_perf(n ,order, iscirc)
    t = getbspline(big(order),0).(1:order)
    A = topl(n, t, iscirc)
    ku = div(order,2)
    kl = order-1-ku
    @time  decLU(A)
    
    spA = LuSpline(A, ku, kl; iscirc=iscirc, isLU=true)
    @time spB = LuSpline(n, t; iscirc=iscirc, isLU=true)

    dividende = 10000
    b = big.((rand(Int, n) .% dividende) .// dividende)
        x, y = sol(A, b)
        x2, y2 = sol(spB, b)
        @test y == y2
        @test x == x2
    
end

function test_interface()
    bsp = B_SplineLU(25,105,big"0.")
    @test 105 == get_n(bsp)
    @test 25 == get_order(bsp)
    @test "B_SplineLU{BigFloat, true}" == get_type(bsp)
end
    

@testset "test decLU that is a tool for test" begin
    test_decLU(30)
end

@testset "test LuSpline" begin
Random.seed!(1899221)
test_splu(30, 9, true, false)
test_splu(30, 10, true, false)
test_splu(30, 9, false, false)
test_splu(30, 10, false, false)
test_splu(30, 9, false, true)
test_splu(30, 10, false, true)
test_splu(30, 9, true, true)
test_splu(31, 10, true, true)
test_splu(30, 9, true, true, type=BigFloat, tol=1e-70 )
test_splu(31, 10, true, true, type=BigFloat, tol=1e-70 )
test_splu(30, 11, false, true, type=BigFloat, tol=1e-70 )
test_splu(31, 10, false, true, type=BigFloat, tol=1e-70 )
test_splu(30, 3, true, false)
test_splu(30, 4, true, false)
test_splu(30, 3, false, false)
test_splu(30, 4, false, false)
test_splu(30, 3, false, true)
test_splu(30, 4, false, true)
test_splu(30, 3, true, true)
test_splu(31, 4, true, true)
test_splu(30, 3, true, true, type=BigFloat, tol=1e-70 )
test_splu(31, 4, true, true, type=BigFloat, tol=1e-70 )
test_splu(30, 5, false, true, type=BigFloat, tol=1e-70 )
test_splu(31, 4, false, true, type=BigFloat, tol=1e-70 )
test_interface()
end

@testset "test perf" begin
    test_splu(1000,9,true, true, type=Float64, tol=1e-12, perf=true)    
end
