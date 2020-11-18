
include("../src/mesh.jl")

using Test
@testset "test tupleshape" begin

    @test (3,4,1) == totuple([3,4,1])
    @test [5,1,9,2] == tovector((5,1,9,2))
    
    @test (1,1,71,1) == tupleshape(3,4,71)
    @test reshape((-5//1):1//10:49//10,(1,100,1,1)) == tupleshape(2,4,UniformMesh(-5//1,5//1,100, endpoint=false))
    t_deb =[-1//1,-10//1,-3//1, -1//1]
    t_end = [3//1, 6//1,5//1,1//1]
    t_sz = [20, 10, 8, 16]
    t_step = (t_end - t_deb) ./ t_sz

    result = Array{Rational{Int},4}(undef, totuple(t_sz))
    ind = zeros(Int,4)
    for ind[1]=1:t_sz[1], ind[2]=1:t_sz[2], ind[3]=1:t_sz[3], ind[4]=1:t_sz[4]
        res = 1//1
        for j = 1:4
            res *= t_deb[j] + (ind[j]-1)*t_step[j]
        end
        result[ind[1],ind[2],ind[3],ind[4]] = res  
    end
#    tt_mesh = [UniformMesh(t_deb[i],t_end[i],t_sz[i]; endpoint=false) for i=1:4]
    
    tt_mesh = UniformMesh.(t_deb,t_end,t_sz; endpoint=false)
    t_mesh = totuple(tt_mesh)
#     println("typeof(t_mesh)=$(typeof(t_mesh))")
    @time @test result == dotprod(t_mesh)
    @time @test result == dotprodother(t_mesh)

    @test prod(t_step) == prod(step, t_mesh)

    @test totuple(t_sz) == length.(t_mesh)
    @test t_sz == length.(tt_mesh)
    @test totuple(t_step) == step.(t_mesh)
    @test t_step == step.(tt_mesh)
    
end

function test_vec_k_fft(vbeg::T, vend::T, len) where{T}
    mid=div(len,2)
    ref = circshift((-mid):(mid-1), mid)
    ref *= 2T(pi)/(vend-vbeg)
    @test ref == vec_k_fft(UniformMesh(vbeg, vend, len; endpoint=false))
end 


@testset "test vec_k_fft" begin
    test_vec_k_fft(-1.0, 1.0, 64)
    test_vec_k_fft(-big"1.0", big"1.0", 64)
end


