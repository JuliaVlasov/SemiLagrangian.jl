
using SemiLagrangian: UniformMesh, totuple, dotprod, dotprodother, length, step, prod, points, vec_k_fft 

@testset "test tupleshape Mesh" begin 
    t_deb =[-1//1,-10//1,-3//1, -1//1]
    t_end = [3//1, 6//1,5//1,1//1]
    t_sz = [20, 10, 8, 16]
    t_width = t_end - t_deb
    t_step = t_width ./ t_sz

    result = Array{Rational{Int},4}(undef, totuple(t_sz))
    ind = zeros(Int,4)
    for ind[1]=1:t_sz[1], ind[2]=1:t_sz[2], ind[3]=1:t_sz[3], ind[4]=1:t_sz[4]
        res = 1//1
        for j = 1:4
            res *= t_deb[j] + (ind[j]-1)*t_step[j]
        end
        result[ind[1],ind[2],ind[3],ind[4]] = res  
    end
    
    tt_mesh = UniformMesh.(t_deb,t_end,t_sz)
    for (i, mesh) in enumerate(tt_mesh)
        @test t_width[i] == mesh.width
    end
    t_mesh = totuple(tt_mesh)
#     println("typeof(t_mesh)=$(typeof(t_mesh))")
    @time @test result == dotprod(points.(t_mesh))
    @time @test result == dotprodother(points.(t_mesh))

    t_v = points.(t_mesh)
    @time @test result == dotprod(t_v)
    @time @test result == dotprodother(t_v)


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
    @test ref == vec_k_fft(UniformMesh(vbeg, vend, len))
end 


@testset "test vec_k_fft" begin
    test_vec_k_fft(-1.0, 1.0, 64)
    test_vec_k_fft(-big"1.0", big"1.0", 64)
end

