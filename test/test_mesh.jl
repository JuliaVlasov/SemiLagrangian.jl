
using SemiLagrangian:
    UniformMesh, totuple, dotprod, dotprodother, length, step, prod, points, width, vec_k_fft, newval, traitmodbegin!, traitmodend!,
    meshtostd, stdtomesh

@testset "test tupleshape Mesh" begin
    t_deb = [-1 // 1, -10 // 1, -3 // 1, -1 // 1]
    t_end = [3 // 1, 6 // 1, 5 // 1, 1 // 1]
    t_sz = [20, 10, 8, 16]
    t_width = t_end - t_deb
    t_step = t_width ./ t_sz

    result = Array{Rational{Int},4}(undef, totuple(t_sz))
    ind = zeros(Int, 4)
    for ind[1] = 1:t_sz[1], ind[2] = 1:t_sz[2], ind[3] = 1:t_sz[3], ind[4] = 1:t_sz[4]
        res = 1 // 1
        for j = 1:4
            res *= t_deb[j] + (ind[j] - 1) * t_step[j]
        end
        result[ind[1], ind[2], ind[3], ind[4]] = res
    end

    tt_mesh = UniformMesh.(t_deb, t_end, t_sz)
    for (i, mesh) in enumerate(tt_mesh)
        @test t_width[i] == mesh.width == width(mesh)
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

function test_vec_k_fft(vbeg::T, vend::T, len) where {T}
    mid = div(len, 2)
    ref = circshift((-mid):(mid-1), mid)
    ref *= 2T(pi) / (vend - vbeg)
    @test ref == vec_k_fft(UniformMesh(vbeg, vend, len))
end


@testset "test vec_k_fft" begin
    test_vec_k_fft(-1.0, 1.0, 64)
    test_vec_k_fft(-big"1.0", big"1.0", 64)
end

function test_newval()
    @test newval(1.0, 1.2, 0.5, 1.0) ≈ 1.2
    @test newval(1.0, 0.2, 0.5, 1.0) ≈ 1.2
    @test newval(1.0, -0.8, 0.5, 1.0) ≈ 1.2
    @test newval(1.0, -10.8, 0.5, 1.0) ≈ 1.2
    @test newval(1.0, 3.3, 0.5, 1.0) ≈ 1.3
    @test newval(-2.4, 1.2, 0.5, 1.0) ≈ -2.8
    @test newval(-2.2, 1.2, 0.5, 1.0) ≈ -1.8    
    @test newval(0.5, 1.2, 0.5, 1.0) ≈ 0.2  
    @test newval(9.9,-9.7, 10.0, 20.0) ≈ 10.3
    @test newval(29.9,-9.7, 10.0, 20.0) ≈ 30.3
    @test newval(-9.7, 29.9, 10.0, 20.0) ≈ -10.1
    @test newval(0.0,-3.2,0.5,1.0) ≈ -0.2
    @test newval(0.0,-2.8,0.5,1.0) ≈ 0.2
    @test newval(0.0,-2.2,0.5,1.0) ≈ -0.2
    @test newval(0.0,-1.8,0.5,1.0) ≈ 0.2
    @test newval(0.0,-1.2,0.5,1.0) ≈ -0.2
    @test newval(0.0,-0.8,0.5,1.0) ≈ 0.2
    @test newval(0.0,-0.2,0.5,1.0) ≈ -0.2
    @test newval(0.0,0.2,0.5,1.0) ≈ 0.2

end

function test_traitmod(c1,c2,c3)

    tab = [ mod(i*c1+j*c2, 1.0)+c3 for i=1:10,j=1:10]

    mesh = UniformMesh(0.0,1.0,64)

    res = copy(tab)
    ret = traitmodbegin!(mesh,res)
    @test ret
    res2 = copy(res)
    ret2 = traitmodbegin!(mesh,res2)
    @test !ret2
    @test res2 == res

    traitmodend!(mesh, res)
    result1 = mod.(tab,1.0)
    result2 = mod.(res,1.0)
    @test result1 == result2 == res

    stdtab = meshtostd(mesh, tab)
    lg=float(length(mesh))
    stdres = copy(stdtab)
    ret = traitmodbegin!(lg, stdres)
    @test ret
    stdres2 = copy(stdres)
    ret2 = traitmodbegin!(lg, stdres2)
    @test !ret2
    @test stdres2 == stdres

    traitmodend!(lg, stdres)
    result1 = mod.(stdtab,lg)
    result2 = mod.(stdres,lg)
    @test result1 == result2 == stdres

    resres = stdtomesh(mesh,stdres)

    @test resres == res

end
@testset "test newval" begin
    test_newval()
end
@testset "test traimtmod" begin
    test_traitmod(0.127615, 0.0871745, -24.89897)
    test_traitmod(0.127615, -0.0871745, -24.89897)
    test_traitmod(-0.345, 0.1871745, 24.89897)
end
