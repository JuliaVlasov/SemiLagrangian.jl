
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

@testset "test compute_ee" begin
    t_deb =[big"-1"//1,-10//1,-3//1, -1//1]
    t_end = [big"3"//1, 6//1,5//1,1//1]
    t_sz = [4, 4, 8, 16]
    t_step = (t_end - t_deb) ./ t_sz
    tt_mesh = UniformMesh.(t_deb,t_end,t_sz; endpoint=false)
    t_mesh = totuple(tt_mesh)
    t = rationalize.(BigInt, rand(Float16,totuple(t_sz)))

    resref = prod(step, t_mesh) * sum(t.^2)
    @test resref == compute_ee(t_mesh, t)
end


@testset "test compute_ke" begin
    t_debx =[big"-1"//1,-10//1,-3//1, -1//1]
    t_endx = [big"3"//1, 6//1,5//1,1//1]
    t_szx = [2, 4, 8, 4]
    t_stepx = (t_endx - t_debx) ./ t_szx
    tt_meshx = UniformMesh.(t_debx,t_endx,t_szx; endpoint=false)
    t_meshx = totuple(tt_meshx)

    t_debv =[big"-3"//1,-9//1,1//1, -1//1]
    t_endv = [big"1"//1, 7//1,5//1,3//1]
    t_szv = [4, 8, 4, 2]
    t_stepv = (t_endv - t_debv) ./ t_szv
    tt_meshv = UniformMesh.(t_debv,t_endv,t_szv; endpoint=false)
    t_meshv = totuple(tt_meshv)

    fvx = rationalize.(BigInt, rand(Float16,(totuple(t_szv)...,totuple(t_szx)...)))

    N =4
    N2 = 8
    T = Rational{BigInt}
    dx = prod(step, t_meshx)
    dv = prod(step, t_meshv)
    sum_x = Array{T,N}(undef,size(fvx)[1:N])
    sum_x .= reshape(sum(fvx, dims = totuple((N+1):N2)), size(fvx)[1:N] )
    refres =  dx * dv * sum( dotprod(t_meshv) .^ 2 .* sum_x)

    @test refres == compute_ke(t_meshv, t_meshx, fvx)

end
@testset "test compute_charge!" begin
    t_debx =[big"-1"//1,-10//1,-3//1, -1//1]
    t_endx = [big"3"//1, 6//1,5//1,1//1]
    t_szx = [8, 4, 2, 4]
    t_stepx = (t_endx - t_debx) ./ t_szx
    tt_meshx = UniformMesh.(t_debx,t_endx,t_szx; endpoint=false)
    t_meshx = totuple(tt_meshx)

    t_debv =[big"-3"//1,-9//1,1//1, -1//1]
    t_endv = [big"1"//1, 7//1,5//1,3//1]
    t_szv = [4, 16, 4, 2]
    t_stepv = (t_endv - t_debv) ./ t_szv
    tt_meshv = UniformMesh.(t_debv,t_endv,t_szv; endpoint=false)
    t_meshv = totuple(tt_meshv)

    fvx = rationalize.(BigInt, rand(Float16,(totuple(t_szv)...,totuple(t_szx)...)))

     N =4
    N2 = 8
    T = Rational{BigInt}
    rho = Array{T,N}(undef,totuple(t_szx))
    rhoref = Array{T,N}(undef,totuple(t_szx))
    dv = prod(step, t_meshv)
    tuplebegin = ntuple(x -> x, N)
    rhoref .= dv * reshape(sum(fvx, dims = tuplebegin), totuple(t_szx))
    rhoref .-= sum(rhoref)/length(rhoref)
    compute_charge_vx!(rho, t_meshv, fvx)
    
    @test rhoref == rho

    fill!(rho, zero(rho[1]))

    fxv = permutedims(fvx,vcat((N+1:2N),1:N))
    compute_charge_xv!(rho, t_meshv, fxv)

    @test rhoref == rho





end

@testset "test compute_elfield!" begin
 
    t_debx = Float64.([big"-1"//1,-10//1,-3//1, -1//1])
    t_endx = Float64.([big"3"//1, 6//1,5//1,1//1])
    t_szx = [16, 32, 128, 64]
    t_stepx = (t_endx - t_debx) ./ t_szx
    tt_meshx = UniformMesh.(t_debx,t_endx,t_szx; endpoint=false)
    t_meshx = totuple(tt_meshx)

    # t_debv =Float64.([big"-3"//1,-9//1,1//1, -1//1])
    # t_endv = Float64.([big"1"//1, 7//1,5//1,3//1])
    # t_szv = [4, 8, 4, 4]
    # t_stepv = (t_endv - t_debv) ./ t_szv
    # tt_meshv = UniformMesh.(t_debv,t_endv,t_szv; endpoint=false)
    # t_meshv = totuple(tt_meshv)

    rho = rand(totuple(t_szx)...)
    rho .-= sum(rho)/length(rho)
    println("size(rho)=$(size(rho))")

    fct_k(ind,v)= v[ind] == 0 ? 0im : im*v[ind]/sum(v.^2)

    v_k = vec_k_fft.(t_meshx)

    pfft = PrepareFftBig(length.(t_meshx),Float64;numdims=4,dims=ntuple(x->x,4))

#     buf = fftgen(pfft, rho)

#     array_k = collect(Iterators.product(v_k...))
    
#     println("size(buf)=$(size(buf)) size(array_k)=$(size(array_k))")

     N = 4

#     e_ref = ntuple( 
#     x -> real(ifftgen(pfft, fct_k.(x, array_k) .* buf )),
#     N
# )
    @time e = compute_elfield(t_meshx, rho, pfft )
    @time e_other = compute_elfieldother(t_meshx, rho, pfft )
    @time e2 = compute_elfield(t_meshx, rho, missing )
    @time e2_other = compute_elfieldother(t_meshx, rho, missing )

    for i=1:N
        @test isapprox(e_other[i], e[i], rtol=1e-10, atol=1e-10)
        @test isapprox(e2_other[i], e2[i], rtol=1e-10, atol=1e-10)
        @test isapprox(e2[i], e[i], rtol=1e-10, atol=1e-10)
    end

end
@testset "test compute_elfield! BigFloat" begin
 
    t_debx = float.([big"-1"//1,-10//1,-3//1, -1//1])
    t_endx = float.([big"3"//1, 6//1,5//1,1//1])
    t_szx = [16, 8, 32, 8]
    t_stepx = (t_endx - t_debx) ./ t_szx
    tt_meshx = UniformMesh.(t_debx,t_endx,t_szx; endpoint=false)
    t_meshx = totuple(tt_meshx)

    # t_debv =Float64.([big"-3"//1,-9//1,1//1, -1//1])
    # t_endv = Float64.([big"1"//1, 7//1,5//1,3//1])
    # t_szv = [4, 8, 4, 4]
    # t_stepv = (t_endv - t_debv) ./ t_szv
    # tt_meshv = UniformMesh.(t_debv,t_endv,t_szv; endpoint=false)
    # t_meshv = totuple(tt_meshv)

    rho = rand(BigFloat, totuple(t_szx)...)
    rho .-= sum(rho)/length(rho)
    println("size(rho)=$(size(rho))")

    fct_k(ind,v)= v[ind] == 0 ? 0im : im*v[ind]/sum(v.^2)

    v_k = vec_k_fft.(t_meshx)

    pfft = PrepareFftBig(length.(t_meshx),BigFloat;numdims=4,dims=ntuple(x->x,4))

#     buf = fftgen(pfft, rho)

#     array_k = collect(Iterators.product(v_k...))
    
#     println("size(buf)=$(size(buf)) size(array_k)=$(size(array_k))")

     N = 4

#     e_ref = ntuple( 
#     x -> real(ifftgen(pfft, fct_k.(x, array_k) .* buf )),
#     N
# )
    @time e = compute_elfield(t_meshx, rho, pfft )
    @time e_other = compute_elfieldother(t_meshx, rho, pfft )

    for i=1:N
        @test isapprox(e_other[i], e[i], rtol=1e-40, atol=1e-40)
    end

end

# @testset "test selectdim for tuple" begin


#     A = rand(5,5,5,5)
#     B= selectdim(A,(1,2,3), (2,3,4))
#     @test B == A[2, 3, 4, :]

#     A = rand(5,5,5,5)
#     B= selectdim(A,(1,3), (3,4))
#     @test B == A[3, :, 4, :]

#     A = rand(5,5,5,5)
#     B= selectdim(A,(1,2), (4,3))
#     @test B == A[4, 3, :, :]

# end
