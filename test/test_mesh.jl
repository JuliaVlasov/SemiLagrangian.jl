
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
    t_debsp =[big"-1"//1,-10//1,-3//1]
    t_endsp = [big"3"//1, 6//1,5//1]
    t_szsp = [2, 4, 8]
    t_stepsp = (t_endsp - t_debsp) ./ t_szsp
    tt_meshsp = UniformMesh.(t_debsp, t_endsp, t_szsp; endpoint=false)
    t_meshsp = totuple(tt_meshsp)
    szsp=totuple(t_szsp)

    t_debv =[big"-3"//1,-9//1,1//1, -1//1]
    t_endv = [big"1"//1, 7//1,5//1,3//1]
    t_szv = [4, 8, 4, 2]
    t_stepv = (t_endv - t_debv) ./ t_szv
    tt_meshv = UniformMesh.(t_debv,t_endv,t_szv; endpoint=false)
    t_meshv = totuple(tt_meshv)
    szv=totuple(t_szv)

    fxv = rationalize.(BigInt, rand(Float64,(szsp...,szv...)))

    Nsp =length(szsp)
    Nv = length(szv)
    Nsum=Nsp+Nv
    T = Rational{BigInt}
    dx = prod(step, t_meshsp)
    dv = prod(step, t_meshv)
    sum_sp = Array{T,Nv}(undef,szv)
    sum_sp .= reshape(sum(fxv, dims = ntuple(x->x,Nsp)), szv )
    refres =  dx * dv * sum( dotprod(t_meshv) .^ 2 .* sum_sp)

    @test refres == compute_ke(t_meshsp, t_meshv, fxv)

end
# @testset "test compute_charge!" begin
#     t_debx =[big"-1"//1,-10//1,-3//1, -1//1]
#     t_endx = [big"3"//1, 6//1,5//1,1//1]
#     t_szx = [8, 4, 2, 4]
#     t_stepx = (t_endx - t_debx) ./ t_szx
#     tt_meshx = UniformMesh.(t_debx,t_endx,t_szx; endpoint=false)
#     t_meshx = totuple(tt_meshx)

#     t_debv =[big"-3"//1,-9//1,1//1, -1//1]
#     t_endv = [big"1"//1, 7//1,5//1,3//1]
#     t_szv = [4, 16, 4, 2]
#     t_stepv = (t_endv - t_debv) ./ t_szv
#     tt_meshv = UniformMesh.(t_debv,t_endv,t_szv; endpoint=false)
#     t_meshv = totuple(tt_meshv)

#     fvx = rationalize.(BigInt, rand(Float16,(totuple(t_szv)...,totuple(t_szx)...)))

#      N =4
#     N2 = 8
#     T = Rational{BigInt}
#     rho = Array{T,N}(undef,totuple(t_szx))
#     rhoref = Array{T,N}(undef,totuple(t_szx))
#     dv = prod(step, t_meshv)
#     tuplebegin = ntuple(x -> x, N)
#     rhoref .= dv * reshape(sum(fvx, dims = tuplebegin), totuple(t_szx))
#     rhoref .-= sum(rhoref)/length(rhoref)
#     compute_charge_vx!(rho, t_meshv, fvx)
    
#     @test rhoref == rho

#     fill!(rho, zero(rho[1]))

#     fxv = permutedims(fvx,vcat((N+1:2N),1:N))
#     compute_charge_xv!(rho, t_meshv, fxv)

#     @test rhoref == rho
# end



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
