
include("../src/poisson.jl")

using Test
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
    T = BigFloat
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
    Base.GC.gc(false)
    @time e = compute_elfield(t_meshx, rho, pfft )
    Base.GC.gc(false)
    @time e_other = compute_elfieldother(t_meshx, rho, pfft )

    e_mem = ntuple(x-> zeros(Complex{T},totuple(t_szx)...), N)
    Base.GC.gc(false)
    @time compute_elfield!(e_mem, t_meshx, rho, pfft )
    e_memother = ntuple(x-> zeros(Complex{T},totuple(t_szx)...), N)
    Base.GC.gc(false)
    @time compute_elfieldother!(e_memother, t_meshx, rho, pfft )
    e_memother2 = ntuple(x-> zeros(Complex{T},totuple(t_szx)...), N)
    Base.GC.gc(false)
    @time compute_elfieldother2!(e_memother2, t_meshx, rho, pfft )

    for i=1:N
        @test isapprox(e_other[i], e[i], rtol=1e-40, atol=1e-40)
        @test isapprox(e_other[i], e_mem[i], rtol=1e-40, atol=1e-40)
        @test isapprox(e_other[i], e_memother[i], rtol=1e-40, atol=1e-40)
        @test isapprox(e_other[i], e_memother2[i], rtol=1e-40, atol=1e-40)
    end

end
