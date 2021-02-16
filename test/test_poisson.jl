

using DoubleFloats

using SemiLagrangian: UniformMesh, compute_ee, compute_ke, PoissonConst, PoissonVar, initcoef!, compute_charge!

function compute_elfield(
    t_mesh_x::NTuple{N,UniformMesh{T}},
    rho::Array{T, N},
    pfft
) where {T <: AbstractFloat, N}

    fct_k(v)= im/sum(v.^2)

    v_k = vec_k_fft.(t_mesh_x)
    sz = length.(t_mesh_x)

    buf = fftgenall(pfft, rho) .* fct_k.(collect(Iterators.product(v_k...)))
    buf[1] = 0im

    return ntuple( 
    x -> real(ifftgenall(pfft, reshape(v_k[x],tupleshape(x,N,sz[x])) .* buf )),
    N
)
end

function initmesh(t_deb, t_end, t_size)
    t_step = (t_end - t_deb) ./ t_size
    return totuple(UniformMesh.(t_deb,t_end,t_size)), t_step
end


function test_poisson(T::DataType, isfft=true)
    
    t_debsp = T.([-1//1,-10//1,-3//1])
    t_endsp = T.([3//1, 6//1,5//1])
    t_szsp = (8, 4, 16)
    Nsp = 3
    t_debv = T.([-3//1, -9//1, 1//1])
    t_endv = T.([ 1//1, 7//1, 5//1])
    t_szv = (4, 8, 4)
    base_dt = one(T)/80

    t_meshsp, t_stepsp = initmesh(t_debsp, t_endsp,t_szsp)
    t_meshv, t_stepv = initmesh(t_debv, t_endv, t_szv)

    interp = Lagrange(3, T)
    adv = Advection(t_meshsp, t_meshv, ntuple(x->interp,3), ntuple(x->interp,3), base_dt)

    tab = rand(T, t_szsp..., t_szv...)

    pc = PoissonConst(adv; isfftbig=isfft)

#    println("t_perms=$(pc.t_perms)")
    @test pc.t_perms == (
    [1, 2, 3, 6, 5, 4], [2, 1, 3, 4, 6, 5], [3, 2, 1, 4, 5, 6], # space states
    [4, 5, 6, 1, 2, 3], [5, 4, 6, 1, 2, 3], [6, 5, 4, 1, 2, 3]  # velocity states
)

    pvar = PoissonVar(pc)

    advdata = AdvectionData(adv, tab, pvar)

    advdata.state_coef = 2
    advdata.state_dim = 1
    initcoef!(pvar, advdata)
    rhoref = zeros(T, t_szsp)
    compute_charge!(rhoref, t_meshv, tab)

    @test rhoref == pvar.rho

    pfft = if isfft 
            PrepareFftBig(t_szsp, T; numdims=Nsp, dims=ntuple(x->x,Nsp))
    else
        missing
    end

    refelfield = compute_elfield(t_meshsp, rhoref, pfft)

#    prec = (T==BigFloat) ? 1e-70 : 1e-15

    for i=1:Nsp
#        @test isapprox(refelfield[i], pvar.t_elfield[i], atol=prec, rtol=prec)
        @test isapprox(refelfield[i], pvar.t_elfield[i])
    end

#    @test isapprox(compute_ee(t_meshsp, refelfield), compute_ee(advdata), atol=prec, rtol=prec)
    @test isapprox(compute_ee(t_meshsp, refelfield), compute_ee(advdata))

end


@testset "test compute_ee" begin
    t_deb =[big"-1"//1,-10//1,-3//1, -1//1]
    t_end = [big"3"//1, 6//1,5//1,1//1]
    t_sz = [4, 4, 8, 16]
    t_step = (t_end - t_deb) ./ t_sz
    tt_mesh = UniformMesh.(t_deb,t_end,t_sz)
    t_mesh = totuple(tt_mesh)
    N=4
    t = ntuple( x->rationalize.(BigInt, rand(Float64,totuple(t_sz))), N)
    
    resref = sum( map( x->prod(step, t_mesh) * sum(x.^2), t ) )

    @test resref == compute_ee(t_mesh, t)

end

@testset "Poisson Float64" begin
    test_poisson(Float64, true)
    test_poisson(Float64, false)
    test_itr(Float64)
end


@testset "Poisson BigFloat" begin
    test_poisson(BigFloat)
    test_itr(BigFloat)
end

@testset "Poisson Double64" begin
    test_poisson(Double64)
    test_itr(Double64)
end

