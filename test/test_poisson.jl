
include("../src/poisson.jl")
include("../src/lagrange.jl")

using Test

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
    return totuple(UniformMesh.(t_deb,t_end,t_size; endpoint=false)), t_step
end


function test_poisson(T::DataType, isfft=true)
    
    t_debsp = T.([-1//1,-10//1,-3//1])
    t_endsp = T.([3//1, 6//1,5//1])
    t_szsp = (8, 4, 16)
    Nsp = 3
    t_debv = T.([-3//1, -9//1, 1//1, -1//1])
    t_endv = T.([ 1//1, 7//1, 5//1, 3//1])
    t_szv = (4, 8, 4, 4)
    base_dt = one(T)/80

    t_meshsp, t_stepsp = initmesh(t_debsp, t_endsp,t_szsp)
    t_meshv, t_stepv = initmesh(t_debv, t_endv, t_szv)

    interp = Lagrange(T,3)
    adv = Advection(t_meshsp, t_meshv, ntuple(x->interp,3), ntuple(x->interp,4), base_dt)

    tab = rand(T, t_szsp..., t_szv...)

    pc = PoissonConst(adv; isfftbig=isfft)

    pvar = PoissonVar(pc)

    advdata = AdvectionData(adv, tab, pvar)

    advdata.state_coef = 2
    advdata.state_dim = 1
    init!(advdata)
    rhoref = zeros(T, t_szsp)
    compute_charge!(rhoref, t_meshv, tab)

    @test rhoref == pvar.rho

    pfft = if isfft 
                PrepareFftBig(t_szsp, T; numdims=Nsp, dims=ntuple(x->x,Nsp))
    else
        missing
    end

    refelfield = compute_elfield(t_meshsp, rhoref, pfft)

    for i=1:Nsp
        @test refelfield[i] == pvar.t_elfield[i]
    end

end

@testset "Poisson Float64" begin
    test_poisson(Float64, true)
    test_poisson(Float64, false)
end
@testset "Poisson Float64" begin
    test_poisson(BigFloat)
end

