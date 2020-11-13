include("../src/advection.jl")
include("../src/lagrange.jl")

using Test


function initmesh(t_deb, t_end, t_size)
    t_step = (t_end - t_deb) ./ t_size
    return totuple(UniformMesh.(t_deb,t_end,t_size; endpoint=false)), t_step
end

@testset "test get_kl_ku" begin
    kl, ku = get_kl_ku(5)
    @test kl == 2 && ku == 2

    kl, ku = get_kl_ku(6)
    @test kl == 2 && ku == 3
end

@testset "test Advection Float" begin
T = Float64
t_debsp = T.([-1,-10,-3,-1])
t_endsp = T.([3, 6, 5, 1])
t_szsp = (16, 8, 32, 8)
t_debv = T.([-3,-9,1,-1])
t_endv = T.([1, 7, 1, 3])
t_szv = (4, 8, 4, 4)

t_meshsp, t_stepsp = initmesh(t_debsp, t_endsp,t_szsp)
t_meshv, t_stepv = initmesh(t_debv, t_endv, t_szv)

interp = Lagrange(T,3)
adv = Advection(t_meshsp, t_meshv, ntuple(x->interp,4),ntuple(x->interp,4),one(T)/80)

sref = (t_szsp..., t_szv...)
@test sref == sizeall(adv)

refitr = ntuple(x-> 1:sref[x],size(sref,1))
@test refitr == sizeitr(adv)

tab = rand(T, sizeall(adv))
end
