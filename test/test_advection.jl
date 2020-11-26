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

@testset "test tools" begin

    v = collect(1:53)
    t = splitvec(5,v)
    @test t[1] == collect(1:11)
    @test t[2] == collect(12:22)
    @test t[3] == collect(23:33)
    @test t[4] == collect(34:43)
    @test t[5] == collect(44:53)

    @test transperm(1,2,5) == [2,1,3,4,5]
    @test transperm(4, 2, 7) == [1, 4, 3, 2, 5, 6, 7]


end

function test_adv(T::DataType)
    t_debsp = T.([-1,-10,-3])
    t_endsp = T.([3, 6, 5])
    t_szsp = (16, 8, 32)
    t_debv = T.([-3,-9,1,-1])
    t_endv = T.([1, 7, 1, 3])
    t_szv = (4, 8, 4, 4)
    base_dt = one(T)/80

    Nsum = 7

    t_meshsp, t_stepsp = initmesh(t_debsp, t_endsp,t_szsp)
    t_meshv, t_stepv = initmesh(t_debv, t_endv, t_szv)

    interp = Lagrange(T,3)
    adv = Advection(
    t_meshsp, t_meshv, 
    ntuple(x->Lagrange(T,3),3), ntuple(x->Lagrange(T,3),4), 
    base_dt
)

    sref = (t_szsp..., t_szv...)
    @test sref == sizeall(adv)

    refitr = ntuple(x-> 1:sref[x],size(sref,1))
    @test refitr == sizeitr(adv)

    tab = rand(T, sizeall(adv))

    advd = AdvectionData(adv, tab, "missing")


    @test compute_ke(t_meshsp, t_meshv, tab) == compute_ke(advd)

    @test advd.state_coef == 1
    @test advd.state_dim == 1

    @test advd.parext == getext(advd)
    @test advd.parext == "missing"
    @test advd.data == getdata(advd)
    @test base_dt * adv.tab_coef[2] == getcur_t(adv, 2)
    @test 1 == getstate_dim(advd)

    @test !isvelocitystate(advd)
    @test getcur_t(advd) == base_dt * advd.adv.tab_coef[1]


    @test addcolon(3,(1,2,3,4,5)) == (1, 2, :, 3, 4, 5)



    t_coef =[1,1,1,2,2,2,2,3,3,3,1]
    t_dim = [1,2,3,1,2,3,4,1,2,3,1]
    t_indice = [1,2,3,4,5,6,7,1,2,3,1]
    t_v=[false,false,false,true,true,true,true,false,false,false,false]
    t_result=[true,true,true,true,true,true,true,true,true,false,true]


    
 
    for i=1:length(t_coef) 
        @test advd.state_coef == t_coef[i] 
        @test advd.state_dim == t_dim[i] 
        @test t_indice[i] == _getcurrentindice(advd) 
        @test isvelocitystate(advd) == t_v[i]
        @test advd.state_dim == getstate_dim(advd)
        @test getcur_t(advd) == base_dt * adv.tab_coef[t_coef[i]]
        @test getbufslgn(advd) == advd.t_buf[t_indice[i]]
        t = isvelocitystate(advd) ? adv.t_interp_v : adv.t_interp_sp
        @test t[t_dim[i]] == getinterp(advd)

        x = t_indice[i]
        @time @test addcolon.(x, Iterators.product(refitr[vcat(1:(x-1),(x+1):Nsum)]...)) == getitr(advd)
        
        ret = nextstate!(advd)
        @test ret == t_result[i]
    end
 
    

end

@testset "test Advection Float" begin

    test_adv(Float64)

end
@testset "test Advection BigFloat" begin

    test_adv(BigFloat)

end
function test_ke(T::DataType)
    t_debsp = T.([-1//1,-10//1,-3//1])
    t_endsp = T.([3//1, 6//1,5//1])
    t_szsp = [2, 4, 8]
    t_stepsp = (t_endsp - t_debsp) ./ t_szsp
    tt_meshsp = UniformMesh.(t_debsp, t_endsp, t_szsp; endpoint=false)
    t_meshsp = totuple(tt_meshsp)
    szsp=totuple(t_szsp)

    t_debv = T.([-3//1,-9//1,1//1, -1//1])
    t_endv = T.([1//1, 7//1,5//1,3//1])
    t_szv = [4, 8, 4, 2]
    t_stepv = (t_endv - t_debv) ./ t_szv
    tt_meshv = UniformMesh.(t_debv,t_endv,t_szv; endpoint=false)
    t_meshv = totuple(tt_meshv)
    szv=totuple(t_szv)

    fxv = if T <: Rational
        rationalize.(BigInt, rand(Float64,(szsp...,szv...)))
    else
        rand(T, (szsp...,szv...))
    end

    Nsp =length(szsp)
    Nv = length(szv)
    Nsum=Nsp+Nv
    dx = prod(step, t_meshsp)
    dv = prod(step, t_meshv)
    sum_sp = Array{T,Nv}(undef,szv)
    sum_sp .= reshape(sum(fxv, dims = ntuple(x->x,Nsp)), szv )
    refres =  (dx * dv ) * sum( dotprod(t_meshv) .^ 2 .* sum_sp)

    @test refres == compute_ke(t_meshsp, t_meshv, fxv)
end

@testset "test compute_ke" begin
    test_ke(Rational{BigInt})
    test_ke(Float64)
    test_ke(BigFloat)
end
