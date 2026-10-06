
using SemiLagrangian:
    transposition,
    totuple,
    tovector,
    tupleshape,
    getextarray,
    CircEdge,
    AbstractInterpolation,
    interpolatemod!

@testset "test tools" begin
    @test transposition(1, 2, 5) == [2, 1, 3, 4, 5]
    @test transposition(4, 2, 7) == [1, 4, 3, 2, 5, 6, 7]
end

@testset "transposition" begin
    @test transposition(1, 2, 5) == [2, 1, 3, 4, 5]
    @test transposition(1, 2, 5) == [2, 1, 3, 4, 5]
    @test transposition(3, 7, 7) == [1, 2, 7, 4, 5, 6, 3]
    @test transposition(7, 3, 7) == [1, 2, 7, 4, 5, 6, 3]
end

@testset "totuple tovector" begin
    @test (3, 4, 1) == totuple([3, 4, 1])
    @test [5, 1, 9, 2] == tovector((5, 1, 9, 2))
end

@testset "tupleshape" begin
    @test (1, 1, 71, 1) == tupleshape(3, 4, 71)
    @test (53, 1) == tupleshape(1, 2, 53)
    @test (1, 1, 67) == tupleshape(3, 3, 67)
end

function test_extarray(T, sz, decbeg, decend)
    tabor = rand(T, sz)
    tabext = getextarray(tabor, decbeg, decend)
    for ind in CartesianIndices(tabext)
        @test tabext[ind] == tabor[CartesianIndex(mod.(ind.I .- decbeg .- 1, sz) .+ 1)]
    end
    @test size(tabext) == sz .+ decbeg .+ decend
end
@testset "ExtArray test" begin
    test_extarray(Float64, (20, 10), (3, 5), (2, 7))
    test_extarray(Float64, (4, 10, 7), (2, 3, 5), (1, 2, 4))
    test_extarray(Float64, (20,), (3,), (7,))
end
