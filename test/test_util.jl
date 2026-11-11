using SemiLagrangian:
    splititr,
    totuple,
    tupleshape,
    getextarray

@testset "split util" begin
    data = [
        (3, 15, [1:5, 6:10, 11:15])
        (3, 11, [1:4, 5:8, 9:11], 13)
        (5, 24, [1:5, 6:10, 11:15, 16:20, 21:24])
    ]
    for d in data
        @test d[3] == splititr(d[1], d[2])
    end
end

@testset "totuple" begin
    @test (3, 4, 1) == totuple([3, 4, 1])
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
    return @test size(tabext) == sz .+ decbeg .+ decend
end
@testset "ExtArray test" begin
    test_extarray(Float64, (20, 10), (3, 5), (2, 7))
    test_extarray(Float64, (4, 10, 7), (2, 3, 5), (1, 2, 4))
    test_extarray(Float64, (20,), (3,), (7,))
end
