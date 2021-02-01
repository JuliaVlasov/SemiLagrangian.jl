


using SemiLagrangian: splititr, splitvec, transperm, totuple, tovector, tupleshape

@testset "split util" begin
    data = [ 
        (3, 15, [1:5, 6:10, 11:15])
        (3, 11, [1:4, 5:8, 9:11], 13)
        (5, 24, [1:5, 6:10, 11:15, 16:20, 21:24]) ]
    for d in data
        @test d[3] == splititr(d[1], d[2])
    end
    data2 = [
        (47, 3, 15, [47:51, 52:56, 57:61])
        (34, 4, 5, [34:35, 36:36, 37:37, 38:38])
        (23, 7, 13, [23:24, 25:26, 27:28, 29:30, 31:32, 33:34, 35:35])
    ]
    for d in data2
        @test d[4] == splitvec(d[2], collect(d[1]:d[1]+d[3]-1))
    end
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

@testset "transperm" begin
    @test transperm(1,2,5) == [2,1,3,4,5]
    @test transperm(1,2,5) == [2,1,3,4,5]
    @test transperm(3,7,7) == [1,2,7,4,5,6,3]
    @test transperm(7,3,7) == [1,2,7,4,5,6,3]
end

@testset "totuple tovector" begin
    @test (3,4,1) == totuple([3,4,1])
    @test [5,1,9,2] == tovector((5,1,9,2))
end

@testset "tupleshape" begin
    @test (1,1,71,1) == tupleshape(3,4,71)
    @test (53, 1) == tupleshape(1,2,53)
    @test (1, 1, 67) == tupleshape(3,3,67)
end

