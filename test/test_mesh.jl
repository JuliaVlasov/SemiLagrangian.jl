include("../src/mesh.jl")

@testset "test selectdim for tuple" begin


    A = rand(5,5,5,5)
    B= selectdim(A,(1,2,3), (2,3,4))
    @test B == A[2, 3, 4, :]

    A = rand(5,5,5,5)
    B= selectdim(A,(1,3), (3,4))
    @test B == A[3, :, 4, :]

    A = rand(5,5,5,5)
    B= selectdim(A,(1,2), (4,3))
    @test B == A[4, 3, :, :]

end
