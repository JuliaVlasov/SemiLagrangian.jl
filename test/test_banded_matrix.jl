@testset "Banded Matrix solver" begin

    n, kl, ku = 9, 2, 3

    matrix = BandedMatrix( n, kl, ku )

    for i in 1:n

        set_element!(matrix, i, min(i+3,n), 4.)
        set_element!(matrix, i, min(i+2,n), 3.)
        set_element!(matrix, i, min(i+1,n), 2.)
        set_element!(matrix, i, max(i-2,1), 3.)
        set_element!(matrix, i, max(i-1,1), 2.)
        set_element!(matrix, i, i, 1.)

    end
    

    factorize!( matrix )

    b = ones(Float64, (n,3))
    b[1:end,2] .= [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0]
    b[1:end,3] .= 1.0:1.0:9.0

    solve!( matrix, b )

    bexact = [ 0.629   1.149   4.713 ;  
              -0.025   1.870  -0.726 ; 
               0.523   0.615   2.946 ; 
              -0.286  -1.434  -2.774 ; 
              -0.104  -0.524  -1.067 ; 
              -0.118  -0.591  -0.970 ; 
               0.220  -0.896   2.027 ;  
              -0.079   1.604  -0.340 ;  
               0.496   0.480   3.599 ]

    @test isapprox(sum(abs.(b.-bexact))/3n, 0.0, atol=1e-3)

end
