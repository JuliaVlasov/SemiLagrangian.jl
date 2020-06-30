

@testset " Get cell number and offset " begin

     n       = 10
     degree  = 5
     bspline = Spline1D( n, degree, -1.0, 1.0  )
     tau     = -1 .+ 2 .* rand(n)

     err = 0.0
     for x in tau 
         cell, offset = get_cell_and_offset( bspline, x )
         err += abs( x - (bspline.start + (cell-1+offset) * bspline.step))
     end
     @test err ≈ 0.0 atol = 1e-15

end
