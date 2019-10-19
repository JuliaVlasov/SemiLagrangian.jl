export Geometry

struct Geometry

    n1 :: Int64
    n2 :: Int64

    x1min  :: Float64
    x2min  :: Float64

    delta1 :: Float64
    delta2 :: Float64
    
    x1grid :: Array{Float64,1}
    x2grid :: Array{Float64,1}

    bc :: Symbol

"""
initialize geometry from origin cordinates, 
space steps and number of points
"""
    function Geometry(n1, n2, x1min, x2min, delta1, delta2)

        x1grid = collect(range(x1min, length=n1, step=delta1))
        x2grid = collect(range(x2min, length=n2, step=delta2))

        new(n1, n2, x1min, x2min, delta1, delta2, x1grid, x2grid, :none) 

    end 

"""
initialize geometry from origin and end coordinates and 
number of points
"""
    function Geometry(n1, n2, x1, x2, bc::Symbol)

	x1min, x1max = x1
	x2min, x2max = x2
   
        if ( (bc == :perxy) || (bc == :perx) ) 
           x1grid = collect(range(x1min, stop=x1max, length=n1+1)[1:end-1])
        else
           x1grid = collect(range(x1min, stop=x1max, length=n1))
        end

        if ( (bc == :perxy) || (bc == :pery) ) 
           x2grid = collect(range(x2min, stop=x2max, length=n2+1)[1:end-1])
        else
           x2grid = collect(range(x2min, stop=x2max, length=n2))
        end

        delta1 = x1grid[2] - x1grid[1]
        delta2 = x2grid[2] - x2grid[1]

        new(n1, n2, x1min, x2min, delta1, delta2, x1grid, x2grid, bc) 

    end 

end

