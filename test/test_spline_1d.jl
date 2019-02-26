using IterTools

mutable struct Spline1d

    degree :: Int64
    ncells :: Int64
    nbasis :: Int64
    start  :: Float64
    stop   :: Float64
    step   :: Float64
    bcoef  :: Vector{Float64}

    function Spline1d( n :: Int64, p :: Int64, xmin, xmax  )
      
        degree   = p
        ncells   = n
        nbasis   = ncells+degree
        start    = xmin
        stop     = xmax
        step     = (xmax-xmin) / ncells
        bcoef    = zeros( 1:n+p )

    end

end

function eval( self, x )

    values = zeros(self.degree+1)

    cell, offset = get_cell_and_offset( self, x )

    jmin = icell

    values[1] = 1.0
    for j = 1:self.degree
        xx     = -offset
        saved  = 0.0
        for r = 1:j
            xx        = xx + 1.0
            temp      = values[r] / j
            values[r] = saved + xx * temp
            saved     = (j - 1 - xx) * temp
         end do
         values[j] = saved
      end do

    jmax = jmin + self.degree

    self.bcoef[jmin:jmax] * values

end 

function get_cell_and_offset( self, x )

    offset = (x - self.start) / self.step  
    cell   = trunc( Int64, x_offset )
    offset = offset - cell
    cell   = min( cell+1, self.stop)

    cell, offset

end 

mutable struct SplineInterpolator1d

    bspl     :: Spline1d
    tau      :: Vector{Float64}
    matrix   :: Array{Float64, 2}
    bc_start :: Symbol
    bc_stop  :: Symbol

    function SplineInterpolator1d( bspl     :: Spline1D, 
                                   bc_start :: Symbol, 
                                   bc_stop  :: Symbol )

        if bc_start == :hermite
            nbc_start = degree÷2
        else
            nbc_start = 0
        end

        if bc_stop == :hermite
            nbc_stop = degree÷2
        else
            nbc_stop = 0
        end

        ntau = nbasis - nbc_start - nbc_stop

        tau = zeros(ntau)

        if ( bc_xmin == :periodic )

            if bspl.degree & 1 > 0  
                tau .= [xmin + (i-1.0)*dx for i=1:ntau]
            else
                tau .= [xmin + (i-0.5)*dx for i=1:ntau]
            end

        else

            # Non-periodic case: create array of temporary knots (integer shifts only)
            # in order to compute interpolation points using Greville-style averaging:
            # tau(i) = xmin + average(knots(i+1-degree:i)) * dx
            iknots = zeros(Float64, 2-degree:ntau)

            # Additional knots near x=xmin
            r = 2-degree
            s = -nbc_xmin
            if (bc_start == :greville) iknots[r:s] = 0 end
            if (bc_stop  == :hermite ) iknots[r:s] = [i for i=r-s-1:-1]

            # Knots inside the domain
            r = -nbc_xmin+1
            s = -nbc_xmin+1+ncells
            iknots[r:s] = [i for i=0:ncells]

            # Additional knots near x=xmax
            r = -nbc_xmin+1+ncells+1
            s = ntau
            if bc_stop = :greville  
                iknots[r:s] .= ncells 
            end
            if bc_stop = :hemite    
                iknots[r:s] .= [i for i=ncells+1:ncells+1+s-r]
            end

            # Compute interpolation points using Greville-style averaging
            for i = 1:ntau
                isum = sum( iknots(i+1-degree:i) )
                if isum % degree == 0
                    tau[i] = xmin + isum / degree * dx
                else
                    tau[i] = xmin + isum / degree  * dx
                end
            end

            if degree & 1 > 0
                tau[1]    = xmin
                tau[ntau] = xmax
            end

        end 

    end

end

function test_spline_1d()
    
    tol = 1e-14

    bc_kinds = [:periodic, :hermite, :greville]

    ncells = 10 # number of cells in grid

    for degree = 1:9 # Cycle over spline degree

        for bcs in Iterators.product(bc_kinds, bc_kinds)

            bc_xmin, bc_xmax = bcs

            if !(any( [bc_xmax,bc_xmin] .== :periodic ) & (bc_xmin != bc_xmax)) 
                println( "degree : $degree,  $bc_xmin, $bc_xmax ")
            end

            bspline = Spline1d( ncells, degree, -2π, 2π  )


        end 

    end

end 

test_spline_1d()
