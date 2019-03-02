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
         end
         values[j] = saved
      end

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

        nbc_start = degree÷2
        nbc_stop  = degree÷2

        ntau = nbasis - nbc_start - nbc_stop

        tau = zeros(ntau)

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

function test_spline_1d()
    
    tol = 1e-14

    ncells = 10 # number of cells in grid

    for degree = 1:9 # Cycle over spline degree

        bspline = Spline1d( ncells, degree, -2π, 2π  )

    end

end 

function compute_num_cells( degree , nipts )


    nbc_xmin = degree / 2
    nbc_xmax = degree / 2

    nipts + nbc_xmin + nbc_xmax - degree

end

function init( self, bspl )

    nbc_xmin = degree/2
    nbc_xmax = degree/2
    odd      = modulo( bspl%degree  , 2 )

    call compute_interpolation_points_uniform( self, self.tau )
    call compute_num_diags_uniform( self, kl, ku )

    spline_matrix_new( self.matrix, nbasis, kl, ku )

    build_system( self, self.matrix )

    factorize!(self.matrix)

end 

function compute_interpolant( self, spline, gtau, derivs_xmin, derivs_xmax )

      bcoef[1:nbc_xmin] = [derivs_xmin[i]*self.dx^(i+self.odd-1) for i=nbc_xmin:-1:1]

      # Interpolation points
      bcoef[nbc_xmin+1+g:nbasis-nbc_xmax+g] .= gtau

      bcoef[nbasis-nbc_xmax+1:nbasis] .= [derivs_xmax(i)*self%dx^(i+self%odd-1) for i=1:nbc_xmax]

      # Solve linear system and compute coefficients
      solve!( matrix, bcoef(1+g:nbasis+g) )

end 

function build_system( self, matrix )

    derivs = zeros(Float64, (0:degree/2, 1:degree+1))

    x = self.bspl.xmin
    eval_basis_and_n_derivs( x, nbc_xmin, derivs, jmin )

    h = [self.dx^i for i=1:derivs[end]]
    for j = lbound(derivs,2):ubound(derivs,2)
        derivs[1:,j] = derivs[1:,j] * h[1:]
    end

    for i = 1:nbc_xmin
        order = nbc_xmin-i+self.odd
        for j = 1, degree
           set_element( matrix, i, j, derivs(order,j) )
        end
    end

    for i = nbc_xmin+1, nbasis-nbc_xmax
       x = self%tau(i-nbc_xmin)
       bspl.eval_basis( x, values, jmin )
       for s = 1:degree+1
         j = mod( jmin-self.offset+s-2, nbasis ) + 1
         set_element( matrix, i, j, values[s] )
       end
    end

    x = self.bspl.xmax
    bspl.eval_basis_and_n_derivs( x, nbc_xmax, derivs, jmin )

    h = [(self%dx**i, i=1, ubound(derivs,1))]
    for j = lbound(derivs,2), ubound(derivs,2)
       derivs(1:,j) = derivs(1:,j) * h(1:)
    end

    for i = nbasis-nbc_xmax+1, nbasis
       order = i-(nbasis-nbc_xmax+1)+self%odd
       j0 = nbasis-degree
       d0 = 1
       for s = 1:degree
          j = j0 + s
          d = d0 + s
          matrix.set_element( i, j, derivs(order,d) )
       end
    end

      end

  end 

  function compute_interpolation_points_uniform( self, tau )


      # Determine size of tau and allocate tau
      ntau = nbasis - nbc_xmin - nbc_xmax
      tau = zeros(Float64, ntau )

      # Non-periodic case: create array of temporary knots (integer shifts only)
      # in order to compute interpolation points using Greville-style averaging:
      # tau(i) = xmin + average(knots(i+1-degree:i)) * dx
      iknots = zeros(2-degree:ntau)

      # Additional knots near x=xmin
      r = 2-degree
      s = -nbc_xmin
      iknots[r:s] = [i for i=r-s-1:-1]

      # Knots inside the domain
      r = -nbc_xmin+1
      s = -nbc_xmin+1+ncells
      iknots[r:s] = [i for i=0:ncells]

      # Additional knots near x=xmax
      r = -nbc_xmin+1+ncells+1, s = ntau 
      iknots(r:s) = [i for i=ncells+1:ncells+1+s-r]

      # Compute interpolation points using Greville-style averaging
      inv_deg = 1.0 / degree
      for i = 1:ntau
         isum = sum( iknots[i+1-degree:i] )
         if (mod( isum, degree ) == 0)
             tau[i] = xmin + isum ÷ degree * dx
         else
             tau[i] = xmin + isum * inv_deg * dx
         end
      end

      # Non-periodic case, odd degree: fix round-off issues
      if ( self%odd == 1 )
        tau(1)    = xmin
        tau(ntau) = xmax
      end

  end 

  function compute_num_diags_uniform( self, kl, ku )

      ku = max( (degree+1)/2, degree-1 )
      kl = max( (degree+1)/2, degree-1 )

  end 


end 

test_spline_1d()
