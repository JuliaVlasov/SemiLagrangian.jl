using IterTools
using OffsetArrays
using LinearAlgebra.LAPACK

struct BandedMatrix

    n    :: Int64
    kl   :: Int64
    ku   :: Int64
    ipiv :: Vector{Int64}
    q    :: Array{Float64,2}

    function BandedMatrix( n, kl, ku )

        ipiv = zeros(Int64,n)
        q = zeros(Float64, (2*kl+ku+1,n))

        new( n, kl, ku, ipiv, q)

    end

end 

function set_element!( self, i, j, a_ij )

    self.q[self.kl+self.ku+1+i-j,j] = a_ij

end

function get_icell_and_offset( self, x, icell, x_offset )
 
    if (x == self.xmin) 
        icell = 1          
        x_offset = 0.0
    elseif (x == self.xmax) 
        icell = self.ncells
        x_offset = 1.0
    else
        x_offset = (x-self.xmin) * self.inv_dx  
        icell    = trunc( Int64, x_offset )
        x_offset = x_offset - icell
        icell    = icell + 1
    end

    if (icell == self.ncells+1 && x_offset == 0.0)
       icell    = self.ncells
       x_offset = 1.0
    end

end 

function eval_basis_and_n_derivs( self, x, n, derivs, jmin )

    ndu = OffsetArray( undef, 0:self%degree, 0:self%degree)
    a   = OffsetArray( undef, 0:1          , 0:self%degree)

    get_icell_and_offset( self, x, icell, x_offset )

    jmin = icell

    ndu[0,0] = 1.0
    for j = 1:self.degree
       j_real = Float64(j)
       xx     = -x_offset
       saved  = 0.0
       for r = 0:j-1
          xx       = xx + 1.0
          temp     = ndu[r,j-1] * inv[j]
          ndu[r,j] = saved + xx * temp
          saved    = (j_real - xx) * temp
       end
       ndu[j,j] = saved
    end

    derivs[0,:] = ndu[:,spline_degree]

    for r = 0:deg
       s1 = 0
       s2 = 1
       a[0,0] = 1.0
       for k = 1:n
          d  = 0.0
          rk = r-k
          pk = deg-k
          if (r >= k)
             a[s2,0] = a[s1,0] * inv[pk+1]
             d = a[s2,0] * ndu[rk,pk]
          end
          if (rk > -1)
             j1 = 1
          else
             j1 = -rk
          end
          if (r-1 <= pk)
             j2 = k-1
          else
             j2 = deg-r
          end
          for j = j1:j2
             a[s2,j] = (a[s1,j] - a[s1,j-1]) * inv[pk+1]
             d = d + a[s2,j] * ndu[rk+j,pk]
          end
          if (r <= pk)
             a[s2,k] = - a[s1,k-1] * inv[pk+1]
             d = d + a[s2,k] * ndu[r,pk]
          end
          derivs[k,r] = d
          j  = s1
          s1 = s2
          s2 = j
       end
    end

    d = deg * self.inv_dx
    for k = 1:n
       derivs[k,:] = derivs[k,:] * d
       d = d * (deg-k) * self.inv_dx
    end

end 


function factorize!( self )

     m = 2*kl+ku+1
     self.q, self.ipiv = gbtrf!( self.kl, self.ku, m, self.q )

end

function solve!( self, b )

    gbtrs!( "N", self.kl, self.ku, 2*self.kl+self.ku+1, self.q, self.ipiv, b )

end 

mutable struct Spline1D

    degree :: Int64
    ncells :: Int64
    nbasis :: Int64
    start  :: Float64
    stop   :: Float64
    step   :: Float64
    bcoef  :: OffsetArray{Float64,1,Array{Float64,1}}

    function Spline1D( ncells :: Int64, degree :: Int64, start, stop  )
      
        nbasis   = ncells+degree
        step     = (stop-start) / ncells
        bcoef    = OffsetArray{Float64}(undef, 1:ncells+degree)

        new( degree, ncells, nbasis, start, stop, step, bcoef)

    end

end

function eval( spline, x )

    values = zeros(Float64, spline.degree+1)

    offset = (x - spline.start) / spline.step  
    cell   = trunc( Int64, offset )
    offset = offset - cell
    cell   = min( cell+1, spline.stop)

    jmin = icell

    values[1] = 1.0
    for j = 1:spline.degree
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

    jmax = jmin + spline.degree

    self.bcoef[jmin:jmax] .* values

end 

mutable struct InterpolatorSpline1D

    bspl     :: Spline1D
    tau      :: Vector{Float64}
    matrix   :: BandedMatrix

    function InterpolatorSpline1D( bspl :: Spline1D )

        nbc_xmin = bspl.degree ÷ 2
        nbc_xmax = bspl.degree ÷ 2

        @show ntau = bspl.nbasis - nbc_xmin - nbc_xmax

        tau = zeros(Float64, ntau)

        iknots = OffsetArray{Float64}(undef, 2-bspl.degree:ntau)

        r = 2 - bspl.degree
        s = - nbc_xmin
        iknots[r:s] .= collect(r-s-1:-1)

        r = - nbc_xmin + 1
        s = - nbc_xmin + 1 + bspl.ncells
        iknots[r:s] .= collect(0:bspl.ncells)

        r = - nbc_xmin + 1 + bspl.ncells + 1
        s = ntau
        iknots[r:s] .= collect(bspl.ncells+1:bspl.ncells+1+s-r)

        for i = 1:ntau
            isum = sum( iknots[i+1-bspl.degree:i] )
            tau[i] = bspl.start + isum / bspl.degree * bspl.step
        end

        if isodd(bspl.degree)
            tau[1]    = bspl.start
            tau[ntau] = bspl.stop
        end

        ku = max( (bspl.degree+1)÷2, bspl.degree-1 )
        kl = max( (bspl.degree+1)÷2, bspl.degree-1 )

        matrix = BandedMatrix( bspl.nbasis, kl, ku )

        derivs = OffsetArray{Float64}(undef, 0:bspl.degree÷2, 1:bspl.degree+1)

        x = bspl.start
        eval_basis_and_n_derivs( x, nbc_xmin, derivs, jmin )

        h = [bspl.step^i for i=1:derivs[end]]
        for j = lbound(derivs,2):ubound(derivs,2)
            derivs[1:end,j] = derivs[1:end,j] * h[1:end]
        end

        for i = 1:nbc_xmin
            order = nbc_xmin - i + (self.degree & 1)
            for j = 1:bspl.degree
               set_element( matrix, i, j, derivs[order,j] )
            end
        end

        for i = nbc_xmin+1:nbasis-nbc_xmax
            x = self.tau[i-nbc_xmin]
            bspl.eval_basis( x, values, jmin )
            for s = 1:degree+1
              j = mod( jmin-self.offset+s-2, nbasis ) + 1
              set_element( matrix, i, j, values[s] )
            end
        end

        x = bspl.stop
        eval_basis_and_n_derivs( x, nbc_xmax, derivs, jmin )

        h = [self.step^i for i=1:ubound(derivs,1)]
        for j = lbound(derivs,2):ubound(derivs,2)
            derivs[1:end,j] .= derivs[1:end,j] .* h[1:end]
        end

        for i = nbasis-nbc_xmax+1:nbasis
            order = i-(nbasis-nbc_xmax+1)+ (bspl.degree & 1)
            j0 = nbasis-degree
            d0 = 1
            for s = 1:degree
               j = j0 + s
               d = d0 + s
               set_element!( matrix, i, j, derivs[order,d] )
            end
        end

        factorize!(matrix)

        new( bspl, tau, matrix )

    end

end


function compute_num_cells( degree , nipts )


    nbc_xmin = degree / 2
    nbc_xmax = degree / 2

    nipts + nbc_xmin + nbc_xmax - degree

end

function compute_interpolant( self, spline, gtau, derivs_xmin, derivs_xmax )

      bcoef[1:nbc_xmin] = [derivs_xmin[i]*self.dx^(i+self.odd-1) for i=nbc_xmin:-1:1]

      # Interpolation points
      bcoef[nbc_xmin+1+g:nbasis-nbc_xmax+g] .= gtau

      bcoef[nbasis-nbc_xmax+1:nbasis] .= [derivs_xmax(i)*self%dx^(i+self%odd-1) for i=1:nbc_xmax]

      # Solve linear system and compute coefficients
      solve!( matrix, bcoef(1+g:nbasis+g) )

end 

lbound( array :: OffsetArray, dim ) = axes(array)[dim].indices[1]

ubound( array :: OffsetArray, dim ) = axes(array)[dim].indices[end]


function test_spline_1d()
    
    ncells = 10 # number of cells in grid

    for degree = 1:9 # Cycle over spline degree

        bspline = Spline1D( ncells, degree, -2π, 2π  )
        interpolator = InterpolatorSpline1D( bspline )

    end

end 

test_spline_1d()
