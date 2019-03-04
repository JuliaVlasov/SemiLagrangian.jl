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

    function Spline1D( n :: Int64, p :: Int64, xmin, xmax  )
      
        degree   = p
        ncells   = n
        nbasis   = ncells+degree
        start    = xmin
        stop     = xmax
        step     = (xmax-xmin) / ncells
        bcoef    = OffsetArray{Float64}(undef, 1:n+p)

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
    matrix   :: Array{Float64, 2}

    function InterpolatorSpline1D( bspl :: Spline1D )

        nbc_xmin = bspl.degree ÷ 2
        nbc_xmax = bspl.degree ÷ 2
        odd      = bspl.degree & 1

        ntau = bspl.nbasis - bspl.degree

        tau = zeros(Float64, ntau)

        iknots = OffsetArray{Float64}(undef, 2-bspl.degree:ntau)

        r = 2 - bspl.degree
        s = - nbc_xmin
        iknots[r:s] .= [i for i=r-s-1:-1]

        r = - nbc_xmin + 1
        s = - nbc_xmin + 1 + bspl.ncells
        iknots[r:s] .= [i for i=0:bspl.ncells]

        r = - nbc_xmin + 1 + bspl.ncells + 1
        s = ntau
        iknots[r:s] .= [i for i=bspl.ncells+1:bspl.ncells+1+s-r]

        # Compute interpolation points using Greville-style averaging
        for i = 1:ntau
            isum = sum( iknots[i+1-bspl.degree:i] )
            if isum % degree == 0
                tau[i] = bspl.xmin + isum / bspl.degree * dx
            else
                tau[i] = bspl.xmin + isum / bspl.degree  * dx
            end
        end

        if Bool(odd)
            tau[1]    = bspl.xmin
            tau[ntau] = bspl.xmax
        end

        ku = max( (degree+1)÷2, degree-1 )
        kl = max( (degree+1)÷2, degree-1 )

        spline_matrix_new( self.matrix, nbasis, kl, ku )

        build_system( self, self.matrix )

        factorize!(self.matrix)

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

function build_system( self, matrix )

    derivs = OffsetArray{Float64}(undef, 0:degree/2, 1:degree+1)

    x = self.bspl.xmin
    eval_basis_and_n_derivs( x, nbc_xmin, derivs, jmin )

    h = [self.dx^i for i=1:derivs[end]]
    for j = lbound(derivs,2):ubound(derivs,2)
        derivs[1:end,j] = derivs[1:end,j] * h[1:end]
    end

    for i = 1:nbc_xmin
        order = nbc_xmin-i+self.odd
        for j = 1:self.degree
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

    x = self.bspl.xmax
    bspl.eval_basis_and_n_derivs( x, nbc_xmax, derivs, jmin )

    h = [self.dx^i for i=1:ubound(derivs,1)]
    for j = lbound(derivs,2):ubound(derivs,2)
        derivs[1:end,j] .= derivs[1:end,j] .* h[1:end]
    end

    for i = nbasis-nbc_xmax+1:nbasis
        order = i-(nbasis-nbc_xmax+1)+self%odd
        j0 = nbasis-degree
        d0 = 1
        for s = 1:degree
           j = j0 + s
           d = d0 + s
           matrix.set_element( i, j, derivs[order,d] )
        end
    end

end 

function test_spline_1d()
    
    ncells = 10 # number of cells in grid

    for degree = 1:9 # Cycle over spline degree

        bspline = Spline1D( ncells, degree, -2π, 2π  )
        interpolator = InterpolatorSpline1D( bspline )

    end

end 

test_spline_1d()
