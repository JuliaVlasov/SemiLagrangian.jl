using IterTools
using OffsetArrays
import LinearAlgebra.LAPACK: gbtrf!, gbtrs!

mutable struct BandedMatrix

    n    :: Int64
    kl   :: Int64
    ku   :: Int64
    ipiv :: AbstractVector
    q    :: AbstractMatrix

    function BandedMatrix( n, kl, ku )

        ipiv = zeros(Int64,n)
        q = zeros(Float64, (2*kl+ku+1,n))

        new( n, kl, ku, ipiv, q)

    end

end 

function set_element!( matrix, i, j, a_ij )

    matrix.q[matrix.kl+matrix.ku+1+i-j,j] = a_ij

end

function factorize!( matrix )

     matrix.q, matrix.ipiv = gbtrf!( matrix.kl, matrix.ku, 
         matrix.n, matrix.q )

end

function solve!( matrix, b )

    gbtrs!( 'N', matrix.kl, matrix.ku, matrix.n, matrix.q, matrix.ipiv, b )

end 

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


mutable struct Spline1D

    degree :: Int64
    ncells :: Int64
    nbasis :: Int64
    start  :: Float64
    stop   :: Float64
    step   :: Float64
    bcoef  :: Vector{Float64}

    function Spline1D( ncells :: Int64, degree :: Int64, start, stop  )
      
        nbasis = ncells+degree
        step   = (stop-start) / ncells
        bcoef  = zeros(Float64, ncells+degree)

        new( degree, ncells, nbasis, start, stop, step, bcoef)

    end

end

function get_cell_and_offset( bspl, x )
 
    if x == bspl.start
        return 1, 0.0
    elseif x == bspl.stop 
        return bspl.ncells, 1.0
    else
        offset = (x-bspl.start) / bspl.step  
        icell    = min(trunc( Int64, offset ), bspl.ncells-1)
        return icell + 1, min(offset - icell, 1.0)
    end

end 

@testset " Get cell number and offset " begin

     n       = 10
     degree  = 5
     bspline = Spline1D( n, degree, -1, 1  )
     tau     = -1 .+ 2 .* rand(n)

     err = 0.0
     for x in tau 
         cell, offset = get_cell_and_offset( bspline, x )
         err += abs( x - (bspline.start + (cell-1+offset) * bspline.step))
     end
     @test err ≈ 0.0 atol = 1e-15

end

function eval_basis_and_n_derivs!( derivs, bs :: Spline1D, 
                                   x :: Float64, n :: Int64 )

    ndu = OffsetArray{Float64}( undef, 0:bs.degree, 0:bs.degree)
    a   = OffsetArray{Float64}( undef, 0:1        , 0:bs.degree)

    icell, offset = get_cell_and_offset( bs, x )

    jmin = icell

    ndu[0,0] = 1.0
    for j = 1:bs.degree
       xx     = -offset
       saved  = 0.0
       for r = 0:j-1
          xx       = xx + 1.0
          temp     = ndu[r,j-1] / j
          ndu[r,j] = saved + xx * temp
          saved    = (j - xx) * temp
       end
       ndu[j,j] = saved
    end

    derivs[0,1:bs.degree+1] .= ndu[0:bs.degree,bs.degree]

    for r = 0:bs.degree
       s1 = 0
       s2 = 1
       a[0,0] = 1.0
       for k = 1:n
          d  = 0.0
          rk = r-k
          pk = bs.degree-k
          if (r >= k)
             a[s2,0] = a[s1,0] / (pk+1)
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
             j2 = bs.degree - r
          end
          for j = j1:j2
             a[s2,j] = (a[s1,j] - a[s1,j-1]) / (pk+1)
             d = d + a[s2,j] * ndu[rk+j,pk]
          end
          if (r <= pk)
             a[s2,k] = - a[s1,k-1] / (pk+1)
             d = d + a[s2,k] * ndu[r,pk]
          end
          derivs[k,r+1] = d
          j  = s1
          s1 = s2
          s2 = j
       end
    end

    d = bs.degree / bs.step
    for k = 1:n
       derivs[k,:] .= derivs[k,:] .* d
       d = d * (bs.degree-k) / bs.step
    end

end 


"""
Evaluate value at x of all basis functions with support in local cell
values[j] = B_j(x) for jmin <= j <= jmin+degree

"""
function eval_basis!( self, x, values )

    jmin, offset = get_cell_and_offset( self, x )

    values[1] = 1.0
    for j = 1:self.degree
        xx     = -offset
        saved  = 0.0
        for r = 0:j-1
            xx          = xx + 1.0
            temp        = values[r+1] / j
            values[r+1] = saved + xx * temp
            saved       = (j - xx) * temp
        end
        values[j+1] = saved
    end

    jmin

end 

mutable struct InterpolatorSpline1D

    bspl     :: Spline1D
    tau      :: Vector{Float64}
    matrix   :: BandedMatrix

    function InterpolatorSpline1D( bspl :: Spline1D )

        ntau = bspl.nbasis - bspl.degree÷2 - bspl.degree÷2

        tau    = zeros(Float64, ntau)
        values = zeros(Float64, bspl.degree+1)

        iknots = OffsetArray{Float64}(undef, 2-bspl.degree:ntau)

        r = 2 - bspl.degree
        s = - bspl.degree÷2
        iknots[r:s] .= collect(r-s-1:-1)

        r = - bspl.degree÷2 + 1
        s = - bspl.degree÷2 + 1 + bspl.ncells
        iknots[r:s] .= collect(0:bspl.ncells)

        r = - bspl.degree÷2 + 1 + bspl.ncells + 1
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

        derivs = OffsetArray{Float64}(undef, 
                                  0:bspl.degree÷2, 1:bspl.degree+1)

        x = bspl.start
        jmin = eval_basis_and_n_derivs!( derivs, bspl, x, bspl.degree÷2 )

        h = [bspl.step^i for i=1:ubound(derivs,1)]
        for j = lbound(derivs,2):ubound(derivs,2)
            derivs[1:end,j] .= derivs[1:end,j] .* h[1:end]
        end

        for i = 1:bspl.degree÷2
            order = bspl.degree÷2 - i + (bspl.degree & 1)
            for j = 1:bspl.degree
               set_element!( matrix, i, j, derivs[order,j] )
            end
        end

        for i = bspl.degree÷2+1:bspl.nbasis-bspl.degree÷2
            x = tau[i-bspl.degree÷2]
            jmin  = eval_basis!( bspl, x, values )
            for s = 1:bspl.degree+1
              j = mod( jmin+s-2, bspl.nbasis ) + 1
              set_element!( matrix, i, j, values[s] )
            end
        end

        x = bspl.stop
        jmin = eval_basis_and_n_derivs!( derivs, bspl, x, bspl.degree÷2 )

        h = [bspl.step^i for i=1:ubound(derivs,1)]
        for j = lbound(derivs,2):ubound(derivs,2)
            derivs[1:end,j] .= derivs[1:end,j] .* h[1:end]
        end

        for i = bspl.nbasis-bspl.degree÷2+1:bspl.nbasis
            order = i-(bspl.nbasis-bspl.degree÷2+1) + (bspl.degree & 1)
            j0 = bspl.nbasis-bspl.degree
            d0 = 1
            for s = 1:bspl.degree
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

    nipts + degree÷2 + degree÷2 - degree

end

function compute_interpolant( self, gtau, derivs_xmin, derivs_xmax )

      degree = self.bspl.degree
      odd    = degree & 1
      dx     = self.bspl.step
      bcoef[1:degree÷2] = [derivs_xmin[i]*dx^(i+odd-1) for i=degree÷2:-1:1]

      bcoef[degree÷2+1:nbasis-degree÷2] .= gtau

      bcoef[nbasis-degree÷2+1:nbasis] .= [derivs_xmax(i)*dx^(i+odd-1) for i=1:degree÷2]

      # Solve linear system and compute coefficients
      solve!( matrix, bcoef(1:nbasis) )

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
