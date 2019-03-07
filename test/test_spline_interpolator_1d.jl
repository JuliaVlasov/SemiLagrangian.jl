using OffsetArrays

lbound( array :: OffsetArray, dim ) = axes(array)[dim].indices[1]

ubound( array :: OffsetArray, dim ) = axes(array)[dim].indices[end]

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

function compute_interpolant( self, gtau, derivs_xmin, derivs_xmax )

      degree = self.bspl.degree
      odd    = degree & 1
      dx     = self.bspl.step
      bcoef[1:degree÷2] = [derivs_xmin[i]*dx^(i+odd-1) for i=degree÷2:-1:1]

      bcoef[degree÷2+1:nbasis-degree÷2] .= gtau

      bcoef[nbasis-degree÷2+1:nbasis] .= [derivs_xmax(i)*dx^(i+odd-1) for i=1:degree÷2]

      solve!( matrix, bcoef(1:nbasis) )

end 


