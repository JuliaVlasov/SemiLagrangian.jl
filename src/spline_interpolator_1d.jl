using OffsetArrays

lbound( array :: OffsetArray, dim ) = first(axes(array)[dim])
ubound( array :: OffsetArray, dim ) = last(axes(array)[dim])

export InterpolatorSpline1D

struct InterpolatorSpline1D

    bspl::Spline1D
    tau::Vector{Float64}
    matrix::BandedMatrix

    function InterpolatorSpline1D(degree, xmin, xmax, ncells)

        bspl = Spline1D(ncells, degree, xmin, xmax)

        ntau = bspl.nbasis - degree ÷ 2 - degree ÷ 2

        tau = zeros(Float64, ntau)
        values = zeros(Float64, degree + 1)

        iknots = OffsetArray{Float64}(undef, 2-degree:ntau)

        r = 2 - degree
        s = -degree ÷ 2
        iknots[r:s] .= collect(r-s-1:-1)

        r = -degree ÷ 2 + 1
        s = -degree ÷ 2 + 1 + bspl.ncells
        iknots[r:s] .= collect(0:bspl.ncells)

        r = -degree ÷ 2 + 1 + bspl.ncells + 1
        s = ntau
        iknots[r:s] .= collect(bspl.ncells+1:bspl.ncells+1+s-r)

        for i = 1:ntau
            isum = sum(iknots[i+1-degree:i])
            tau[i] = bspl.start + isum / degree * bspl.step
        end

        if isodd(degree)
            tau[1] = bspl.start
            tau[ntau] = bspl.stop
        end

        ku = max((degree + 1) ÷ 2, degree - 1)
        kl = max((degree + 1) ÷ 2, degree - 1)

        matrix = BandedMatrix(bspl.nbasis, kl, ku)

        derivs = OffsetArray{Float64}(undef, 0:degree÷2, 1:degree+1)

        x = bspl.start
        jmin = eval_basis_and_n_derivs!(derivs, bspl, x, degree ÷ 2)

        h = [bspl.step^i for i = 1:ubound(derivs, 1)]
        for j = lbound(derivs, 2):ubound(derivs, 2)
            derivs[1:end, j] .= derivs[1:end, j] .* h[1:end]
        end

        for i = 1:degree÷2
            order = degree ÷ 2 - i + (degree & 1)
            for j = 1:degree
                set_element!(matrix, i, j, derivs[order, j])
            end
        end

        for i = degree÷2+1:bspl.nbasis-degree÷2
            x = tau[i-degree÷2]
            jmin = eval_basis!(bspl, x, values)
            for s = 1:degree+1
                j = mod(jmin + s - 2, bspl.nbasis) + 1
                set_element!(matrix, i, j, values[s])
            end
        end

        x = bspl.stop
        jmin = eval_basis_and_n_derivs!(derivs, bspl, x, degree ÷ 2)

        h = [bspl.step^i for i = 1:ubound(derivs, 1)]
        for j = lbound(derivs, 2):ubound(derivs, 2)
            derivs[1:end, j] .= derivs[1:end, j] .* h[1:end]
        end

        for i = bspl.nbasis-degree÷2+1:bspl.nbasis
            order = i - (bspl.nbasis - degree ÷ 2 + 1) + (degree & 1)
            j0 = bspl.nbasis - degree
            d0 = 1
            for s = 1:degree
                j = j0 + s
                d = d0 + s
                set_element!(matrix, i, j, derivs[order, d])
            end
        end

        factorize!(matrix)

        new(bspl, tau, matrix)

    end

end

export compute_interpolant!

function compute_interpolant!(self, gtau, derivs_left, derivs_right)

    degree = self.bspl.degree
    odd = degree & 1
    dx = self.bspl.step

    self.bspl.bcoef[1:degree÷2] .= [derivs_left[i] * dx^(i + odd - 1) for i = degree÷2:-1:1]

    self.bspl.bcoef[end-degree÷2+1:end] .= [derivs_right[i] * dx^(i + odd - 1) for i = 1:degree÷2]

    for i in eachindex(gtau)
        self.bspl.bcoef[i+degree÷2] = gtau[i]
    end

    solve!(self.matrix, self.bspl.bcoef)

end
