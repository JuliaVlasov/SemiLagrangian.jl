# include("lapack.jl")
using LinearAlgebra

mutable struct SplinePP{T}

    geom::Geometry{T}
    a1x::T
    a2x::T
    a3x::T
    a4x::T
    a1y::T
    a2y::T
    a3y::T
    a4y::T
    axd::Vector{T}
    ayd::Vector{T}
    axod::Vector{T}
    ayod::Vector{T}
    axm1gamma::Array{T,2}
    aym1gamma::Array{T,2}
    coef::Array{T,2}
    bcoef::Array{T,2}

    function SplinePP(geom::Geometry{T}) where{T<:AbstractFloat}

        n1 = geom.n1
        n2 = geom.n2
        n1p1 = n1 + 1
        n1p2 = n1 + 2
        n1p3 = n1 + 3
        n2p1 = n2 + 1
        n2p2 = n2 + 2
        n2p3 = n2 + 3

        axm1gamma = zeros(T, (n1p1, 2))
        aym1gamma = zeros(T, (n2p1, 2))
        axd = 4 * ones(T, n1p1)
        axod = ones(T, n1)
        ayd = 4 * ones(T, n2p1)
        ayod = ones(T, n2)
        coef = zeros(T, (n1p3, n2p3))
        bcoef = zeros(T, (n1p3, n2))

        pttrfgen!(axd, axod)
        pttrfgen!(ayd, ayod)

        axm1gamma[1, 2] = 1
        axm1gamma[n1p1, 1] = 1

        pttrsgen!(axd, axod, axm1gamma)

        aym1gamma[1, 2] = 1
        aym1gamma[n2p1, 1] = 1

        pttrsgen!(ayd, ayod, aym1gamma)

        aa1x = 3 / geom.delta1
        aa1y = 3 / geom.delta2
        aa2x = 6 / (geom.delta1 * geom.delta1)
        aa2y = 6 / (geom.delta2 * geom.delta2)

        a1x = -aa1x * (1.0 + axm1gamma[2, 1] + axm1gamma[n1, 1])
        a2x = -aa1x * (1.0 + axm1gamma[2, 2] + axm1gamma[n1, 2])
        a3x = aa2x *
              (-1.0 + 2 * axm1gamma[1, 1] - axm1gamma[2, 1] + axm1gamma[n1, 1] -
               2 * axm1gamma[n1p1, 1])
        a4x = aa2x *
              (1.0 + 2 * axm1gamma[1, 2] - axm1gamma[2, 2] + axm1gamma[n1, 2] -
               2 * axm1gamma[n1p1, 2])

        a1y = -aa1y * (1.0 + aym1gamma[2, 1] + aym1gamma[n2, 1])
        a2y = -aa1y * (1.0 + aym1gamma[2, 2] + aym1gamma[n2, 2])
        a3y = aa2y *
              (-1.0 + 2 * aym1gamma[1, 1] - aym1gamma[2, 1] + aym1gamma[n2, 1] -
               2 * aym1gamma[n2p1, 1])
        a4y = aa2y *
              (1.0 + 2 * aym1gamma[1, 2] - aym1gamma[2, 2] + aym1gamma[n2, 2] -
               2 * aym1gamma[n2p1, 2])
        new{T}(
            geom,
            a1x,
            a2x,
            a3x,
            a4x,
            a1y,
            a2y,
            a3y,
            a4y,
            axd,
            ayd,
            axod,
            ayod,
            axm1gamma,
            aym1gamma,
            coef,
            bcoef,
        )
    end
end

function interpolate!(
    spline::SplinePP,
    f::Array{T,2},
    x::Array{T,2},
    y::Array{T,2},
) where {T<:AbstractFloat}

    per_x!(spline, f)

    per_y!(spline)

    evaltab!(spline, x, y, f)

end

function interpolate!(
    spline::SplinePP,
    f::Array{T,2},
    depx::T,
    depy::T,
) where {T<:AbstractFloat}

    per_x!(spline, f)
    per_y!(spline)
    evaldep!(spline, depx, depy, f)

end

function per_x!(spline::SplinePP, gtau::Array{T,2}) where {T<:AbstractFloat}

    n1 = spline.geom.n1
    n2 = spline.geom.n2
    n1p2 = n1 + 2
    n1p3 = n1 + 3
    det = spline.a1x * spline.a4x - spline.a2x * spline.a3x

    axm1f = zeros(T, (spline.geom.n1 + 1, spline.geom.n2))

    for j = 1:n2
        for i = 1:n1
            axm1f[i, j] = 6 * gtau[i, j]
        end
        axm1f[n1+1, j] = 6 * gtau[1, j]
    end

    pttrsgen!(spline.axd, spline.axod, axm1f)

    for j = 1:n2

        gamma1 = -(3.0 / spline.geom.delta1) *
                 (axm1f[2, j] + axm1f[spline.geom.n1, j])
        gamma2 = (6.0 / (spline.geom.delta1)^2) *
                 (2 * axm1f[1, j] - axm1f[2, j] + axm1f[spline.geom.n1, j] -
                  2 * axm1f[spline.geom.n1+1, j])

        spline.bcoef[n1p3, j] = (gamma1 * spline.a4x - gamma2 * spline.a2x) /
                                det
        spline.bcoef[1, j] = (gamma2 * spline.a1x - gamma1 * spline.a3x) / det

        for i = 2:n1p2
            spline.bcoef[i, j] = axm1f[i-1, j]
            -spline.axm1gamma[i-1, 1] * spline.bcoef[n1p3, j]
            -spline.axm1gamma[i-1, 2] * spline.bcoef[1, j]
        end

    end

end

function per_y!(spline::SplinePP{T}) where{T<:AbstractFloat}

    n1 = spline.geom.n1
    n2 = spline.geom.n2
    n1p3 = spline.geom.n1 + 3
    n2p3 = spline.geom.n2 + 3
    det = spline.a1y * spline.a4y - spline.a2y * spline.a3y

    delta1 = spline.geom.delta1
    delta2 = spline.geom.delta2

    aym1f = zeros(T, (n2 + 1, n1 + 3))

    for i = 1:n1p3
        for j = 1:n2
            aym1f[j, i] = 6 * spline.bcoef[i, j]
        end
        aym1f[n2+1, i] = 6 * spline.bcoef[i, 1]
    end

    pttrsgen!(spline.ayd, spline.ayod, aym1f)

    for i = 1:n1p3

        gamma1 = -(3.0 / delta2) * (aym1f[2, i] + aym1f[n2, i])
        gamma2 = (6.0 / delta2^2) *
                 (2 * aym1f[1, i] - aym1f[2, i] + aym1f[n2, i] -
                  2 * aym1f[n2+1, i])

        spline.coef[i, n2p3] = (gamma1 * spline.a4y - gamma2 * spline.a2y) / det
        spline.coef[i, 1] = (gamma2 * spline.a1y - gamma1 * spline.a3y) / det

        for j = 2:spline.geom.n2+2
            spline.coef[i, j] = aym1f[j-1, i]
            -spline.aym1gamma[j-1, 1] * spline.coef[i, n2p3]
            -spline.aym1gamma[j-1, 2] * spline.coef[i, 1]
        end
    end

end

function evaltab!(
    spline::SplinePP{T},
    xd::Array{T,2},
    yd::Array{T,2},
    fout::Array{T,2},
) where {T<:AbstractFloat}

    n1 = spline.geom.n1
    n2 = spline.geom.n2

    delta1 = spline.geom.delta1
    delta1x = delta1 * delta1
    delta1xx = delta1x * delta1
    delta1xx6 = 1 / (6 * delta1xx)

    delta2 = spline.geom.delta2
    delta2y = delta2 * delta2
    delta2yy = delta2y * delta2
    delta2yy6 = 1 / (6 * delta2yy)

    idelta1 = 1 / delta1
    idelta2 = 1 / delta2

    lx = (n1 - 1) * delta1
    ly = (n2 - 1) * delta2

    for j = 1:n2

        for i = 1:n1

            i1 = floor(Int64, (xd[i, j] - spline.geom.x1min) * idelta1)
            j1 = floor(Int64, (yd[i, j] - spline.geom.x2min) * idelta2)

            xdp1 = spline.geom.x1grid[i1+2] - xd[i, j]
            bvalx1 = xdp1 * xdp1 * xdp1
            bvalx2 = (delta1xx + 3 * delta1x * xdp1 + 3 * delta1 * xdp1 * xdp1 -
                      3 * xdp1 * xdp1 * xdp1)
            xd1 = xd[i, j] - spline.geom.x1grid[i1+1]
            bvalx3 = (delta1xx + 3 * delta1x * xd1 + 3 * delta1 * xd1 * xd1 -
                      3 * xd1 * xd1 * xd1)
            bvalx4 = xd1 * xd1 * xd1
            ydp1 = spline.geom.x2grid[j1+2] - yd[i, j]
            bvaly1 = ydp1 * ydp1 * ydp1
            bvaly2 = delta2yy + 3 * ydp1 * (delta2y + ydp1 * (delta2 - ydp1))
            yd1 = yd[i, j] - spline.geom.x2grid[j1+1]
            bvaly3 = delta2yy + 3 * yd1 * (delta2y + yd1 * (delta2 - yd1))
            bvaly4 = yd1 * yd1 * yd1

            sval = 0.0
            sval1 = spline.coef[i1+1, j1+1] * bvaly1
            sval1 = sval1 + spline.coef[i1+1, j1+2] * bvaly2
            sval1 = sval1 + spline.coef[i1+1, j1+3] * bvaly3
            sval1 = sval1 + spline.coef[i1+1, j1+4] * bvaly4
            sval = sval + sval1 * bvalx1

            sval2 = spline.coef[i1+2, j1+1] * bvaly1
            sval2 = sval2 + spline.coef[i1+2, j1+2] * bvaly2
            sval2 = sval2 + spline.coef[i1+2, j1+3] * bvaly3
            sval2 = sval2 + spline.coef[i1+2, j1+4] * bvaly4
            sval = sval + sval2 * bvalx2

            sval3 = spline.coef[i1+3, j1+1] * bvaly1
            sval3 = sval3 + spline.coef[i1+3, j1+2] * bvaly2
            sval3 = sval3 + spline.coef[i1+3, j1+3] * bvaly3
            sval3 = sval3 + spline.coef[i1+3, j1+4] * bvaly4
            sval = sval + sval3 * bvalx3

            sval4 = spline.coef[i1+4, j1+1] * bvaly1
            sval4 = sval4 + spline.coef[i1+4, j1+2] * bvaly2
            sval4 = sval4 + spline.coef[i1+4, j1+3] * bvaly3
            sval4 = sval4 + spline.coef[i1+4, j1+4] * bvaly4
            sval = sval + sval4 * bvalx4

            fout[i, j] = delta1xx6 * delta2yy6 * sval

        end
    end

end

function evaldep!(
    spline::SplinePP{T},
    alphax::T,
    alphay::T,
    fout::Array{T,2},
) where {T<:AbstractFloat}

    n1 = spline.geom.n1
    n2 = spline.geom.n2

    delta1 = spline.geom.delta1
    delta1x = delta1 * delta1
    delta1xx = delta1x * delta1
    delta1xx6 = 1 / (6 * delta1xx)

    delta2 = spline.geom.delta2
    delta2y = delta2 * delta2
    delta2yy = delta2y * delta2
    delta2yy6 = 1 / (6 * delta2yy)

    if (alphax > 0)
        intaxsdelta1 = trunc(Int64, -alphax / delta1 + eps(T)) - 1
    else
        intaxsdelta1 = trunc(Int64, -alphax / delta1)
    end

    if (alphay > 0)
        intaysdelta2 = trunc(Int64, -alphay / delta2 + eps(T)) - 1
    else
        intaysdelta2 = trunc(Int64, -alphay / delta2)
    end

    xd1 = -alphax - intaxsdelta1 * delta1
    xdp1 = delta1 - xd1
    yd1 = -alphay - intaysdelta2 * delta2
    ydp1 = delta2 - yd1

    bvalx1 = xdp1 * xdp1 * xdp1
    bvalx2 = delta1xx + 3 * delta1x * xdp1 + 3 * delta1 * xdp1 * xdp1 -
             3 * xdp1 * xdp1 * xdp1
    bvalx3 = delta1xx + 3 * delta1x * xd1 + 3 * delta1 * xd1 * xd1 -
             3 * xd1 * xd1 * xd1
    bvalx4 = xd1 * xd1 * xd1
    bvaly1 = ydp1 * ydp1 * ydp1
    bvaly2 = delta2yy + 3 * ydp1 * (delta2y + ydp1 * (delta2 - ydp1))
    bvaly3 = delta2yy + 3 * yd1 * (delta2y + yd1 * (delta2 - yd1))
    bvaly4 = yd1 * yd1 * yd1

    for j = 1:n2

        j1 = mod(n2 + j - 1 + intaysdelta2, n2)

        for i = 1:n1

            i1 = mod(n1 + i - 1 + intaxsdelta1, n1)

            fout[i, j] = delta1xx6 * delta2yy6 *
                         (bvalx1 * (spline.coef[i1+1, j1+1] * bvaly1 +
                           spline.coef[i1+1, j1+2] * bvaly2 +
                           spline.coef[i1+1, j1+3] * bvaly3 +
                           spline.coef[i1+1, j1+4] * bvaly4) +
                          bvalx2 * (spline.coef[i1+2, j1+1] * bvaly1 +
                           spline.coef[i1+2, j1+2] * bvaly2 +
                           spline.coef[i1+2, j1+3] * bvaly3 +
                           spline.coef[i1+2, j1+4] * bvaly4) +
                          bvalx3 * (spline.coef[i1+3, j1+1] * bvaly1 +
                           spline.coef[i1+3, j1+2] * bvaly2 +
                           spline.coef[i1+3, j1+3] * bvaly3 +
                           spline.coef[i1+3, j1+4] * bvaly4) +
                          bvalx4 * (spline.coef[i1+4, j1+1] * bvaly1 +
                           spline.coef[i1+4, j1+2] * bvaly2 +
                           spline.coef[i1+4, j1+3] * bvaly3 +
                           spline.coef[i1+4, j1+4] * bvaly4))
        end
    end

end
