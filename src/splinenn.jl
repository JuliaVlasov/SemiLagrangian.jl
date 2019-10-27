using LinearAlgebra
import LinearAlgebra.LAPACK: pttrf!, pttrs!

mutable struct SplineNN

    geom::Geometry
    a1x::Float64
    a2x::Float64
    a3x::Float64
    a4x::Float64
    a1y::Float64
    a2y::Float64
    a3y::Float64
    a4y::Float64
    axd::Vector{Float64}
    ayd::Vector{Float64}
    axod::Vector{Float64}
    ayod::Vector{Float64}
    aym1gamma1::Vector{Float64}
    aym1gamma2::Vector{Float64}
    axm1gamma1::Vector{Float64}
    axm1gamma2::Vector{Float64}
    coef::Array{Float64,2}
    bcoef::Array{Float64,2}

    function SplineNN(geom::Geometry)

        n1 = geom.n1
        n2 = geom.n2
        n1p1 = n1 + 1
        n1p2 = n1 + 2
        n2p1 = n2 + 1
        n2p2 = n2 + 2

        delta1 = geom.delta1
        delta2 = geom.delta2

        axm1gamma1 = zeros(Float64, n1)
        axm1gamma2 = zeros(Float64, n1)
        aym1gamma1 = zeros(Float64, n2)
        aym1gamma2 = zeros(Float64, n2)
        axd = 4 * ones(Float64, n1)
        axod = ones(Float64, n1 - 1)
        ayd = 4 * ones(Float64, n2)
        ayod = ones(Float64, n2 - 1)

        coef = zeros(Float64, (n1p2, n2p2))
        bcoef = zeros(Float64, (n2p2, n1p2))

        pttrf!(axd, axod)
        pttrf!(ayd, ayod)

        axm1gamma2[1] = 1
        axm1gamma1[n1] = 1

        pttrs!(axd, axod, axm1gamma1)
        pttrs!(axd, axod, axm1gamma2)

        # compute Ay-1.gamma
        aym1gamma2[1] = 1
        aym1gamma1[n2] = 1

        pttrs!(ayd, ayod, aym1gamma1)
        pttrs!(ayd, ayod, aym1gamma2)

        aa1x = 3 / delta1
        aa1y = 3 / delta2
        aa2x = 6 / (delta1 * delta1)
        aa2y = 6 / (delta2 * delta2)

        # assemblage de la matrice 2x2 pour la spline dans la direction Ox

        a1x = aa2x * (1.0 - axm1gamma1[n1-1] + 2 * axm1gamma1[n1])
        a2x = aa2x * (-axm1gamma2[n1-1] + 2 * axm1gamma2[n1])
        a3x = aa2x * (2 * axm1gamma1[1] - axm1gamma1[2])
        a4x = aa2x * (1.0 + 2 * axm1gamma2[1] - axm1gamma2[2])

        # assemblage de la matrice 2x2 pour spline naturels (direction Oy)

        a1y = aa2y * (1.0 - aym1gamma1[n2-1] + 2 * aym1gamma1[n2])
        a2y = aa2y * (-aym1gamma2[n2-1] + 2 * aym1gamma2[n2])
        a3y = aa2y * (2 * aym1gamma1[1] - aym1gamma1[2])
        a4y = aa2y * (1.0 + 2 * aym1gamma2[1] - aym1gamma2[2])

        new(
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
            aym1gamma1,
            aym1gamma2,
            axm1gamma1,
            axm1gamma2,
            coef,
            bcoef,
        )

    end

end

function interpolate!(
    spline::SplineNN,
    f::Array{Float64,2},
    x1::Array{Float64,2},
    x2::Array{Float64,2},
)

    nat_x!(spline, f)
    nat_y!(spline)
    evaltab!(spline, x1, x2, f)

end

"""
interpolation par spline periodique dans les deux directions.
Les points d'interpolation sont definis grace a depx et depy
qui definissent le deplacement par rapport au maillage.
 - f contient les valeurs de la fonction de distribution
 - depx et depy : deplacements par rapport au maillage
     des points dans les quels on veut evaluer la spline.
"""
function interpolate!(
    spline::SplineNN,
    f::Array{Float64,2},
    depx::Float64,
    depy::Float64,
)

    nat_x!(spline, f)
    nat_y!(spline)
    evaldep!(spline, depx, depy, f)

end

"""
    natural splines

"""
function nat_x!(spline::SplineNN, gtau::Array{Float64,2})

    n1 = spline.geom.n1
    n2 = spline.geom.n2
    delta1 = spline.geom.delta1
    delta2 = spline.geom.delta2

    axm1f = similar(gtau)

    n1p2 = n1 + 2
    n1p1 = n1 + 1
    det = spline.a1x * spline.a4x - spline.a2x * spline.a3x

    for j = 1:n2
        for i = 1:n1
            axm1f[i, j] = 6 * gtau[i, j]
        end
    end

    pttrs!(spline.axd, spline.axod, axm1f)

    for j = 1:n2
       # assemblage du second membre du systeme 2x2
        gamma1 = (6.0 / (delta1)^2) * (-axm1f[n1-1, j] + 2 * axm1f[n1, j])
        gamma2 = (6.0 / (delta1)^2) * (2 * axm1f[1, j] - axm1f[2, j])

        coefnp2 = (gamma1 * spline.a4x - gamma2 * spline.a2x) / det
        coef1 = (gamma2 * spline.a1x - gamma1 * spline.a3x) / det
        spline.bcoef[j, n1p2] = coefnp2
        spline.bcoef[j, 1] = coef1

        for i = 2:n1p1
            spline.bcoef[j, i] = axm1f[i-1, j]
            -spline.axm1gamma1[i-1] * coefnp2
            -spline.axm1gamma2[i-1] * coef1
        end
    end
end

function nat_y!(spline::SplineNN)

    n1 = spline.geom.n1
    n2 = spline.geom.n2
    delta1 = spline.geom.delta1
    delta2 = spline.geom.delta2

    n1p2 = n1 + 2
    n2p2 = n2 + 2
    det = spline.a1y * spline.a4y - spline.a2y * spline.a3y

    aym1f = zeros(Float64, (n2, n1p2))
    for i = 1:n1p2
        for j = 1:n2
            aym1f[j, i] = 6 * spline.bcoef[j, i]
        end
    end

    pttrs!(spline.ayd, spline.ayod, aym1f)

    for i = 1:n1p2

        gamma1 = (6.0 / (delta2)^2) * (-aym1f[n2-1, i] + 2 * aym1f[n2, i])
        gamma2 = (6.0 / (delta2)^2) * (2 * aym1f[1, i] - aym1f[2, i])

        coefnp2 = (gamma1 * spline.a4y - gamma2 * spline.a2y) / det
        coef1 = (gamma2 * spline.a1y - gamma1 * spline.a3y) / det

        spline.coef[i, n2p2] = coefnp2
        spline.coef[i, 1] = coef1
        for j = 2:n2+1
            spline.coef[i, j] = aym1f[j-1, i]
            -spline.aym1gamma1[j-1] * coefnp2
            -spline.aym1gamma2[j-1] * coef1
        end

    end
end

function evaltab!(
    spline::SplineNN,
    xd::Array{Float64,2},
    yd::Array{Float64,2},
    f::Array{Float64,2},
)

    n1 = spline.geom.n1
    n2 = spline.geom.n2
    delta1 = spline.geom.delta1
    delta2 = spline.geom.delta2
    delta1x = delta1 * delta1
    delta1xx = delta1x * delta1
    delta1xx6 = 1 / (6 * delta1xx)
    delta2y = delta2 * delta2
    delta2yy = delta2y * delta2
    delta2yy6 = 1 / (6 * delta2yy)
    idelta1 = 1 / delta1
    idelta2 = 1 / delta2

    for j = 2:n2-1
        for i = 2:n1-1

            i1 = trunc(Int64, (xd[i, j] - spline.geom.x1min) * idelta1)
            j1 = trunc(Int64, (yd[i, j] - spline.geom.x2min) * idelta2)

            xdp1 = spline.geom.x1grid[i1+2] - xd[i, j]
            bvalx1 = xdp1 * xdp1 * xdp1
            bvalx2 = delta1xx + 3 * delta1x * xdp1
            +3 * delta1 * xdp1 * xdp1 - 3 * xdp1 * xdp1 * xdp1
            xd1 = xd[i, j] - spline.geom.x1grid[i1+1]
            bvalx3 = delta1xx + 3 * delta1x * xd1
            +3 * delta1 * xd1 * xd1 - 3 * xd1 * xd1 * xd1
            bvalx4 = xd1 * xd1 * xd1

            ydp1 = spline.geom.x2grid[j1+2] - yd[i, j]
            bvaly1 = ydp1 * ydp1 * ydp1
            bvaly2 = delta2yy + 3 * ydp1 * (delta2y + ydp1 * (delta2 - ydp1))
            yd1 = yd[i, j] - geom.x2grid[j1+1]
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

            f[i, j] = delta1xx6 * delta2yy6 * sval

        end
    end

    fill!(view(f, 1, :), 0)
    fill!(view(f, n1, :), 0)

    fill!(view(f, :, 1), 0)
    fill!(view(f, :, n2), 0)

end

function evaldep!(
    spline::SplineNN,
    alphax::Float64,
    alphay::Float64,
    f::Array{Float64,2},
)

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
        intaxsdx = trunc(Int64, -alphax / delta1 + eps(Float64)) - 1
        ideb = trunc(Int64, alphax / delta1) + 2
        ifin = n1 - 1
    else
        intaxsdx = trunc(Int64, -alphax / delta1)
        ideb = 2
        ifin = -trunc(Int64, -alphax / delta1) + n1 - 1
    end

    if (alphay > 0)
        intaysdy = trunc(Int64, -alphay / delta2 + eps(Float64)) - 1
        jdeb = trunc(Int64, alphay / delta2) + 2
        jfin = n2 - 1
    else
        intaysdy = trunc(Int64, -alphay / delta2)
        jdeb = 2
        jfin = -trunc(Int64, -alphay / delta2) + n2 - 1
    end

    xd1 = -alphax - intaxsdx * delta1
    xdp1 = delta1 - xd1
    yd1 = -alphay - intaysdy * delta2
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

    for j = jdeb:jfin
        j1 = j - 1 + intaysdy
        for i = ideb:ifin
            i1 = i - 1 + intaxsdx

            f[i, j] = delta1xx6 * delta2yy6 *
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

    fill!(view(f, 1:ideb-1, :), 0)
    fill!(view(f, ifin+1:n1, :), 0)

    fill!(view(f, :, 1:jdeb-1), 0)
    fill!(view(f, :, jfin+1:n2), 0)

end
