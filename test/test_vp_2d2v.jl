struct VlasovPoisson

    xmin::Vector{Float64}
    xmax::Vector{Float64}
    ncells::Vector{Int64}

    interp_x1::InterpolatorSpline1D
    interp_x2::InterpolatorSpline1D
    interp_x3::InterpolatorSpline1D
    interp_x4::InterpolatorSpline1D

    ex::Array{Float64,2}
    ey::Array{Float64,2}
    rho::Array{Float64,2}

    fxv::Array{Float64,4}
    fvx::Array{Float64,4}

    function VlasovPoisson(degree, xmin, xmax, ncells)

        inter_x1 = InterpolatorSpline1D(degree, xmin[1], xmax[1], ncells[1])
        inter_x2 = InterpolatorSpline1D(degree, xmin[2], xmax[2], ncells[2])
        inter_x3 = InterpolatorSpline1D(degree, xmin[3], xmax[3], ncells[3])
        inter_x4 = InterpolatorSpline1D(degree, xmin[4], xmax[4], ncells[4])

        ex = zeros(Float64, (ncells[1], ncells[2]))
        ey = zeros(Float64, (ncells[1], ncells[2]))
        rho = zeros(Float64, (ncells[1], ncells[2]))

        fxv = zeros(ComplexF64, (ncells...))
        fvx = permutedims(fxv, [3, 4, 1, 2])

    end

end

function advection_x1!(vlasov::VlasovPoisson, dt::Float64)

    for l = 1:vlasov.ncells[4], k = 1:vlasov.ncells[3]
        alpha = (vlasov.xmin[3] + (k - 1) * vlasov.delta[3]) * dt
        for j = 1:vlasov.ncells[2]
            interpolate!(vlasov.fxv[:, j, k, l], vlasov.interp_x1, alpha)
        end
    end

end

function advection_x2!(vlasov::VlasovPoisson, dt::Float64)

    for l = 1:vlasov.ncells[4]
        alpha = (xmin[4] + (l - 1) * vlasov.delta[4]) * dt
        for k = 1:vlasov.ncells[3], i = 1:vlasov.ncells[1]
            interpolate!(vlasov.fxv[i, :, k, l], vlasov.interp_x2, alpha)
        end
    end

end

function advection_x3!(vlasov::VlasovPoisson, dt::Float64)

    for l = 1:vlasov.ncells[4], j = 1:vlasov.ncells[2], i = 1:vlasov.ncells[1]
        alpha = vlasov.ex[i, j] * dt
        interpolate!(vlasov.fvx[:, l, i, j], vlasov.interp_x3, alpha)
    end

end

function advection_x4!(vlasov::VlasovPoisson, dt::Float64)

    for k = 1:vlasov.ncells[3], j = 1:vlasov.ncells[2], i = 1:vlasov.ncells[1]
        alpha = vlasov.ey[i, j] * dt
        interpolate!(vlasov.fvx[k, :, i, j], vlasov.interp_x4, alpha)
    end

end

function transposexv(vlasov)

    permutedims!(vlasov.fvx, vlasov.fxv, [3, 4, 1, 2])

end

function transposevx(vlasov)

    permutedims!(vlasov.fxv, vlasov.fvx, [3, 4, 1, 2])

end

function landau!(vlasov; eps = 0.05)

    eps = 0.05
    kx = 2π / (vlasov.xmax[1] - vlasov.xmin[1])
    ky = 2π / (vlasov.xmax[2] - vlasov.xmin[2])

    dx = (xmax .- xmin) ./ ncells

    x4 = x4_min
    for i4 = 1:vlasov.ncells[4]
        x3 = xmin[3]
        for i3 = 1:vlasov.ncells[3]
            x2 = xmin[2]
            v2 = x3 * x3 + x4 * x4
            for i2 = 1:vlasov.ncells[2]
                x1 = xmin[1]
                for i1 = 1:vlasov.ncells[1]
                    vlasov.fxv[i1, i2, i3, i4] = (1 + eps * cos(kx * x1) *
                                                  cos(ky * x2)) / (2π) *
                                                 exp(-.5 * v2)
                    x1 += dx[1]
                end
                x2 += dx[2]
            end
            x3 += dx[3]
        end
        x4 += dx[4]
    end


end


function run_vp2d2v(dt::Float64)

    degree = 5

    xmin = [0.0, 0.0, -6.0, -6.0]
    xmax = [4π, 4π, 6.0, 6.0]
    ncells = [64, 64, 64, 64]

    vlasov = VlasovPoisson(degree, xmin, xmax, ncells)
    landau!(vlasov, eps = 0.05)

    advection_x1!(vlasov, 0.5dt)
    advection_x2!(vlasov, 0.5dt)

    transposexv!(vlasov)

    compute_charge!(vlasov)

    poisson!(vlasov)

    advection_x3!(vlasov, dt)
    advection_x4!(vlasov, dt)

    transposevx!(vlasov)

    advection_x1!(vlasov, dt)
    advection_x2!(vlasov, dt)
    true


end


@test run_vp2d2v(0.01)
