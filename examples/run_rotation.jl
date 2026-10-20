using LinearAlgebra
using SemiLagrangian

function exact_rotation!(f, mesh_x, mesh_y, t)
    s, c = sincos(t)
    for (i, x) in enumerate(mesh_x.points), (j, y) in enumerate(mesh_y.points)
        x0, y0 = c * x - s * y, s * x + c * y
        f[i, j] = exp(-13 * (x0^2 + (y0 + 1.2)^2))
    end
    return f
end

function run_rotation(; n = 256, nsteps = 50)
    mesh_x = UniformMesh(-5.0, 5.0, n)
    mesh_y = UniformMesh(-5.0, 5.0, n)
    interp = Lagrange(11)
    dt = 2pi / nsteps
    x, y = mesh_x.points, mesh_y.points
    data = [exp(-13 * (x0^2 + (y0 + 1.2)^2)) for x0 in x, y0 in y]
    exact = similar(data)

    for step in 1:nsteps
        for (i, coefficient) in enumerate(magicsplit(dt))
            if isodd(i)
                advection!(data, mesh_x, y, coefficient, interp; dim = 1)
            else
                advection!(data, mesh_y, -x, coefficient, interp; dim = 2)
            end
        end
        exact_rotation!(exact, mesh_x, mesh_y, step * dt)
        @info "rotation step $step" error = norm(data - exact, Inf)
    end
    return data
end

run_rotation()
