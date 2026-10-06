import Statistics: mean

function compute_rho(
        meshv::UniformMesh,
        f::Array{Float64, 2}
    )

    dv = meshv.dx
    rho = dv * sum(f, dims = 2)
    return vec(rho .- mean(rho))
end

