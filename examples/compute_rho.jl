import Statistics: mean

function compute_rho(
        meshv::UniformMesh,
        f::Array{Float64, 2}
    )

    dv = step(meshv)
    rho = dv * sum(f, dims = 2)
    return vec(rho .- mean(rho))
end
