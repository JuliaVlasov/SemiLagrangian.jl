using FFTW

function compute_e(meshx::UniformMesh, rho::Vector{Float64})
    nx = meshx.nx
    k = 2pi / (meshx.xmax - meshx.xmin)
    modes = zeros(Float64, nx)
    modes .= k * fftfreq(nx, nx)
    modes[1] = 1.0
    rhok = fft(rho) ./ modes
    rhok .*= -1im
    ifft!(rhok)
    return real(rhok)
end

