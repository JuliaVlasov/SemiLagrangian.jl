using FFTW

function compute_e(meshx::UniformMesh, rho::Vector{Float64})
    nx = length(meshx)
    k = 2pi / meshx.width
    modes = zeros(Float64, nx)
    modes .= k * fftfreq(nx, nx)
    modes[1] = 1.0
    rhok = fft(rho) ./ modes
    rhok .*= -1im
    ifft!(rhok)
    return real(rhok)
end
