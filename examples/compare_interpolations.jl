using FastInterpolations
using SemiLagrangian

n = 128
alpha = 0.25
xi = 1.0:n
u = sin.(2π .* (0:n-1) ./ n) .+ 0.2 .* cos.(6π .* (0:n-1) ./ n)
expected = sin.(2π .* ((0:n-1) .+ alpha) ./ n) .+
    0.2 .* cos.(6π .* ((0:n-1) .+ alpha) ./ n)
result = zeros(n)

methods = [
    "Spectral" => (out -> interpolate!(out, Spectral(n), u, alpha)),
    "Lagrange FFT" => (out -> interpolate!(out, LagrangeFFT(n, 7), u, alpha)),
    "Fast Lagrange" => (out -> interpolate!(out, FastLagrange(7), u, alpha)),
    "FastInterpolations cubic" => (
        out -> (out .= cubic_interp(
            xi,
            u,
            xi .+ alpha;
            bc=FastInterpolations.PeriodicBC(endpoint=:exclusive),
        ))
    ),
]

for (name, interpolate) in methods
    interpolate(result)
    println(rpad(name, 28), " max error = ", maximum(abs, result .- expected))
end
