```@meta
CurrentModule = SemiLagrangian
```

# B-spline

Both interpolants use periodic boundaries on a uniform grid. Their constructors
take the interpolation order and the grid length, followed by an optional
numeric type (defaulting to `Float64`). The grid excludes the repeated periodic
endpoint. In both examples, the shift is measured in grid cells.

## BSplineLU

`BSplineLU(order, n)` computes the spline coefficients using an LU-factorized
periodic system. Both odd and even orders are supported.

```@example bspline_lu
using SemiLagrangian

n = 128
x = (0:n-1) ./ n
u = sin.(2pi .* x)
shift = 0.3
expected = sin.(2pi .* (x .+ shift / n))

interpolant = BSplineLU(4, n)
output = similar(u)
interpolate!(output, u, shift, interpolant)

println("Maximum error: ", maximum(abs.(output .- expected)))
```

```@docs
BSplineLU
```

## BSplineFFT

`BSplineFFT(order, n)` uses FFTs to solve the periodic spline system.

```@example bspline_fft
using SemiLagrangian

n = 128
x = (0:n-1) ./ n
u = sin.(2pi .* x)
shift = 0.3
expected = sin.(2pi .* (x .+ shift / n))

interpolant = BSplineFFT(5, n)
output = similar(u)
interpolate!(output, u, shift, interpolant)

println("Maximum error: ", maximum(abs.(output .- expected)))
```

```@docs
BSplineFFT
```
