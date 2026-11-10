"""
$(README)

# Dependencies

$(IMPORTS)

"""
module SemiLagrangian

using DocStringExtensions
using FFTW
using Base.Threads


include("util.jl")
include("polynomial.jl")
include("cplxlagrange.jl")
include("fftbig.jl")
include("mesh.jl")
include("interpolation.jl")
include("lagrange.jl")
include("hermite.jl")
include("spline.jl")
include("bspline.jl")
include("bsplinelu.jl")
include("bsplinefft.jl")
include("splitting.jl")
include("advection.jl")
include("spectral.jl")
include("fast_lagrange.jl")
include("lagrange_fft.jl")
include("cubic_interp.jl")

export UniformMesh, start, stop, AbstractInterpolation, get_order
export AbstractBoundaryCondition, PeriodicBC, ExtrapolationBC
export nosplit, standardsplit, strangsplit, sheartransport
export triplejumpsplit, order6split, hamsplit_3_11, table2split
export Lagrange, Hermite, BSplineLU, BSplineFFT, interpolate!
export Polynomial, derivative, degree
export advection!
export AdvectionOperator, split_advection!
export dotprod

export OpTuple
export Spectral, LocalLagrange, CubicSpline, LagrangeFFT

end
