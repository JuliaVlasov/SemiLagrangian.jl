"""
$(README)

# Dependencies

$(IMPORTS)

"""
module SemiLagrangian

using DocStringExtensions
using FFTW
using Base.Threads
using Requires


include("util.jl")
include("polynomial.jl")
include("cplxlagrange.jl")
include("fftbig.jl")
include("mesh.jl")
include("interpolation.jl")

function __init__()
    @require MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195" include("mpiinterface.jl")
    @require MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195" include("mpiinterpolation.jl")
end

include("lagrange.jl")
include("lagrange_interp.jl")
include("hermite.jl")
include("hermite_interp.jl")
include("spline.jl")
include("spline_interp.jl")
include("bspline.jl")
include("bsplinelu.jl")
include("bsplinefft.jl")
include("splitting.jl")
include("advection.jl")
include("rotation.jl")
include("translation.jl")
include("periodic_interpolation/bsplines.jl")
include("periodic_interpolation/spectral.jl")
include("periodic_interpolation/fast_lagrange.jl")
include("periodic_interpolation/lagrange.jl")
include("periodic_interpolation/cubic_interp.jl")

export UniformMesh, start, stop, AbstractInterpolation, get_order
export Advection, AdvectionData
export nosplit, standardsplit, strangsplit, triplejumpsplit, order6split, hamsplit_3_11, table2split
export Lagrange, Hermite, BSplineLU, BSplineFFT, interpolate!
export LagrangeCache, lagrange_interp, lagrange_interp!
export HermiteCache, hermite_interp, hermite_interp!
export SplineCache, spline_interp, spline_interp!
export Polynomial, derivative, degree
export advection!
export dotprod, getrotationvar, gettranslationvar
export TimeOptimization,
    NoTimeOpt,
    SimpleThreadsOpt,
    SplitThreadsOpt,
    MPIOpt,
    TimeAlgorithm,
    NoTimeAlg,
    ABTimeAlg_ip,
    ABTimeAlg_init,
    ABTimeAlg_new

export sizeall, getdata, OpTuple
export Spectral, FastLagrange, CubicSpline, PeriodicLagrange, PeriodicBSpline

end
