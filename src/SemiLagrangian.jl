"""
$(EXPORTS)

$(IMPORTS)

$(README)

"""
module SemiLagrangian

using DocStringExtensions
using Polynomials: StandardBasisPolynomial
using Polynomials
using FFTW
using Base.Threads
using Requires


include("util.jl")
include("cplxlagrange.jl")
include("fftbig.jl")
include("mesh.jl")
include("interpolation.jl")

function __init__()
    @require MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195" include("mpiinterface.jl")
    @require MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195" include("mpiinterpolation.jl")
end

include("lagrange.jl")
include("hermite.jl")
include("spline.jl")
include("bspline.jl")
include("bsplinelu.jl")
include("bsplinefft.jl")
include("splitting.jl")
include("advection.jl")
include("util_poisson.jl")
include("poisson.jl")
include("rotation.jl")
include("translation.jl")
include("quasigeostrophic.jl")

export UniformMesh, start, stop, AbstractInterpolation, get_order
export Advection, AdvectionData
export nosplit, standardsplit, strangsplit, triplejumpsplit, order6split, hamsplit_3_11, table2split
export Lagrange, Hermite, BSplineLU, BSplineFFT, interpolate!
export compute_charge!,
    compute_elfield!, compute_elfield, compute_ee, compute_ke, advection!
export dotprod, getpoissonvar, getrotationvar, gettranslationvar
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

export TypePoisson, StdPoisson, StdPoisson2d, StdABp

export get_type, sizeall, getdata, OpTuple, getgeovar, initdata!, getenergyall

end
