module SemiLagrangian

# import Base.Threads: @sync, @spawn, nthreads, threadid



using Polynomials: StandardBasisPolynomial
using Polynomials
using FFTW
using Base.Threads

# include("mpoly/MultiPoly.jl")
include("mpiinterface.jl")
include("clockobs.jl")

include("util.jl")
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
include("advection.jl")
include("util_poisson.jl")
include("poisson.jl")
include("rotation.jl")
include("translation.jl")
include("quasigeostrophic.jl")

export UniformMesh, start, stop, AbstractInterpolation, get_order
export Advection, AdvectionData
export nosplit, standardsplit, strangsplit, triplejumpsplit, order6split, hamsplit_3_11, ymsplit, table2split
export Lagrange, Hermite, B_SplineLU, B_SplineFFT, interpolate!
export compute_charge!,
    compute_elfield!, compute_elfield, compute_ee, compute_ke, advection!
export dotprod, getpoissonvar, getrotationvar, gettranslationvar
export TimeOptimization, NoTimeOpt, SimpleThreadsOpt, SplitThreadsOpt, MPIOpt, TimeAlgorithm, NoTimeAlg, ABTimeAlg
export TypePoisson, StdPoisson, StdPoisson2d, StdOrder2_1, StdOrder2_2, StdAB, StdAB2, StdRK4, StdABinit, StdABp
export get_type, sizeall, getdata, OpTuple, getgeovar

end
