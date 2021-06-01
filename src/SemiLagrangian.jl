module SemiLagrangian

# import Base.Threads: @sync, @spawn, nthreads, threadid



using Polynomials
using FFTW
using Base.Threads

include("mpiinterface.jl")
include("clockobs.jl")

include("util.jl")
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

export UniformMesh, start, stop, AbstractInterpolation, get_order
export Advection, AdvectionData
export standardsplit, strangsplit, triplejumpsplit, order6split, hamsplit_3_11
export Lagrange, Hermite, B_SplineLU, B_SplineFFT, interpolate!
export compute_charge!,
    compute_elfield!, compute_elfield, compute_ee, compute_ke, advection!
export dotprod, getpoissonvar, getrotationvar, gettranslationvar
export TimeOptimization, NoTimeOpt, SimpleThreadsOpt, SplitThreadsOpt, MPIOpt
export get_type, sizeall, getdata

end
