module SemiLagrangian

# import Base.Threads: @sync, @spawn, nthreads, threadid



using Polynomials
using FFTW
using Base.Threads

include("nompiinterface.jl")

include("util.jl")
include("fftbig.jl")
include("mesh.jl")
include("interpolation.jl")
include("lagrange.jl")
include("spline.jl")
include("bspline.jl")
include("bsplinelu.jl")
include("bsplinefft.jl")
include("advection.jl")
include("util_poisson.jl")
include("poisson.jl")
include("rotation.jl")
include("translation.jl")

export UniformMesh, Advection, AdvectionData, AbstractInterpolation, AbstractInterpolation2d
export Lagrange, B_SplineLU, B_SplineFFT, interpolate!
export compute_charge!, compute_elfield!, compute_ee, compute_ke, advection!
export dotprod, getpoissonvar, getrotationvar, gettranslationvar 
export TimeOptimization, NoTimeOpt, SimpleThreadsOpt, SplitThreadsOpt, MPIOpt 
export get_type, sizeall, getdata

end
