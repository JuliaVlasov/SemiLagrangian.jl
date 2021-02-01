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


end
