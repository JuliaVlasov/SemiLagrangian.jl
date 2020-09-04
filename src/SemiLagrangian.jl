# module SemiLagrangianBig

import Base.Threads: @sync, @spawn, nthreads, threadid


include("lapack.jl")
include("mesh.jl")
include("geometry.jl")
include("splinenn.jl")
include("splinepp.jl")
include("banded_matrix.jl")
include("spline_1d.jl")
include("spline_interpolator_1d.jl")
include("advection.jl")
include("lagrange.jl")
include("lagrange_interpolation.jl")
include("bspline_periodic.jl")


# end
