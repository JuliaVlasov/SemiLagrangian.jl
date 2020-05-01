module SemiLagrangian

import Base.Threads: @sync, @spawn, nthreads, threadid

abstract type InterpolationType end
abstract type AbstractAdvection end

include("mesh.jl")
include("geometry.jl")
include("bspline_periodic.jl")
include("splinenn.jl")
include("splinepp.jl")
include("banded_matrix.jl")
include("spline_1d.jl")
include("spline_interpolator_1d.jl")
include("lagrange_interpolation.jl")
include("advection.jl")

end
