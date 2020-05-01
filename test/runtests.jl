using Test
using LinearAlgebra
using SemiLagrangian

include("test_lagrange_1d.jl")
include("test_vp_1d1v.jl")
#include("test_vp_2d2v.jl")
include("test_bspline_periodic.jl")
include("test_banded_matrix.jl")
include("test_spline_1d.jl")
include("test_spline_interpolator_1d.jl")
include("test_splinepp.jl")
include("test_splinenn.jl")
include("test_rotation_2d.jl")
