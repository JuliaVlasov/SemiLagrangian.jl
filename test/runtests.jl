using LinearAlgebra: size
using Test
using SemiLagrangian

# include("mpoly/runtests.jl")
include("test_interpolation.jl")
include("test_lagrange.jl")
include("test_hermite.jl")
include("test_splinelu.jl")
include("test_bspline.jl")
include("testfftbig.jl")
include("test_util.jl")
include("test_mesh.jl")
include("test_advection.jl")
# include("test_poisson1d.jl")
include("test_poisson.jl")
include("test_poisson2d.jl")
# include("test_rotation1d.jl")
include("test_rotation.jl")
# include("test_translation1d.jl")
include("test_translation.jl")
include("test_swirling.jl")
include("test_quasigeostrophic.jl")
