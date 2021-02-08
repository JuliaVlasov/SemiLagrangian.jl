var documenterSearchIndex = {"docs":
[{"location":"functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"Modules = [SemiLagrangian]\nOrder   = [:function]","category":"page"},{"location":"functions/#SemiLagrangian._getpolylagrange-Union{Tuple{N}, Tuple{Int64,Int64,Int64,N}} where N","page":"Functions","title":"SemiLagrangian._getpolylagrange","text":"_getpolylagrange(k::Int64, order::Int64, origin::Int64, N::DataType)\n\nFunction that return the k-th Lagrange Polynomial of a certain order. Coefficients are rational then the return is exact. The polynomial is equal to : prod_i=0 i neq k^order fracx - i - origink - i\n\nArguments\n\nk::Int64 : number of the Polynomial, k must be between 0 and order (0<= k <= order).\norder::Int64 : order of the polynomial.\norigin::Int64 : origin of the first indice.\nN::DataType : type of the Integer that is the base of the rational type, in fact Int64 or BigInt. BigInt is needed for order greater or equal to 21.\n\nReturns\n\nPolynomial{Rational{N}} : the k-th Lagrange polynomial of order order\n\nThrows\n\nDommaineError : when 0 <= k <= order is false or when N ∉ {BInt64, BigInt}\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.advection!-Union{Tuple{AdvectionData{T,Nsp,Nv,Nsum,timeopt}}, Tuple{timeopt}, Tuple{Nsum}, Tuple{Nv}, Tuple{Nsp}, Tuple{T}} where timeopt where Nsum where Nv where Nsp where T","page":"Functions","title":"SemiLagrangian.advection!","text":"advection!(self::AdvectionData)\n\nAdvection function of a multidimensional function f discretized on mesh\n\nArgument\n\nself::AdvectionData : mutable structure of variables data\n\nReturn value\n\ntrue : means that the advection series must continue\nfalse : means that the advection series is ended.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_charge!-Union{Tuple{AdvectionData{T,Nsp,Nv,Nsum,timeopt} where timeopt}, Tuple{Nsum}, Tuple{Nv}, Tuple{Nsp}, Tuple{T}} where Nsum where Nv where Nsp where T","page":"Functions","title":"SemiLagrangian.compute_charge!","text":"compute_charge!( self::AdvectionData)\n\nCompute charge density\n\nρ(x,t) = ∫ f(x,v,t) dv\n\nArgument\n\nself::AdvectionData : mutable structure of variables data.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_charge!-Union{Tuple{Nsum}, Tuple{Nv}, Tuple{Nsp}, Tuple{T}, Tuple{Array{T,Nsp},Tuple{Vararg{UniformMesh{T},Nv}},Array{T,Nsum}}} where Nsum where Nv where Nsp where T","page":"Functions","title":"SemiLagrangian.compute_charge!","text":"compute_charge!( rho, mesh_v, fvx)\n\nCompute charge density\n\nρ(x,t) = ∫ f(x,v,t) dv\n\nArguments\n\nrho::Array{T,Nsp} : result table correctly sized, it is he output\nt_mesh_v::NTuple{Nv,UniformMesh{T}} : velocity meshes\nf::Array{T,Nsum} : input data\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_ee-Tuple{AdvectionData}","page":"Functions","title":"SemiLagrangian.compute_ee","text":"compute_ee(self::AdvectionData)\n\ncompute electric enegie || E(t,.) ||_L2\n\nArgument\n\nself::AdvectionData : veriable advection data structure.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_ee-Union{Tuple{N}, Tuple{T}, Tuple{Tuple{Vararg{UniformMesh{T},N}},Tuple{Vararg{Array{T,N},N}}}} where N where T","page":"Functions","title":"SemiLagrangian.compute_ee","text":"compute_ee(t_mesh_sp, t_elf)\n\ncompute electric enegie || E(t,.) ||_L2\n\nArguments\n\nt_mesh_sp::NTuple{N,UniformMesh{T}} : tuple of space meshes\nt_elf::NTuple{N,Array{T,N}} : tuple of electric field\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_elfield!-Union{Tuple{AdvectionData{T,Nsp,Nv,Nsum,timeopt} where timeopt}, Tuple{Nsum}, Tuple{Nv}, Tuple{Nsp}, Tuple{T}} where Nsum where Nv where Nsp where T","page":"Functions","title":"SemiLagrangian.compute_elfield!","text":"compute_elfield!( self:AdvectionData)\n\ncomputation of electric field     ∇.e = - ρ\n\nArgument\n\nself::AdvectionData : mutable structure of variables data.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_ke-Union{Tuple{AdvectionData{T,Nsp,Nv,Nsum,timeopt} where timeopt}, Tuple{Nsum}, Tuple{Nv}, Tuple{Nsp}, Tuple{T}} where Nsum where Nv where Nsp where T","page":"Functions","title":"SemiLagrangian.compute_ke","text":"compute_ke(self::AdvectionData)\n\nkinetic Energie\n\n1/2∫∫ v^2 f(x,v,t) dv dx\n\nArguments\n\nself::AdvectionData : mutable structure of variables data.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_ke-Union{Tuple{Nsum}, Tuple{Nv}, Tuple{Nsp}, Tuple{T}, Tuple{Tuple{Vararg{UniformMesh{T},Nsp}},Tuple{Vararg{UniformMesh{T},Nv}},Array{T,Nsum}}} where Nsum where Nv where Nsp where T","page":"Functions","title":"SemiLagrangian.compute_ke","text":"compute_ke(t_mesh_sp, t_mesh_v, f)\n\nkinetic Energie\n\n1/2∫∫ v^2 f(x,v,t) dv dx\n\nArguments\n\nt_mesh_sp::NTuple{Nsp, UniformMesh{T}} : tuple of space meshes\nt_mesh_v::NTuple{Nv, UniformMesh{T}} : tuple of velocity meshes\nf::Array{T,Nsum} : function data.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.getalpha-Tuple{SemiLagrangian.RotationVar,AdvectionData,Any}","page":"Functions","title":"SemiLagrangian.getalpha","text":"getalpha(self::AdvectionData{T, Nsp, Nv, Nsum}, ind)\n\nImplementation of the interface function that is called before each interpolation in advection\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.interpolate!-Union{Tuple{order}, Tuple{T}, Tuple{AbstractArray{T,1},AbstractArray{T,1},Int64,Array{T,1},InterpolationType{T,true,order}}, Tuple{AbstractArray{T,1},AbstractArray{T,1},Int64,Array{T,1},InterpolationType{T,true,order},Any}} where order where T","page":"Functions","title":"SemiLagrangian.interpolate!","text":"interpolate!( fp, fi, decint, precal, interp)\n\napply an offset to the function fi interpolate by interp struct, the result is in fp vector, decint and precal are precompute with get_precal method.\n\nArguments\n\nfp::AbstractVector : output vector\nfi::AbstractVector : input vector\ndecint : offset in units of dx\nprecal::Vector : vector of length order+1 precompute with get_precal(interp, dec) (dec is the offset)\ninterp::InterpolationType : interpolation implementation\n\nReturns :\n\nNo return\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.interpolate!-Union{Tuple{order}, Tuple{iscirc}, Tuple{T}, Tuple{Any,Any,Any,InterpolationType{T,iscirc,order}}} where order where iscirc where T","page":"Functions","title":"SemiLagrangian.interpolate!","text":"interpolate!( fp, fi, dec, interp)\n\napply the offset dec to the function fi interpolate by interp struct, the result is in fp Vector\n\nArguments\n\nfp : output vector of length n\nfi : input vector of length n\ndec : offset in units of dx\ninterp::InterpolationType : interpolation implementation\n\nReturns :\n\nNo return\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.nextstate!-Union{Tuple{AdvectionData{T,Nsp,Nv,Nsum,timeopt} where timeopt}, Tuple{Nsum}, Tuple{Nv}, Tuple{Nsp}, Tuple{T}} where Nsum where Nv where Nsp where T","page":"Functions","title":"SemiLagrangian.nextstate!","text":"nextstate!(self::AdvectionData{T, Nsp, Nv, Nsum})\n\nFunction called at the end of advection function to update internal state of AdvectionData structure\n\nArgument\n\nself::AdvectionData{T, Nsp, Nv, Nsum} : object to update\n\nreturn value\n\nret::Bool : true if the series must continue               false at the end of the series.\n\n\n\n\n\n","category":"method"},{"location":"contents/#Contents","page":"Contents","title":"Contents","text":"","category":"section"},{"location":"contents/","page":"Contents","title":"Contents","text":"","category":"page"},{"location":"contents/#Index","page":"Contents","title":"Index","text":"","category":"section"},{"location":"contents/","page":"Contents","title":"Contents","text":"","category":"page"},{"location":"quickstart/#Quickstart","page":"Quickstart","title":"Quickstart","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"First simulation using semi-lagrangian method:","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"\nusing SemiLagrangian\n\n","category":"page"},{"location":"#SemiLagrangian.jl-Documentation","page":"Home","title":"SemiLagrangian.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Let us consider an abstract scalar advection equation of the form","category":"page"},{"location":"","page":"Home","title":"Home","text":"$ \\frac{∂f}{∂t}+ a(x, t) ⋅ ∇f = 0. $","category":"page"},{"location":"","page":"Home","title":"Home","text":"The characteristic curves associated to this equation are the solutions of  the ordinary differential equations","category":"page"},{"location":"","page":"Home","title":"Home","text":"$ \\frac{dX}{dt} = a(X(t), t) $","category":"page"},{"location":"","page":"Home","title":"Home","text":"We shall denote by X(t x s) the unique solution of this equation  associated to the initial condition X(s) = x.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The classical semi-Lagrangian method is based on a backtracking of  characteristics. Two steps are needed to update the distribution function  f^n+1 at t^n+1 from its value f^n at time t^n :","category":"page"},{"location":"","page":"Home","title":"Home","text":"For each grid point x_i compute X(t^n x_i t^n+1) the value  of the characteristic at t^n which takes the value x_i at  t^n+1.\nAs the distribution solution of first equation verifies f^n+1(x_i) = f^n(X(t^n x_i t^n+1)) we obtain the desired value of f^n+1(x_i) by computing  f^n(X(t^nx_it^n+1) by interpolation as X(t^n x_i t^n+1)  is in general not a grid point.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Eric Sonnendrücker - Numerical methods for the Vlasov equations","category":"page"},{"location":"types/#Types","page":"Types","title":"Types","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"Modules = [SemiLagrangian]\nOrder   = [:type]","category":"page"},{"location":"types/#SemiLagrangian.Advection","page":"Types","title":"SemiLagrangian.Advection","text":"Advection{T, Nsp, Nv, Nsum, timeopt}\nAdvection(\n    t_mesh_sp::NTuple{Nsp, UniformMesh{T}},\n    t_mesh_v::NTuple{Nv, UniformMesh{T}},\n    t_interp_sp::NTuple{Nsp, InterpolationType{T}},\n    t_interp_v::NTuple{Nv, InterpolationType{T}},\n    dt_base::T;\n    tab_coef=[1//2, 1//1, 1//2],\n    tab_fct=[identity, identity, identity],\n    timeopt::TimeOptimization=NoTimeOpt\n)\n\nImmutable structure that contains constant parameters for multidimensional advection\n\nType parameters\n\nT::DataType : type of data\nNsp : number of space dimensions\nNv : number of velocity dimensions\nNsum : the total number of dimensions (Nsum = Nsp + Nv)\ntimeopt::TimeOptimization : time optimization\n\nArguments\n\nt_mesh_sp::NTuple{Nsp, UniformMesh{T}} : tuple of space meshes (one per space dimension)\nt_mesh_v::NTuple{Nv, UniformMesh{T}} : tuple of velocity meshes (one per velocity dimension)\nt_interp_sp::NTuple{Nsp, InterpolationType{T}} : tuple of space interpolations (one per space dimension)\nt_interp_v::NTuple{Nv, InterpolationType{T}} : tuple of velocity interpolations(one per velocity dimension)\ndt_base::T : time delta for one advection series\n\nKeywords\n\ntab_coef=[1//2, 1//1, 1//2] : coefficient table for one advection series, the   coefficients at odd indexes is for space advection series, the coefficients at even indexes is for velocity advection series\ntab_fct=[identity, identity, identity] : function table for one advection series, with the same indexes than tab_coef\n\nImplementation\n\nsizeall : tuple of the sizes of all dimensions (space before velocity)\nsizeitr : tuple of iterators of indexes of each dimension\nt_mesh_sp : tuple of space meshes\nt_mesh_v : tuple of velocity meshes\nt_interp_sp : tuple of space interpolation types\nt_interp_v : tuple of velocity interpolation types\ndt_base::T : time unit of an advection series\ntab_coef : coefficient table\ntab_fct : function table\nv_square : precompute for ke\nnbsplit : number of slices for split\nmpiid : MPI id\n\nThrows\n\nArgumentError : Nsp must be less or equal to Nv.\n\n\n\n\n\n","category":"type"},{"location":"types/#SemiLagrangian.AdvectionData","page":"Types","title":"SemiLagrangian.AdvectionData","text":"AdvectionData{T,Nsp,Nv,Nsum,timeopt}\nAdvectionData(\nadv::Advection{T,Nsp,Nv,Nsum,timeopt}, \ndata::Array{T,Nsum},\nparext)\n\nMutable structure that contains variable parameters of advection series\n\nType parameters\n\nT::DataType : type of data\nNsp : number of space dimensions\nNv : number of velocity dimensions\nNsum : the total number of dimensions (Nsum = Nsp + Nv)\ntimeopt::TimeOptimization : time optimization\n\nArguments\n\nadv::Advection{T,Nsp,Nv,Nsum} : link to the constant data of this advection\ndata::Array{T,Nsum} : Initial data of this advection\nparext : external data of this advection to compute alpha of each interpolations\n\nImplementation\n\nadv::Advection{T,Nsp,Nv,Nsum,timeopt} : link to the constant data of this advection\nstate_coef::Int : state that is the index of tab_coef, it is from one to lenth(tab_coef)\nstate_dim::Int : the dimension index, from 1 to Nsp in space states, from one to Nv in velocity state\ndata:Array{T,Nsum} : it is the working buffer\nbufdata::Vector{T} : vector of the same size of the working buffer\nt_buf::NTuple{Nsum, Array{T,2}} : tuple of buffer that is used to get the linear data for interpolation, one buffer per thread\ncache_alpha::Union{T,Nothing} : cache for precal, the precal is compute only when the alpha or decint values change\ncache_decint::Int64 : for precal cache\ncache_precal::Vector{T} : for precal cache\nparext::ExtDataAdv : external data of this advection to compute alpha of each interpolations\n\nMethods to define\n\ninitcoef!(parext::AbstractExtDataAdv, self::AdvectionData) : this method called at the beginning of each advection to initialize parext data. The self.parext mutable structure is the only data that init! can modify otherwise it leads to unpredictable behaviour.\ngetalpha(parext::AbstractExtDataAdv, self::AdvectionData, ind) : return the alpha number that is used for interpolation.\ngetperm(parext::AbstractExtDataAdv, advd::AdvectionData) : get the permutation of the dimension as a function of the current state\n\n\n\n\n\n","category":"type"},{"location":"types/#SemiLagrangian.B_Spline","page":"Types","title":"SemiLagrangian.B_Spline","text":"B_Spline{T,iscirc, order} <: InterpolationType{T,iscirc, order}\n\nAbstract supertype for all bspline interpolation type\n\n\n\n\n\n","category":"type"},{"location":"types/#SemiLagrangian.InterpolationType","page":"Types","title":"SemiLagrangian.InterpolationType","text":"InterpolationType{T, iscirc, order}\n\nAbstract supertype for all interpolation type\n\n\n\n\n\n","category":"type"},{"location":"types/#SemiLagrangian.Lagrange","page":"Types","title":"SemiLagrangian.Lagrange","text":"Lagrange{T, iscirc, order, N} <: InterpolationType{T, iscirc, order}\nLagrange(order, T::DataType=Float64; iscirc::Bool=true)\n\nType containing Lagrange Polynomials coefficients for Lagrange interpolation\n\nType parameters\n\nT : the type of data that is interpolate\niscirc::Bool : true if function is circular\norder::Int: order of lagrange interpolation\nN : type of integer, in fact Int64 or BigInt that is used to store lagrange polynomial\n\nImplementation :\n\nfact_order::N : factorial of the order\nlagpol:Vector{Polynomial{N}} : vector of all lagrange polynomial, per example the k-th Lagrange polynomial for the designed order is lagpol[k+1]/fact_order\n\nArguments :\n\norder : the order of interpolation\n[T::DataType=Float64] : The type values to interpolate \n\n\n\n\n\n","category":"type"},{"location":"types/#SemiLagrangian.LuSpline","page":"Types","title":"SemiLagrangian.LuSpline","text":"struct LuSpline{T}\n\nStructure of a LU decomposition of circular banded matrix, a LU decomposition can be stored in a Matrix which is equal to L + U - I. For a circular Banded matrix all non zero coefficients are in the band and in the last columns and lines\n\nImplementation\n\nband::Matrix{T} : matrix of size (kl+ku+1, n-kl)\nku : size of the band upper the diagonal\nkl : size of the band lower the diagonal\niscirc : true if and only if original matrix is circular\nisLU : true if LU decomposition has been perform\nlastcols : only in circular case, Matrix of size (n-ku, kl) that represents teh last columns of matrix\nlastrows : only in circular case, Matrix of size (n, ku) that represents the last rows of matrix\n\n\n\n\n\n","category":"type"},{"location":"types/#SemiLagrangian.PoissonConst","page":"Types","title":"SemiLagrangian.PoissonConst","text":"PoissonConst{T, Nsp, Nv}\nPoissonConst(adv::Advection{T, Nsp, Nv, Nsum}; isfftbig=true)\n\nConstant data for the computation of poisson coefficients\n\nArguments\n\nadv::Advection{T, Nsp, Nv, Nsum, timeopt} : Advection constant data\nisfftbig=true: if true compute the fttbig structure\n\nImplementation\n\nadv : Advection constant data\n`v_k' : vector of vector of fourier coefficents for integration for each space dimension\nfctv_k : Array of space dimensions of the inverse of the norm of fourier coefficients\npfftbig : Fourier data for space dimensions\n\n\n\n\n\n","category":"type"},{"location":"types/#SemiLagrangian.PoissonVar","page":"Types","title":"SemiLagrangian.PoissonVar","text":"PoissonVar{T, Nsp, Nv} <: AbstractExtDataAdv\nPoissonVar(pc::PoissonConst{T, Nsp, Nv})\n\nmutable structure of variable data for the poisson computation\n\nArguments\n\npc::PoissonConst{T, Nsp, Nv} : poisson constant data\n\nImplementation\n\npc::PoissonConst{T, Nsp, Nv} : poisson constant data\nrho::Array{T, Nsp} : result of the compute_charge that is the sum along velocity dimensions\nt_elfield::NTuple{Nsp,Array{Complex{T}, Nsp}} : electric fields initialized at each beginning of velocity advection subseries\n\n\n\n\n\n","category":"type"},{"location":"types/#SemiLagrangian.UniformMesh","page":"Types","title":"SemiLagrangian.UniformMesh","text":"UniformMesh(start, stop, length)\n\n1D uniform mesh data.\n\nlength   : number of points length-1 : number of cells\n\nIf you want remove the last point for periodic domain, set endpoint=false\n\nArguments\n\n- `step`  : size step\n- `points`: Array with node positions\n- `width` : Distance between left and right edges.\n\n\n\n\n\n","category":"type"}]
}
