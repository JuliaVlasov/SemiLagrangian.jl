var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Documentation",
    "title": "Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "#SemiLagrangian.jl-Documentation-1",
    "page": "Documentation",
    "title": "SemiLagrangian.jl Documentation",
    "category": "section",
    "text": "Modules = [SemiLagrangian] Order   = [:type,:function]"
},

{
    "location": "advections/#",
    "page": "Advections",
    "title": "Advections",
    "category": "page",
    "text": ""
},

{
    "location": "advections/#SemiLagrangian.interpolate-Tuple{Array{Float64,1},Int64,Float64,Float64,Float64}",
    "page": "Advections",
    "title": "SemiLagrangian.interpolate",
    "category": "method",
    "text": "interpolate( coeffs, n, f)\n\nCompute interpolatted value at x from coefficients get  from compute_interpolants.\n\nThis function is a julia version of a Fortran code written by\n\nEdwin Chacon Golcher (Institute of Physics of the Czech Academy of Sciences)\n\n\n\n\n\n"
},

{
    "location": "advections/#SemiLagrangian.interpolate-Tuple{Int64,Array{Complex{Float64},1},Float64,Float64}",
    "page": "Advections",
    "title": "SemiLagrangian.interpolate",
    "category": "method",
    "text": "interpolate( p, f, delta, alpha)\n\nCompute the interpolating spline of degree p of odd degree of a 1D function f on a periodic uniform mesh, at all points x-alpha input f is the Fourier transform and complex you have to in2erse transform after return\n\n\n\n\n\n"
},

{
    "location": "advections/#SemiLagrangian.interpolate-Tuple{Int64,Array{Float64,1},Float64,Float64}",
    "page": "Advections",
    "title": "SemiLagrangian.interpolate",
    "category": "method",
    "text": "interpolate( p, f, delta, alpha)\n\nCompute the interpolating spline of degree p of odd degree of a 1D function f on a periodic uniform mesh, at all points x-alpha. f type is Vector{Float64}.\n\nDerived from a Python program written by  Eric Sonnendrucker (Max-Planck-Institut fur Plasmaphysik - Garching)\n\n\n\n\n\n"
},

{
    "location": "advections/#SemiLagrangian.advection!-Tuple{Array{Complex{Float64},2},VlasovBase.UniformMesh,Array{Float64,1},Int64,Float64,BSpline}",
    "page": "Advections",
    "title": "SemiLagrangian.advection!",
    "category": "method",
    "text": "advection!(f, mesh, v, n2, dt, interp)\n\nAdvection of a 2d function f along its first dimension with velocity v. Since the fft are computed inplace, the function  must be represented by a Array{Complex{Float64},2}.\n\n\n\n\n\n"
},

{
    "location": "advections/#SemiLagrangian.advection!-Tuple{Array{Float64,2},VlasovBase.UniformMesh,Array{Float64,1},Float64,CubicSpline,Int64}",
    "page": "Advections",
    "title": "SemiLagrangian.advection!",
    "category": "method",
    "text": " advection!( f, mesh1, v,  dt)\n\nSemi-lagrangian advection function of 2D distribution function represented  by array f. The advection operates along axis (=1 is most efficient)  with speed v during dt.\n\nIt uses cubic splines interpolation.\n\n\n\n\n\n"
},

{
    "location": "advections/#SemiLagrangian.advection!-Tuple{Array{Float64,2},VlasovBase.UniformMesh,VlasovBase.UniformMesh,Array{Float64,1},Float64,BSpline,Int64}",
    "page": "Advections",
    "title": "SemiLagrangian.advection!",
    "category": "method",
    "text": "advection!( mesh, f, v, dt, interp, axis)\n\nAdvection of a 2d function f discretized on a 2d mesh along the input axis at velocity v\n\n\n\n\n\n"
},

{
    "location": "advections/#SemiLagrangian.bspline-Tuple{Int64,Int64,Float64}",
    "page": "Advections",
    "title": "SemiLagrangian.bspline",
    "category": "method",
    "text": "bspline(p, j, x)\n\nReturn the value at x in [0,1[ of the B-spline with integer nodes of degree p with support starting at j. Implemented recursively using the de Boor\'s recursion formula\n\nDerived from a Python program written by  Eric Sonnendrucker (Max-Planck-Institut fur Plasmaphysik - Garching))\n\n\n\n\n\n"
},

{
    "location": "advections/#SemiLagrangian.compute_interpolants-Tuple{Int64,Array{Float64,N} where N}",
    "page": "Advections",
    "title": "SemiLagrangian.compute_interpolants",
    "category": "method",
    "text": "compute_interpolants( n, f)\n\nCompute interpolation coefficients\n\nThis function is a julia version of a Fortran code written by Edwin Chacon Golcher (Institute of Physics of the Czech Academy of Sciences)\n\n\n\n\n\n"
},

{
    "location": "advections/#SemiLagrangian.eval_basis!-Tuple{Any,Any,Any}",
    "page": "Advections",
    "title": "SemiLagrangian.eval_basis!",
    "category": "method",
    "text": "Evaluate value at x of all basis functions with support in local cell values[j] = B_j(x) for jmin <= j <= jmin+degree\n\n\n\n\n\n"
},

{
    "location": "advections/#SemiLagrangian.eval_deriv!-Tuple{Any,Any,Any}",
    "page": "Advections",
    "title": "SemiLagrangian.eval_deriv!",
    "category": "method",
    "text": "Evaluate derivative at x of all basis functions with support in local cell derivs[j] = B_j\'(x) for jmin <= j <= jmin+degree\n\n\n\n\n\n"
},

{
    "location": "advections/#SemiLagrangian.eval_deriv-Tuple{Any,Any}",
    "page": "Advections",
    "title": "SemiLagrangian.eval_deriv",
    "category": "method",
    "text": "Evaluate derivative of 1D spline at location x: y=S\'(x)\n\n\n\n\n\n"
},

{
    "location": "advections/#SemiLagrangian.eval_value-Tuple{Any,Any}",
    "page": "Advections",
    "title": "SemiLagrangian.eval_value",
    "category": "method",
    "text": "Evaluate value of 1D spline at location x: y=S(x)\n\n\n\n\n\n"
},

{
    "location": "advections/#SemiLagrangian.interpolate!-Tuple{SemiLagrangian.SplineNN,Array{Float64,2},Float64,Float64}",
    "page": "Advections",
    "title": "SemiLagrangian.interpolate!",
    "category": "method",
    "text": "interpolation par spline periodique dans les deux directions. Les points d\'interpolation sont definis grace a depx et depy qui definissent le deplacement par rapport au maillage.\n\nf contient les valeurs de la fonction de distribution\ndepx et depy : deplacements par rapport au maillage   des points dans les quels on veut evaluer la spline.\n\n\n\n\n\n"
},

{
    "location": "advections/#SemiLagrangian.nat_x!-Tuple{SemiLagrangian.SplineNN,Array{Float64,2}}",
    "page": "Advections",
    "title": "SemiLagrangian.nat_x!",
    "category": "method",
    "text": "natural splines\n\n\n\n\n\n"
},

{
    "location": "advections/#Advection-functions-1",
    "page": "Advections",
    "title": "Advection functions",
    "category": "section",
    "text": "CurrentModule = SemiLagrangianModules = [SemiLagrangian]\nOrder   = [:type,:function]"
},

{
    "location": "bsl/#",
    "page": "BSL",
    "title": "BSL",
    "category": "page",
    "text": ""
},

{
    "location": "bsl/#Semi-Lagrangian-method-1",
    "page": "BSL",
    "title": "Semi-Lagrangian method",
    "category": "section",
    "text": "Let us consider an abstract scalar advection equation of the form$ \\frac{∂f}{∂t}+ a(x, t) ⋅ ∇f = 0. $The characteristic curves associated to this equation are the solutions of  the ordinary differential equations$ \\frac{dX}{dt} = a(X(t), t) $We shall denote by X(t x s) the unique solution of this equation  associated to the initial condition X(s) = x.The classical semi-Lagrangian method is based on a backtracking of  characteristics. Two steps are needed to update the distribution function  f^n+1 at t^n+1 from its value f^n at time t^n :For each grid point x_i compute X(t^n x_i t^n+1) the value  of the characteristic at t^n which takes the value x_i at  t^n+1.\nAs the distribution solution of first equation verifies f^n+1(x_i) = f^n(X(t^n x_i t^n+1)) we obtain the desired value of f^n+1(x_i) by computing  f^n(X(t^nx_it^n+1) by interpolation as X(t^n x_i t^n+1)  is in general not a grid point.Eric Sonnendrücker - Numerical methods for the Vlasov equationsCurrentModule = SemiLagrangian"
},

{
    "location": "bsl/#Advection-functions-1",
    "page": "BSL",
    "title": "Advection functions",
    "category": "section",
    "text": "advection!(::Array{Float64,2}, ::UniformMesh, v, ::CubicSpline, ::Float64, ::Int64)advection!(::Array{Float64,2}, ::UniformMesh, v, ::BSpline, ::Float64, ::Int64)"
},

{
    "location": "contents/#",
    "page": "Contents",
    "title": "Contents",
    "category": "page",
    "text": ""
},

{
    "location": "contents/#Contents-1",
    "page": "Contents",
    "title": "Contents",
    "category": "section",
    "text": ""
},

{
    "location": "contents/#Index-1",
    "page": "Contents",
    "title": "Index",
    "category": "section",
    "text": ""
},

]}
