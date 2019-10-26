var documenterSearchIndex = {"docs":
[{"location":"functions/#Functions-and-types-1","page":"Functions","title":"Functions and types","text":"","category":"section"},{"location":"functions/#","page":"Functions","title":"Functions","text":"Modules = [SemiLagrangian]\nOrder   = [:type,:function]","category":"page"},{"location":"functions/#SemiLagrangian.Advection","page":"Functions","title":"SemiLagrangian.Advection","text":"Advection(interpolation_type, mesh, LBC, RBC)\n\nCreates a 1d backward semi-lagrangian advection.\n\ninterp   : Interpolation type (Bspline(degree), Lagrange(degree))\nmesh     : UniformMesh along advection direction\nLBC, RBC : Boundary conditions type (:periodic, :Hermite)\n\n\n\n\n\n","category":"type"},{"location":"functions/#SemiLagrangian.Advection-Tuple{Array{Float64,2},Array{Float64,1},Float64}","page":"Functions","title":"SemiLagrangian.Advection","text":"advection!(f, v, dt)\n\nAdvection of a 2d function f discretized on a 2d mesh along the input axis at velocity v. This function is created from Advector callable type.\n\n```julia mesh = UniformMesh( -π, π, 64 ) advection! = Advection( mesh, Bspline(3), :periodic )\n\nf = exp.( - mesh.points .^ 2 )\n\ndt = 0.5 v  = ones( Float64, mesh.length)\n\nadvection!( f, v, dt )\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.PeriodicAdvection","page":"Functions","title":"SemiLagrangian.PeriodicAdvection","text":"advection! = PeriodicAdvection( mesh, bspl )\n\nType to perform 1d advection on periodic domain. Bspline interpolation is used combined with fft transform to solve the system. Bspline degree must be odd.\n\n\np = 5\nnx = 128\n\nx1min, x1max = -10, 10\n\nmesh1 = UniformMesh(x1min, x1max, n1; endpoint=false)\n\nadvection! = PeriodicAdvection( mesh1, Bspline(p) )\n\nadvection! computes the interpolating spline of degree p of odd degree of first dimension of array f on a periodic uniform mesh, at all points x-alpha. f type is Array{Float64,2}.\n\n\n\n\n\n","category":"type"},{"location":"functions/#SemiLagrangian.UniformMesh","page":"Functions","title":"SemiLagrangian.UniformMesh","text":"UniformMesh(start, stop, length)\n\n1D uniform mesh data.\n\nlength   : number of points length-1 : number of cells\n\nIf you want remove the last point for periodic domain, set endpoint=false\n\n- `step`  : size step\n- `points`: Array with node positions\n- `width` : Distance between left and right edges.\n\n\n\n\n\n","category":"type"},{"location":"functions/#SemiLagrangian.compute_charge!-Tuple{Array{Float64,1},UniformMesh,Array{Float64,2}}","page":"Functions","title":"SemiLagrangian.compute_charge!","text":"compute_rho!( rho, mesh_v, fvx)\n\nCompute charge density\n\nρ(x,t) = ∫ f(x,v,t) dv\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_e!-Tuple{Array{Float64,1},UniformMesh,Array{Float64,1}}","page":"Functions","title":"SemiLagrangian.compute_e!","text":"compute_e!( e, mesh, ρ)\n\n∇.e = - ρ\n\nInplace computation of electric field. Fields e and rho are already allocated.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.bspline-Tuple{Int64,Int64,Float64}","page":"Functions","title":"SemiLagrangian.bspline","text":"bspline(p, j, x)\n\nReturn the value at x in [0,1[ of the B-spline with integer nodes of degree p with support starting at j. Implemented recursively using the de Boor's recursion formula using the De Boor's Algorithm\n\nB_i0(x) = left\nbeginmatrix\n1  mathrmif  quad t_i  x  t_i+1 \n0  mathrmotherwise\nendmatrix\nright\n\nB_ip(x) = fracx - t_it_i+p - t_i B_ip-1(x)\n+ fract_i+p+1 - xt_i+p+1 - t_i+1 B_i+1p-1(x)\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.eval_basis!-Tuple{Any,Any,Any}","page":"Functions","title":"SemiLagrangian.eval_basis!","text":"eval_basis!( spl, x, values )\n\nEvaluate value at x of all basis functions with support in local cell values[j] = B_j(x) for jmin <= j <= jmin+degree\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.eval_deriv!-Tuple{Any,Any,Any}","page":"Functions","title":"SemiLagrangian.eval_deriv!","text":"eval_deriv!( derivs, spl, x )\n\nEvaluate derivative at x of all basis functions with support in local cell derivs[j] = B_j'(x) for jmin <= j <= jmin+degree\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.eval_deriv-Tuple{Any,Any}","page":"Functions","title":"SemiLagrangian.eval_deriv","text":"eval_deriv( spl, x )\n\nEvaluate derivative of 1D spline at location x: y=S'(x)\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.eval_value-Tuple{Any,Any}","page":"Functions","title":"SemiLagrangian.eval_value","text":"eval_value( spl, x )\n\nEvaluate value of 1D spline at location x: y=S(x)\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.interpolate!-Tuple{SemiLagrangian.SplineNN,Array{Float64,2},Float64,Float64}","page":"Functions","title":"SemiLagrangian.interpolate!","text":"interpolation par spline periodique dans les deux directions. Les points d'interpolation sont definis grace a depx et depy qui definissent le deplacement par rapport au maillage.\n\nf contient les valeurs de la fonction de distribution\ndepx et depy : deplacements par rapport au maillage   des points dans les quels on veut evaluer la spline.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.nat_x!-Tuple{SemiLagrangian.SplineNN,Array{Float64,2}}","page":"Functions","title":"SemiLagrangian.nat_x!","text":"natural splines\n\n\n\n\n\n","category":"method"},{"location":"contents/#Contents-1","page":"Contents","title":"Contents","text":"","category":"section"},{"location":"contents/#","page":"Contents","title":"Contents","text":"","category":"page"},{"location":"contents/#Index-1","page":"Contents","title":"Index","text":"","category":"section"},{"location":"contents/#","page":"Contents","title":"Contents","text":"","category":"page"},{"location":"#SemiLagrangian.jl-Documentation-1","page":"Documentation","title":"SemiLagrangian.jl Documentation","text":"","category":"section"},{"location":"#","page":"Documentation","title":"Documentation","text":"Let us consider an abstract scalar advection equation of the form","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"$ \\frac{∂f}{∂t}+ a(x, t) ⋅ ∇f = 0. $","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"The characteristic curves associated to this equation are the solutions of  the ordinary differential equations","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"$ \\frac{dX}{dt} = a(X(t), t) $","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"We shall denote by X(t x s) the unique solution of this equation  associated to the initial condition X(s) = x.","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"The classical semi-Lagrangian method is based on a backtracking of  characteristics. Two steps are needed to update the distribution function  f^n+1 at t^n+1 from its value f^n at time t^n :","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"For each grid point x_i compute X(t^n x_i t^n+1) the value  of the characteristic at t^n which takes the value x_i at  t^n+1.\nAs the distribution solution of first equation verifies f^n+1(x_i) = f^n(X(t^n x_i t^n+1)) we obtain the desired value of f^n+1(x_i) by computing  f^n(X(t^nx_it^n+1) by interpolation as X(t^n x_i t^n+1)  is in general not a grid point.","category":"page"},{"location":"#","page":"Documentation","title":"Documentation","text":"Eric Sonnendrücker - Numerical methods for the Vlasov equations","category":"page"}]
}