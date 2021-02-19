var documenterSearchIndex = {"docs":
[{"location":"functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"Modules = [SemiLagrangian]\nOrder   = [:function]","category":"page"},{"location":"functions/#Base.length-Tuple{UniformMesh}","page":"Functions","title":"Base.length","text":"Base.length(mesh::UniformMesh)\n\nGet the length of the mesh\n\nArgument\n\nmesh::UniformMesh : the mesh\n\nReturn\n\nlength : the length of the mesh that is the number of points and cells\n\n\n\n\n\n","category":"method"},{"location":"functions/#Base.step-Tuple{UniformMesh}","page":"Functions","title":"Base.step","text":"Base.step(mesh::UniformMesh)\n\nGet the step of the mesh\n\nArgument\n\nmesh::UniformMesh : the mesh\n\nReturn\n\nstep : the step of the mesh that is the difference between two contiguous points\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian._getpolylagrange-Union{Tuple{Int64,Int64,Int64}, Tuple{N}} where N<:Integer","page":"Functions","title":"SemiLagrangian._getpolylagrange","text":"_getpolylagrange(k::Int64, order::Int64, origin::Int64)\n\nFunction that return the k-th Lagrange Polynomial of a certain order. Coefficients are rational then the return is exact. The polynomial is equal to : prod_i=0 i neq k^order fracx - i - origink - i\n\nArguments\n\nk::Int64 : number of the Polynomial, k must be between 0 and order (0<= k <= order).\norder::Int64 : order of the polynomial.\norigin::Int64 : origin of the first indice.\n\nReturns\n\nPolynomial{Rational{BigInt}} : the k-th Lagrange polynomial of order order\n\nThrows\n\nDommaineError : when 0 <= k <= order is false or when N ∉ {BInt64, BigInt}\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.advection!-Union{Tuple{AdvectionData{T,Nsp,Nv,Nsum,timeopt}}, Tuple{timeopt}, Tuple{Nsum}, Tuple{Nv}, Tuple{Nsp}, Tuple{T}} where timeopt where Nsum where Nv where Nsp where T","page":"Functions","title":"SemiLagrangian.advection!","text":"advection!(self::AdvectionData)\n\nAdvection function of a multidimensional function f discretized on mesh\n\nArgument\n\nself::AdvectionData : mutable structure of variables data\n\nReturn value\n\ntrue : means that the advection series must continue\nfalse : means that the advection series is ended.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_charge!-Union{Tuple{AdvectionData{T,Nsp,Nv,Nsum,timeopt} where timeopt}, Tuple{Nsum}, Tuple{Nv}, Tuple{Nsp}, Tuple{T}} where Nsum where Nv where Nsp where T","page":"Functions","title":"SemiLagrangian.compute_charge!","text":"compute_charge!( self::AdvectionData)\n\nCompute charge density\n\nρ(x,t) = ∫ f(x,v,t) dv\n\nArgument\n\nself::AdvectionData : mutable structure of variables data.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_charge!-Union{Tuple{Nsum}, Tuple{Nv}, Tuple{Nsp}, Tuple{T}, Tuple{Array{T,Nsp},Tuple{Vararg{UniformMesh{T},Nv}},Array{T,Nsum}}} where Nsum where Nv where Nsp where T","page":"Functions","title":"SemiLagrangian.compute_charge!","text":"compute_charge!( rho, mesh_v, fvx)\n\nCompute charge density\n\nρ(x,t) = ∫ f(x,v,t) dv\n\nArguments\n\nrho::Array{T,Nsp} : result table correctly sized, it is he output\nt_mesh_v::NTuple{Nv,UniformMesh{T}} : velocity meshes\nf::Array{T,Nsum} : input data\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_ee-Tuple{AdvectionData}","page":"Functions","title":"SemiLagrangian.compute_ee","text":"compute_ee(self::AdvectionData)\n\ncompute electric enegie || E(t,.) ||_L2\n\nArgument\n\nself::AdvectionData : veriable advection data structure.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_ee-Union{Tuple{N}, Tuple{T}, Tuple{Tuple{Vararg{UniformMesh{T},N}},Tuple{Vararg{Array{T,N},N}}}} where N where T","page":"Functions","title":"SemiLagrangian.compute_ee","text":"compute_ee(t_mesh_sp, t_elf)\n\ncompute electric enegie || E(t,.) ||_L2\n\nArguments\n\nt_mesh_sp::NTuple{N,UniformMesh{T}} : tuple of space meshes\nt_elf::NTuple{N,Array{T,N}} : tuple of electric field\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_elfield!-Union{Tuple{AdvectionData{T,Nsp,Nv,Nsum,timeopt} where timeopt}, Tuple{Nsum}, Tuple{Nv}, Tuple{Nsp}, Tuple{T}} where Nsum where Nv where Nsp where T","page":"Functions","title":"SemiLagrangian.compute_elfield!","text":"compute_elfield!( self:AdvectionData)\n\ncomputation of electric field     ∇.e = - ρ\n\nArgument\n\nself::AdvectionData : mutable structure of variables data.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_ke-Union{Tuple{AdvectionData{T,Nsp,Nv,Nsum,timeopt} where timeopt}, Tuple{Nsum}, Tuple{Nv}, Tuple{Nsp}, Tuple{T}} where Nsum where Nv where Nsp where T","page":"Functions","title":"SemiLagrangian.compute_ke","text":"compute_ke(self::AdvectionData)\n\nkinetic Energie\n\n∫∫ v^2 f(x,v,t) dv dx\n\nArguments\n\nself::AdvectionData : mutable structure of variables data.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.compute_ke-Union{Tuple{Nsum}, Tuple{Nv}, Tuple{Nsp}, Tuple{T}, Tuple{Tuple{Vararg{UniformMesh{T},Nsp}},Tuple{Vararg{UniformMesh{T},Nv}},Array{T,Nsum}}} where Nsum where Nv where Nsp where T","page":"Functions","title":"SemiLagrangian.compute_ke","text":"compute_ke(t_mesh_sp, t_mesh_v, f)\n\nkinetic Energie\n\n∫∫ v^2 f(x,v,t) dv dx\n\nArguments\n\nt_mesh_sp::NTuple{Nsp, UniformMesh{T}} : tuple of space meshes\nt_mesh_v::NTuple{Nv, UniformMesh{T}} : tuple of velocity meshes\nf::Array{T,Nsum} : function data.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.get_order-Union{Tuple{AbstractInterpolation{T,edge,order}}, Tuple{order}, Tuple{edge}, Tuple{T}} where order where edge where T","page":"Functions","title":"SemiLagrangian.get_order","text":"get_order(_::AbstractInterpolation{T, edge, order}) where{T, edge, order}\n\nReturn the order of interpolation implementation       \n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.getalpha-Tuple{SemiLagrangian.RotationVar,AdvectionData,Any}","page":"Functions","title":"SemiLagrangian.getalpha","text":"getalpha(pv::RotationVar, self::AdvectionData, ind)\n\nImplementation of the interface function that is called before each interpolation in advection\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.interpolate!-Union{Tuple{order}, Tuple{T}, Tuple{AbstractArray{T,1},AbstractArray{T,1},Int64,Array{Array{T,1},1},AbstractInterpolation{T,SemiLagrangian.InsideEdge,order}}, Tuple{AbstractArray{T,1},AbstractArray{T,1},Int64,Array{Array{T,1},1},AbstractInterpolation{T,SemiLagrangian.InsideEdge,order},Any}} where order where T","page":"Functions","title":"SemiLagrangian.interpolate!","text":"interpolate!( \nfp::AbstractVector{T}, fi::AbstractVector{T}, decint::Int, \nallprecal::Vector{Vector{T}}, \ninterp::AbstractInterpolation{T, InsideEdge, order},\ntabmod=gettabmod(length(fi))\n) where {T, order}\n\napply an offset to the function fi interpolate by interp struct, the result is in fp vector, decint and precal are precompute with get_precal method, the TypeEdge is InsideEdge, it is a marginal case\n\nArguments\n\nfp::AbstractVector : output vector\nfi::AbstractVector : input vector\ndecint : offset in units of dx\nallprecal::Vector{Vector{T}} : vector of vector of length order+1 precompute with get_precal(interp, dec) (dec is the offset)\ninterp::AbstractInterpolation{T, InsideEdge, order} : interpolation implementation, note that TypeEdge is CircEdge\ntabmod=gettabmod(length(fi)) : precompute for \"begin at one\" modulo\n\nReturns :\n\nNo return\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.interpolate!-Union{Tuple{order}, Tuple{T}, Tuple{AbstractArray{T,1},AbstractArray{T,1},Int64,Array{T,1},AbstractInterpolation{T,SemiLagrangian.CircEdge,order}}, Tuple{AbstractArray{T,1},AbstractArray{T,1},Int64,Array{T,1},AbstractInterpolation{T,SemiLagrangian.CircEdge,order},Any}} where order where T","page":"Functions","title":"SemiLagrangian.interpolate!","text":"interpolate!( fp::AbstractVector{T}, \n    fi::AbstractVector{T},\n    decint::Int, \n    precal::Vector{T}, \n    interp::AbstractInterpolation{T, CircEdge, order},\n    tabmod=gettabmod(length(fi)) ) where {T, order}\n\napply an offset to the function fi interpolate by interp struct, the result is in fp vector, decint and precal are precompute with get_precal method, the TypeEdge is CircEdge\n\nArguments\n\nfp::AbstractVector : output vector\nfi::AbstractVector : input vector\ndecint : offset in units of dx\nprecal::Vector : vector of length order+1 precompute with get_precal(interp, dec) (dec is the offset)\ninterp::AbstractInterpolation{T, CircEdge, order} : interpolation implementation, note that TypeEdge is CircEdge\ntabmod=gettabmod(length(fi)) : precompute for \"begin at one\" modulo\n\nReturns :\n\nNo return\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.interpolate!-Union{Tuple{order}, Tuple{edge}, Tuple{T}, Tuple{Any,Any,Any,AbstractInterpolation{T,edge,order}}} where order where edge where T","page":"Functions","title":"SemiLagrangian.interpolate!","text":"interpolate!( fp, fi, dec, interp)\n\napply the offset dec to the function fi interpolate by interp struct, the result is in fp Vector\n\nArguments\n\nfp : output vector of length n\nfi : input vector of length n\ndec : offset in units of dx\ninterp::AbstractInterpolation : interpolation implementation\n\nReturns :\n\nNo return\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.nextstate!-Union{Tuple{AdvectionData{T,Nsp,Nv,Nsum,timeopt} where timeopt}, Tuple{Nsum}, Tuple{Nv}, Tuple{Nsp}, Tuple{T}} where Nsum where Nv where Nsp where T","page":"Functions","title":"SemiLagrangian.nextstate!","text":"nextstate!(self::AdvectionData{T, Nsp, Nv, Nsum})\n\nFunction called at the end of advection function to update internal state of AdvectionData structure\n\nArgument\n\nself::AdvectionData{T, Nsp, Nv, Nsum} : object to update\n\nreturn value\n\nret::Bool : true if the series must continue               false at the end of the series.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.points-Tuple{UniformMesh}","page":"Functions","title":"SemiLagrangian.points","text":"points(mesh::UniformMesh)\n\nGet the points of the mesh\n\nArgument\n\nmesh::UniformMesh : the mesh\n\nReturn\n\npoints : the points of the mesh that is the vector of all points of the mesh except the last\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.sizeall-Tuple{Any}","page":"Functions","title":"SemiLagrangian.sizeall","text":"sizeall(adv::Advection)\n\nReturn a tuple of the sizes of each dimensions\n\nArgument\n\nadv::Advection : Advection structure.\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.vec_k_fft-Union{Tuple{UniformMesh{T}}, Tuple{T}} where T","page":"Functions","title":"SemiLagrangian.vec_k_fft","text":"vec_k_fft(mesh::UniformMesh{T}) where{T}\n\nGet the fft coefficients of the mesh\n\nArgument\n\nmesh::UniformMesh{T} : the mesh\n\nReturn\n\nfft coefficients\n\n\n\n\n\n","category":"method"},{"location":"functions/#SemiLagrangian.width-Tuple{UniformMesh}","page":"Functions","title":"SemiLagrangian.width","text":"width(mesh::UniformMesh)\n\nGet the width of the mesh\n\nArgument\n\nmesh::UniformMesh : the mesh\n\nReturn\n\nwidth : the width that is step*length or distance between left and right edges.\n\n\n\n\n\n","category":"method"},{"location":"api/#Mesh","page":"API","title":"Mesh","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"SemiLagrangian.UniformMesh\nBase.step\nBase.length\nSemiLagrangian.points\nSemiLagrangian.width\nSemiLagrangian.vec_k_fft","category":"page"},{"location":"api/#SemiLagrangian.UniformMesh","page":"API","title":"SemiLagrangian.UniformMesh","text":"UniformMesh{T}\nUniformMesh(start::T, stop::T, length::Int) where {T}\n\n1D uniform mesh data.\n\nArguments\n\nstart::T : beginning of the mesh\nstop::T : end of the mesh\nlength::Int : number of cells of the mesh\n\nImplementation\n\nstep::T  : size step\npoints::Vector{T}: Array with node positions\nwidth::T : Distance between left and right edges.\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.step","page":"API","title":"Base.step","text":"Base.step(mesh::UniformMesh)\n\nGet the step of the mesh\n\nArgument\n\nmesh::UniformMesh : the mesh\n\nReturn\n\nstep : the step of the mesh that is the difference between two contiguous points\n\n\n\n\n\n","category":"function"},{"location":"api/#Base.length","page":"API","title":"Base.length","text":"Base.length(mesh::UniformMesh)\n\nGet the length of the mesh\n\nArgument\n\nmesh::UniformMesh : the mesh\n\nReturn\n\nlength : the length of the mesh that is the number of points and cells\n\n\n\n\n\n","category":"function"},{"location":"api/#SemiLagrangian.points","page":"API","title":"SemiLagrangian.points","text":"points(mesh::UniformMesh)\n\nGet the points of the mesh\n\nArgument\n\nmesh::UniformMesh : the mesh\n\nReturn\n\npoints : the points of the mesh that is the vector of all points of the mesh except the last\n\n\n\n\n\n","category":"function"},{"location":"api/#SemiLagrangian.width","page":"API","title":"SemiLagrangian.width","text":"width(mesh::UniformMesh)\n\nGet the width of the mesh\n\nArgument\n\nmesh::UniformMesh : the mesh\n\nReturn\n\nwidth : the width that is step*length or distance between left and right edges.\n\n\n\n\n\n","category":"function"},{"location":"api/#SemiLagrangian.vec_k_fft","page":"API","title":"SemiLagrangian.vec_k_fft","text":"vec_k_fft(mesh::UniformMesh{T}) where{T}\n\nGet the fft coefficients of the mesh\n\nArgument\n\nmesh::UniformMesh{T} : the mesh\n\nReturn\n\nfft coefficients\n\n\n\n\n\n","category":"function"},{"location":"api/#Interpolation","page":"API","title":"Interpolation","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"SemiLagrangian.AbstractInterpolation\nSemiLagrangian.get_order\nSemiLagrangian.sol(_::AbstractInterpolation, b::AbstractVector)\nSemiLagrangian.interpolate!","category":"page"},{"location":"api/#SemiLagrangian.AbstractInterpolation","page":"API","title":"SemiLagrangian.AbstractInterpolation","text":"AbstractInterpolation{T, edge, order}\n\nAbstract supertype for all interpolation type\n\nImplementation constraint\n\ntabfct::Vector : this attribut must be on the implementation, it is a table of function of size order+1\n\n\n\n\n\n","category":"type"},{"location":"api/#SemiLagrangian.get_order","page":"API","title":"SemiLagrangian.get_order","text":"get_order(_::AbstractInterpolation{T, edge, order}) where{T, edge, order}\n\nReturn the order of interpolation implementation       \n\n\n\n\n\n","category":"function"},{"location":"api/#SemiLagrangian.sol-Tuple{AbstractInterpolation,AbstractArray{T,1} where T}","page":"API","title":"SemiLagrangian.sol","text":"sol(_::AbstractInterpolation, line::AbstractVector)\n\nInterface method to transform the treated line, by default this method does nothing\n\nArguments :\n\n_::AbstractInterpolation : interpolation implementation\nline::AbstractVector : line to transform\n\nReturn :\n\nThe transformed line\n\n\n\n\n\n","category":"method"},{"location":"api/#SemiLagrangian.interpolate!","page":"API","title":"SemiLagrangian.interpolate!","text":"interpolate!( fp::AbstractVector{T}, \n    fi::AbstractVector{T},\n    decint::Int, \n    precal::Vector{T}, \n    interp::AbstractInterpolation{T, CircEdge, order},\n    tabmod=gettabmod(length(fi)) ) where {T, order}\n\napply an offset to the function fi interpolate by interp struct, the result is in fp vector, decint and precal are precompute with get_precal method, the TypeEdge is CircEdge\n\nArguments\n\nfp::AbstractVector : output vector\nfi::AbstractVector : input vector\ndecint : offset in units of dx\nprecal::Vector : vector of length order+1 precompute with get_precal(interp, dec) (dec is the offset)\ninterp::AbstractInterpolation{T, CircEdge, order} : interpolation implementation, note that TypeEdge is CircEdge\ntabmod=gettabmod(length(fi)) : precompute for \"begin at one\" modulo\n\nReturns :\n\nNo return\n\n\n\n\n\ninterpolate!( \nfp::AbstractVector{T}, fi::AbstractVector{T}, decint::Int, \nallprecal::Vector{Vector{T}}, \ninterp::AbstractInterpolation{T, InsideEdge, order},\ntabmod=gettabmod(length(fi))\n) where {T, order}\n\napply an offset to the function fi interpolate by interp struct, the result is in fp vector, decint and precal are precompute with get_precal method, the TypeEdge is InsideEdge, it is a marginal case\n\nArguments\n\nfp::AbstractVector : output vector\nfi::AbstractVector : input vector\ndecint : offset in units of dx\nallprecal::Vector{Vector{T}} : vector of vector of length order+1 precompute with get_precal(interp, dec) (dec is the offset)\ninterp::AbstractInterpolation{T, InsideEdge, order} : interpolation implementation, note that TypeEdge is CircEdge\ntabmod=gettabmod(length(fi)) : precompute for \"begin at one\" modulo\n\nReturns :\n\nNo return\n\n\n\n\n\ninterpolate!( fp, fi, dec, interp)\n\napply the offset dec to the function fi interpolate by interp struct, the result is in fp Vector\n\nArguments\n\nfp : output vector of length n\nfi : input vector of length n\ndec : offset in units of dx\ninterp::AbstractInterpolation : interpolation implementation\n\nReturns :\n\nNo return\n\n\n\n\n\n","category":"function"},{"location":"contents/#Contents","page":"Contents","title":"Contents","text":"","category":"section"},{"location":"contents/","page":"Contents","title":"Contents","text":"","category":"page"},{"location":"contents/#Index","page":"Contents","title":"Index","text":"","category":"section"},{"location":"contents/","page":"Contents","title":"Contents","text":"","category":"page"},{"location":"quickstart/#Quickstart","page":"Quickstart","title":"Quickstart","text":"","category":"section"},{"location":"quickstart/#Rotation","page":"Quickstart","title":"Rotation","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"First simulation using semi-lagrangian method to get a rotation","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"\nusing SemiLagrangian\nusing LinearAlgebra\n\nfunction exact!(f, mesh1::UniformMesh{T}, mesh2::UniformMesh{T}, tf::T) where {T}\n    for (i, x) in enumerate(mesh1.points), (j, y) in enumerate(mesh2.points)\n        s, c = sincos(tf)\n        xn, yn = c * x - s * y, s * x + c * y\n        f[i,j] = exp(-13*((xn)^2+(yn+T(6//5))^2))\n    end\nend\n\nmesh_sp = UniformMesh( -5.0, 5.0, 256)\nmesh_v = UniformMesh( -5.0, 5.0, 256)\nnbdt=50\ndt = 2pi/nbdt\ninterp_sp = Lagrange(11)\ninterp_v = Lagrange(11)\n\nadv = Advection((mesh_sp,), (mesh_v,), (interp_sp,), (interp_v,), dt; tab_fct=[tan, sin, tan])\nsz = sizeall(adv)\ntabref = zeros(sz)\nexact!(tabref, mesh_sp, mesh_v, 0.0)\n\npvar = getrotationvar(adv)\n\nadvdata = AdvectionData(adv, tabref, pvar)\n\ndiffmax = 0\ndata = getdata(advdata)\nfor ind=1:nbdt\n    while advection!(advdata) end\n    exact!(tabref, mesh_sp, mesh_v, dt*ind)\n    diff = norm(data .- tabref, Inf)\n    diffmax = max(diffmax, diff)\n    println(\"ind=$ind sz=$sz interp=$interp_sp, $interp_v nbdt=$nbdt diff,diffmax=$diff,$diffmax\")\nend\n","category":"page"},{"location":"quickstart/#Vlasov/poisson","page":"Quickstart","title":"Vlasov/poisson","text":"","category":"section"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"For Vlasov Poisson equation ","category":"page"},{"location":"quickstart/","page":"Quickstart","title":"Quickstart","text":"using SemiLagrangian\n\nmesh_sp = UniformMesh( 0.0, 4pi, 32)\nmesh_v = UniformMesh( -6.0, 6.0, 32)\n\nnbdt=50\ndt = 0.1\nepsilon = 0.5\n\ninterp = Lagrange(7)\n\nadv = Advection( (mesh_sp, mesh_sp), (mesh_v, mesh_v), (interp, interp), (interp, interp), dt)\n\nfct_sp(x)=epsilon * cos(x/2) + 1\nfct_v(v)=exp( - v^2 / 2)/sqrt(2pi)\nlgn_sp = fct_sp.(mesh_sp.points)\nlgn_v = fct_v.(mesh_v.points)\n\ndata = dotprod((lgn_sp, lgn_sp, lgn_v, lgn_v))\n\npvar = getpoissonvar(adv)\n\nadvdata = AdvectionData(adv, data, pvar)\n\nt=0\ncompute_charge!(advdata)\ncompute_elfield!(advdata)\nee = compute_ee(advdata)\nke = compute_ke(advdata)\nprintln(\"$t\\t$ee\\t$ke\\t$(ee+ke)\")\n\nfor ind=1:nbdt\n    while advection!(advdata) end\n    t = Float32(dt*ind)\n    compute_charge!(advdata)\n    compute_elfield!(advdata)\n    ee = compute_ee(advdata)\n    ke = compute_ke(advdata)\n    println(\"$t\\t$ee\\t$ke\\t$(ee+ke)\")\nend\n","category":"page"},{"location":"extdataadv/#ExtDataAdv","page":"ExtDataAdv","title":"ExtDataAdv","text":"","category":"section"},{"location":"extdataadv/","page":"ExtDataAdv","title":"ExtDataAdv","text":"The AbstractExtDataAdv interface is an interface which allows to define the treatments which will make it possible to obtain the values to be applied to the advections.","category":"page"},{"location":"extdataadv/#Functions-that-needs-to-be-implemented","page":"ExtDataAdv","title":"Functions that needs to be implemented","text":"","category":"section"},{"location":"extdataadv/#Methods-to-define","page":"ExtDataAdv","title":"Methods to define","text":"","category":"section"},{"location":"extdataadv/","page":"ExtDataAdv","title":"ExtDataAdv","text":"initcoef!(parext::AbstractExtDataAdv, self::AdvectionData) : this method called at the beginning of each advection to initialize parext data. The self.parext mutable structure is the only data that initcoef! can modify otherwise it leads to unpredictable behaviour.\ngetalpha(parext::AbstractExtDataAdv, self::AdvectionData, ind) : return the alpha number that is used for interpolation.\ngetperm(parext::AbstractExtDataAdv, advd::AdvectionData) : get the permutation of the dimension as a function of the current state, the dimension where advection occurs must be first, the dimensions used to compute alpha must be at the end.","category":"page"},{"location":"#SemiLagrangian.jl-Documentation","page":"Home","title":"SemiLagrangian.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Let us consider an abstract scalar advection equation of the form","category":"page"},{"location":"","page":"Home","title":"Home","text":"$ \\frac{∂f}{∂t}+ a(x, t) ⋅ ∇f = 0. $","category":"page"},{"location":"","page":"Home","title":"Home","text":"The characteristic curves associated to this equation are the solutions of  the ordinary differential equations","category":"page"},{"location":"","page":"Home","title":"Home","text":"$ \\frac{dX}{dt} = a(X(t), t) $","category":"page"},{"location":"","page":"Home","title":"Home","text":"We shall denote by X(t x s) the unique solution of this equation  associated to the initial condition X(s) = x.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The classical semi-Lagrangian method is based on a backtracking of  characteristics. Two steps are needed to update the distribution function  f^n+1 at t^n+1 from its value f^n at time t^n :","category":"page"},{"location":"","page":"Home","title":"Home","text":"For each grid point x_i compute X(t^n x_i t^n+1) the value  of the characteristic at t^n which takes the value x_i at  t^n+1.\nAs the distribution solution of first equation verifies f^n+1(x_i) = f^n(X(t^n x_i t^n+1)) we obtain the desired value of f^n+1(x_i) by computing  f^n(X(t^nx_it^n+1) by interpolation as X(t^n x_i t^n+1)  is in general not a grid point.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Eric Sonnendrücker - Numerical methods for the Vlasov equations","category":"page"},{"location":"types/#Types","page":"Types","title":"Types","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"Modules = [SemiLagrangian]\nOrder   = [:type]","category":"page"},{"location":"types/#SemiLagrangian.Advection","page":"Types","title":"SemiLagrangian.Advection","text":"Advection{T, Nsp, Nv, Nsum, timeopt}\nAdvection(\n    t_mesh_sp::NTuple{Nsp, UniformMesh{T}},\n    t_mesh_v::NTuple{Nv, UniformMesh{T}},\n    t_interp_sp::NTuple{Nsp, AbstractInterpolation{T}},\n    t_interp_v::NTuple{Nv, AbstractInterpolation{T}},\n    dt_base::T;\n    tab_coef=[1//2, 1//1, 1//2],\n    tab_fct=[identity, identity, identity],\n    timeopt::TimeOptimization=NoTimeOpt\n)\n\nImmutable structure that contains constant parameters for multidimensional advection\n\nType parameters\n\nT::DataType : type of data\nNsp : number of space dimensions\nNv : number of velocity dimensions\nNsum : the total number of dimensions (Nsum = Nsp + Nv)\ntimeopt::TimeOptimization : time optimization\n\nArguments\n\nt_mesh_sp::NTuple{Nsp, UniformMesh{T}} : tuple of space meshes (one per space dimension)\nt_mesh_v::NTuple{Nv, UniformMesh{T}} : tuple of velocity meshes (one per velocity dimension)\nt_interp_sp::NTuple{Nsp, AbstractInterpolation{T}} : tuple of space interpolations (one per space dimension)\nt_interp_v::NTuple{Nv, AbstractInterpolation{T}} : tuple of velocity interpolations(one per velocity dimension)\ndt_base::T : time delta for one advection series\n\nKeywords\n\ntab_coef=[1//2, 1//1, 1//2] : coefficient table for one advection series, the   coefficients at odd indexes is for space advection series, the coefficients at even indexes is for velocity advection series\ntab_fct=[identity, identity, identity] : function table for one advection series, with the same indexes than tab_coef\n\nImplementation\n\nsizeall : tuple of the sizes of all dimensions (space before velocity)\nt_mesh_sp : tuple of space meshes\nt_mesh_v : tuple of velocity meshes\nt_interp_sp : tuple of space interpolation types\nt_interp_v : tuple of velocity interpolation types\ndt_base::T : time unit of an advection series\ntab_coef : coefficient table\ntab_fct : function table\nv_square : precompute for ke\nnbsplit : number of slices for split\nmpiid : MPI id\n\nThrows\n\nArgumentError : Nsp must be less or equal to Nv.\n\n\n\n\n\n","category":"type"},{"location":"types/#SemiLagrangian.AdvectionData","page":"Types","title":"SemiLagrangian.AdvectionData","text":"AdvectionData{T,Nsp,Nv,Nsum,timeopt}\nAdvectionData(\nadv::Advection{T,Nsp,Nv,Nsum,timeopt}, \ndata::Array{T,Nsum},\nparext)\n\nMutable structure that contains variable parameters of advection series\n\nType parameters\n\nT::DataType : type of data\nNsp : number of space dimensions\nNv : number of velocity dimensions\nNsum : the total number of dimensions (Nsum = Nsp + Nv)\ntimeopt::TimeOptimization : time optimization\n\nArguments\n\nadv::Advection{T,Nsp,Nv,Nsum} : link to the constant data of this advection\ndata::Array{T,Nsum} : Initial data of this advection\nparext : external data of this advection to compute alpha of each interpolations\n\nImplementation\n\nadv::Advection{T,Nsp,Nv,Nsum,timeopt} : link to the constant data of this advection\nstate_coef::Int : state that is the index of tab_coef, it is from one to lenth(tab_coef)\nstate_dim::Int : the dimension index, from 1 to Nsp in space states, from one to Nv in velocity state\ndata:Array{T,Nsum} : it is the working buffer\nbufdata::Vector{T} : vector of the same size of the working buffer\nt_buf::NTuple{Nsum, Array{T,2}} : tuple of buffer that is used to get the linear data for interpolation, one buffer per thread\ncache_alpha::Union{T,Nothing} : cache for precal, the precal is compute only when the alpha or decint values change\ncache_decint::Int64 : for precal cache\ncache_precal::Vector{T} : for precal cache\nparext::ExtDataAdv : external data of this advection to compute alpha of each interpolations\n\nMethods to define\n\ninitcoef!(parext::AbstractExtDataAdv, self::AdvectionData) : this method called at the beginning of each advection to initialize parext data. The self.parext mutable structure is the only data that initcoef! can modify otherwise it leads to unpredictable behaviour.\ngetalpha(parext::AbstractExtDataAdv, self::AdvectionData, ind) : return the alpha number that is used for interpolation.\ngetperm(parext::AbstractExtDataAdv, advd::AdvectionData) : get the permutation of the dimension as a function of the current state, the dimension where advection occurs must be first, the dimensions used to compute alpha must be at the end.\n\n\n\n\n\n","category":"type"},{"location":"types/#SemiLagrangian.B_Spline","page":"Types","title":"SemiLagrangian.B_Spline","text":"B_Spline{T, edge, order} <: AbstractInterpolation{T, edge, order}\n\nAbstract supertype for all bspline interpolation type\n\n\n\n\n\n","category":"type"},{"location":"types/#SemiLagrangian.Lagrange","page":"Types","title":"SemiLagrangian.Lagrange","text":"Lagrange{T, edge, order, N} <: AbstractInterpolation{T, edge, order}\nLagrange(order, T::DataType=Float64; edge::EdgeType=CircEdge)\n\nType containing Lagrange Polynomials coefficients for Lagrange interpolation\n\nType parameters\n\nT : the type of data that is interpolate\nedge::EdgeType : type of edge traitment\norder::Int: order of lagrange interpolation\n\nImplementation :\n\ntabfct::Vector{Polynomial{T}} : vector of all lagrange polynomial, per example the k-th Lagrange polynomial for the designed order is tabfct[k+1]\n\nArguments :\n\norder::Int : the order of interpolation\n[T::DataType=Float64] : The type values to interpolate \n\nKeywords arguments :\n\nedge::EdgeType=CircEdge : type of edge traitment\n\n\n\n\n\n","category":"type"},{"location":"types/#SemiLagrangian.LuSpline","page":"Types","title":"SemiLagrangian.LuSpline","text":"struct LuSpline{T}\nLuSpline(n, t::Vector{T}; iscirc=true, isLU=true) where{T}\n\nStructure of a LU decomposition of circular banded matrix, a LU decomposition can be stored in a Matrix which is equal to L + U - I. For a circular Banded matrix all non zero coefficients are in the band and in the last columns and lines\n\nImplementation\n\nband::Matrix{T} : matrix of size (kl+ku+1, n-kl)\nku : size of the band upper the diagonal\nkl : size of the band lower the diagonal\niscirc : true if and only if original matrix is circular\nisLU : true if LU decomposition has been perform\nlastcols : only in circular case, Matrix of size (n-ku, kl) that represents teh last columns of matrix\nlastrows : only in circular case, Matrix of size (n, ku) that represents the last rows of matrix\n\nArguments\n\n'n' : size of the matrix\n't::Vector{T}` : vector of all values, the size is order+1, where order is the order of the spline.\n\n\n\n\n\n","category":"type"},{"location":"types/#SemiLagrangian.PoissonConst","page":"Types","title":"SemiLagrangian.PoissonConst","text":"PoissonConst{T, Nsp, Nv}\nPoissonConst(adv::Advection{T, Nsp, Nv, Nsum}; isfftbig=true)\n\nConstant data for the computation of poisson coefficients\n\nArguments\n\nadv::Advection{T, Nsp, Nv, Nsum, timeopt} : Advection constant data\nisfftbig=true: if true compute the fttbig structure\n\nImplementation\n\nadv : Advection constant data\n`v_k' : vector of vector of fourier coefficents for integration for each space dimension\nfctv_k : Array of space dimensions of the inverse of the norm of fourier coefficients\npfftbig : Fourier data for space dimensions\n\n\n\n\n\n","category":"type"},{"location":"types/#SemiLagrangian.PoissonVar","page":"Types","title":"SemiLagrangian.PoissonVar","text":"PoissonVar{T, Nsp, Nv} <: AbstractExtDataAdv\nPoissonVar(pc::PoissonConst{T, Nsp, Nv})\n\nmutable structure of variable data for the poisson computation\n\nArguments\n\npc::PoissonConst{T, Nsp, Nv} : poisson constant data\n\nImplementation\n\npc::PoissonConst{T, Nsp, Nv} : poisson constant data\nrho::Array{T, Nsp} : result of the compute_charge that is the sum along velocity dimensions\nt_elfield::NTuple{Nsp,Array{Complex{T}, Nsp}} : electric fields initialized at each beginning of velocity advection subseries\n\n\n\n\n\n","category":"type"}]
}
