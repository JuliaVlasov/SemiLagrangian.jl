export CubicSpline

struct CubicSpline end

"""

    compute_interpolants( n, f)

Compute interpolation coefficients

This function is a julia version of a Fortran code written by
Edwin Chacon Golcher (Institute of Physics of the Czech Academy of Sciences)
   
"""
function compute_interpolants( n::Int, f::Array{Float64})
        
    @assert (n > 27)

    num_terms = 27
    
    coeffs = zeros(Float64,n+3)
    d      = zeros(Float64,n)

    a   = sqrt((2.0+sqrt(3.0))/6.0)
    r_a = 1.0/a
    b   = sqrt((2.0-sqrt(3.0))/6.0)
    b_a = b/a

    d1 = f[1]
    coeff_tmp = 1.0
    for i in 0:num_terms-1
        coeff_tmp *= (-b_a)
        d1 += coeff_tmp*f[n-1-i]
    end

    d[1] = d1*r_a
    for i in 2:n-1
        d[i] = r_a*(f[i] - b*d[i-1])
    end
        
    d1        = d[n-1]
    coeff_tmp = 1.0
    for i in 1:num_terms
        coeff_tmp *= (-b_a)
        d1 += coeff_tmp*d[i]
    end

    coeffs[n] = d1*r_a
    
    for i in n-2:-1:1
        coeffs[i+1] = r_a*(d[i] - b*coeffs[i+2])
    end

    coeffs[1]   = coeffs[n]
    coeffs[n+1] = coeffs[2]
    coeffs[n+2] = coeffs[3]
    coeffs[n+3] = coeffs[4]
    coeffs
end

export interpolate

"""

    interpolate( coeffs, n, f)

Compute interpolatted value at `x` from coefficients get 
from `compute_interpolants`.

This function is a julia version of a Fortran code written by
   
Edwin Chacon Golcher (Institute of Physics of the Czech Academy of Sciences)
   
"""
function interpolate( coeffs :: Array{Float64,1}, 
                      n1     :: Int,
                      x1min   :: Float64, 
                      x1max   :: Float64, 
                      x      :: Float64 )

      rh        = (n1-1) / (x1max - x1min)
      t0        = (x-x1min)*rh
      cell      = floor(t0)
      delta1        = t0 - cell
      cdelta1       = 1.0 - delta1
      icell     = trunc(Int, cell+1)
      cim1      = coeffs[icell]
      ci        = coeffs[icell+1]
      cip1      = coeffs[icell+2]
      cip2      = coeffs[icell+3]
      t1        = 3.0*ci
      t3        = 3.0*cip1
      t2        = cdelta1*(cdelta1*(cdelta1*(cim1 - t1) + t1) + t1) + ci
      t4        =  delta1*( delta1*( delta1*(cip2 - t3) + t3) + t3) + cip1

      (1.0/6.0) * (t2 + t4)

end

"""
     advection!( f, mesh1, v,  dt)

Semi-lagrangian advection function of 2D distribution function represented 
by array `f`. The advection operates along `axis` (=1 is most efficient) 
with speed `v` during `dt`.

It uses cubic splines interpolation.

"""
function advection!( f::Array{Float64,2}, 
		     mesh::UniformMesh, 
                     v::Vector{Float64},  
		     dt::Float64,
		     interp::CubicSpline,
		     axis::Int64 )
    

    if (axis == 1) 

        @assert ( mesh.length == size(f)[1] )
        @assert ( size(v)[1]  == size(f)[2] )
        lx = mesh.stop - mesh.start
        @simd for j in 1:size(v)[1]
            coeffs = compute_interpolants(mesh.length, f[:,j])        
            @inbounds for i in 1:mesh.length
                x_new = mesh.points[i] - dt * v[j]
                x_new = mesh.start + mod(x_new - mesh.start, lx)
                f[i,j] = interpolate(coeffs, mesh.length, mesh.start, mesh.stop, x_new)
            end
        end

    else

        @assert ( mesh.length == size(f)[2] )
        @assert ( size(v)[1]  == size(f)[1] )
    
        lv = mesh.stop - mesh.start
        @simd for i in 1:size(v)[1]
            coeffs = compute_interpolants(mesh.length, f[i,:])       
            @inbounds for j in 1:mesh.length
                v_new = mesh.points[j] - dt * v[i] 
                v_new = mesh.start + mod(v_new - mesh.start, lv)
                f[i,j] = interpolate(coeffs, mesh.length, mesh.start, mesh.stop, v_new)
            end
        end

    end
    
end
