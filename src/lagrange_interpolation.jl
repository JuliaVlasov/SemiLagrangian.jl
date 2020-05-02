"""
# Lagrange interpolation

- authors: Klaus Reuter (MPCDF) and Katharina Kormann (IPP).

translated from Fortran to Julia by Pierre Navaro (IRMAR)

functions for 1D Lagrange interpolation on a uniform grid (only odd order)

This is an alternative implementation of the Lagrange interpolation for
equidistant grids. The only function implemented is an interpolation for a
given displacement (interpolate_array_disp). The purpose of this implementation
is to provide a fast alternative that exploits the simplifications in this
special case.

Note: The implementation is based on the formulas in Abramowitz and Stegun:
Handbook of Mathematical Functions, Chapter 25.2
"""
struct Lagrange

    stencil :: Int
    pp :: Vector{Float64}

    Lagrange(stencil) = new(stencil, zeros(stencil))

end


# --- compile-time constants to avoid run-time division

const inv_6       = 1/6
const inv_12      = 1/12
const inv_24      = 1/24
const inv_36      = 1/36
const inv_48      = 1/48
const inv_120     = 1/120
const inv_144     = 1/144
const inv_240     = 1/240
const inv_576     = 1/576
const inv_720     = 1/720
const inv_1440    = 1/1440
const inv_5040    = 1/5040
const inv_14400   = 1/14400
const inv_17280   = 1/17280
const inv_30240   = 1/30240
const inv_40320   = 1/40320
const inv_80640   = 1/80640
const inv_362880  = 1/362880
const inv_3628800 = 1/3628800

"""
    lagr_4pt_coeff(p)

Even order
Compute coefficients for Lagrange interpolation for normalized displacement \a p
    - pp(4) : Lagrange interpolations coefficients
    - p     : displacement in units of grid spacing
"""
@inline function lagr_4pt_coeff!(pp, p)

    pp[1] = -p*(p-1)*(p-2)*inv_6
    pp[2] = (p*p-1)*(p-2)*0.5
    pp[3] = -p*(p+1)*(p-2)*0.5
    pp[4] = p*(p*p-1)*inv_6

end

"""
   lagr_4pt(f0, f1, f2, p)

  single point 4-pt-lagrange interpolation
  - lagr_4pt : interpolated value
  - fm1 : known function values at point -1 (relative to where we want to interpolate)
  - f0  : known function values at point 0 (relative to where we want to interpolate)
  - f1  : known function values at point 1 (relative to where we want to interpolate)
  - f2  : known function values at point 2 (relative to where we want to interpolate)
  - p   : displacement in units of grid spacing

"""
function lagr_4pt(lag, fm1, f0, f1, f2, p)

    lagr_4pt_coeff!(lag.pp, p)

    lag.pp[1] * fm1 + lag.pp[2] * f0 + lag.pp[3] * f1 + lag.pp[4] * f2

end 


"""
    lagr_4pt_vec(fi, fp, p, index_shift)

vectorizable 4-pt-lagrange interpolation

   - fi : known function values
   - fp : interpolated function values
   - p  : displacement in units of grid spacing (between 0 and 1)
   - index_shift : index shift due to displacement

"""
function lagr_4pt_vec(lag, fi, fp, p, index_shift)

    lagr_4pt_coeff!(lag.pp, p)

    n = length(fi)

    @simd for i = max(2-index_shift,1):min(n-2-index_shift, n)
      fp[i] = (lag.pp[1] * fi[i-1+index_shift] 
             + lag.pp[2] * fi[i+index_shift]   
             + lag.pp[3] * fi[i+1+index_shift] 
             + lag.pp[4] * fi[i+2+index_shift])
    end

end


"""
    lagr_6pt_coeff(pp, p)

Compute coefficients for Lagrange interpolation for normalized displacement \a p
- pp : Lagrange interpolations coefficients
- p  : displacement in units of grid spacing
"""
@inline function lagr_6pt_coeff!(pp, p)

    pp[1] = - p*(p*p-1)*(p-2)*(p-3)*inv_120
    pp[2] =   p*(p-1)*(p*p-4)*(p-3)*inv_24
    pp[3] = -(p*p-1)*(p*p-4)*(p-3)*inv_12
    pp[4] = p*(p+1)*(p*p-4)*(p-3)*inv_12
    pp[5] = - p*(p*p-1)*(p+2)*(p-3)*inv_24
    pp[6] = p*(p*p-1)*(p*p-4)*inv_120

end

"""
    lagr_6pt(fm2, fm1, f0, f1, f2, f3, p)

single point 6-pt-lagrange interpolation
- lagr_6pt : interpolated value
- fm2 : known function values at point -2 (relative to where we want to interpolate)
- fm1 : known function values at point -1 (relative to where we want to interpolate)
- f0  : known function values at point 0 (relative to where we want to interpolate)
- f1  : known function values at point 1 (relative to where we want to interpolate)
- f2  : known function values at point 2 (relative to where we want to interpolate)
- f3  : known function values at point 3 (relative to where we want to interpolate)
- p   : displacement in units of grid spacing
"""
function lagr_6pt(lag, fm2, fm1, f0, f1, f2, f3, p)

    lagr_6pt_coeff!(lag.pp, p)
    
    lag.pp[1] * fm2 + lag.pp[2] * fm1 + lag.pp[3] * f0  + lag.pp[4] * f1  + lag.pp[5] * f2 + lag.pp[6] * f3

end

"""
   lagr_6pt_vec(fi, fp, p, index_shift)

vectorizable 6-pt-lagrange interpolation

- fi(:) : known function values
- fp(:) : interpolated function values
- p     : displacement in units of grid spacing (between 0 and 1)
- index_shift : index shift due to displacement

"""
function lagr_6pt_vec(lag, fi, fp, p, index_shift)

    lagr_6pt_coeff!(lag.pp, p)
    n = length(fi)
    @simd for i=max(3-index_shift,1):min(n-3-index_shift, n)
        fp[i] = ( lag.pp[1] * fi[i-2+index_shift]
                + lag.pp[2] * fi[i-1+index_shift]
                + lag.pp[3] * fi[i+index_shift]
                + lag.pp[4] * fi[i+1+index_shift]
                + lag.pp[5] * fi[i+2+index_shift] 
                + lag.pp[6] * fi[i+3+index_shift])
    end
end

"""
   lagr_8pt_coeff(pp, p)

Compute coefficients for Lagrange interpolation for normalized displacement \a p
- pp(8) : Lagrange interpolations coefficients
- p     : displacement in units of grid spacing
"""
@inline function lagr_8pt_coeff!(pp, p)

    pp[1] = -p*(p-3)*(p-4)*(p^2-4)*(p^2-1)*inv_5040
    pp[2] = p*(p-2)*(p-4)*(p^2-9)*(p^2-1)*inv_720
    pp[3] = -p*(p-1)*(p-4)*(p^2-9)*(p^2-4)*inv_240
    pp[4] = (p-4)*(p^2-9)*(p^2-4)*(p^2-1)*inv_144
    pp[5] = -(p+1)*p*(p-4)*(p^2-9)*(p^2-4)*inv_144
    pp[6] = (p+2)*p*(p-4)*(p^2-9)*(p^2-1)*inv_240
    pp[7] = -(p+3)*p*(p-4)*(p^2-4)*(p^2-1)*inv_720
    pp[8] = p*(p^2-9)*(p^2-4)*(p^2-1)*inv_5040
end

"""
    lagr_8pt(fm3, fm2, fm1, f0, f1, f2, f3, f4, p)

single point 8-pt-lagrange interpolation
- lagr_8pt : interpolated value
- fm3 : known function values at point -3 (relative to where we want to interpolate)
- fm2 : known function values at point -2 (relative to where we want to interpolate)
- fm1 : known function values at point -1 (relative to where we want to interpolate)
- f0  : known function values at point 0 (relative to where we want to interpolate)
- f1  : known function values at point 1 (relative to where we want to interpolate)
- f2  : known function values at point 2 (relative to where we want to interpolate)
- f3  : known function values at point 3 (relative to where we want to interpolate)
- f4  : known function values at point 4 (relative to where we want to interpolate)
- p   : displacement in units of grid spacing
"""
function lagr_8pt(lag, fm3, fm2, fm1, f0, f1, f2, f3, f4, p)

    lagr_8pt_coeff!(lag.pp, p)

    ( lag.pp[1] * fm3
    + lag.pp[2] * fm2
    + lag.pp[3] * fm1
    + lag.pp[4] * f0
    + lag.pp[5] * f1
    + lag.pp[6] * f2
    + lag.pp[7] * f3
    + lag.pp[8] * f4 )

end 


"""
    lagr_8pt_vec(fi, fp, p, index_shift)

vectorizable 6-pt-lagrange interpolation
- fi : known function values
- fp : interpolated function values
- p  : displacement in units of grid spacing (between 0 and 1)
- index_shift : index shift due to displacement
"""
function lagr_8pt_vec(lag, fi, fp, p, index_shift)

    lagr_8pt_coeff!(lag.pp, p)
    n = length(fi)
    @simd for i=max(4-index_shift,1):min(n-4-index_shift, n)
        fp[i] = ( lag.pp[1] * fi(i-3+index_shift) 
                + lag.pp[2] * fi(i-2+index_shift) 
                + lag.pp[3] * fi(i-1+index_shift) 
                + lag.pp[4] * fi(i+index_shift)   
                + lag.pp[5] * fi(i+1+index_shift) 
                + lag.pp[6] * fi(i+2+index_shift) 
                + lag.pp[7] * fi(i+3+index_shift) 
                + lag.pp[8] * fi(i+4+index_shift))
    end
end

"""
    lagr_3pt_coeff(pp, p)

  Compute coefficients for Lagrange interpolation for normalized displacement \a p
  pp : Lagrange interpolations coefficients
  p  : displacement in units of grid spacing
"""
@inline function lagr_3pt_coeff!(pp, p)

    pp[1] = p*(p-1)*0.5
    pp[2] = 1 - p*p
    pp[3] = p*(p+1)*0.5

end


"""
    lagr_3pt(fm1, f0, f1, p)

single point 3-pt-lagrange interpolation

- lagr_3pt : interpolated value
- fm1 : known function values at point -1 (relative to where we want to interpolate)
- f0  : known function values at point 0 (relative to where we want to interpolate)
- f1  : known function values at point 1 (relative to where we want to interpolate)
- p   : displacement in units of grid spacing
"""
function lagr_3pt(lag, fm1, f0, f1, p)

    lagr_3pt_coeff!(lag.pp, p)

    lag.pp[1] * fm1 + lag.pp[2] * f0 + lag.pp[3] * f1
    
end

"""
    lagr_3pt_vec(fi, fp, p)

vectorizable 3-pt-lagrange interpolation
- fi : known function values
- fp : interpolated function values
- p  : displacement in units of grid spacing
"""
function lagr_3pt_vec(lag, fi, fp, p)

    lagr_3pt_coeff!(lag.pp, p)
    n = length(fi)
    @simd for i=2:n-1
        fp[i] = lag.pp[1] * fi[i-1] + lag.pp[2] * fi[i] + lag.pp[3] * fi[i+1]
    end

end


"""
    lagr_5pt_coeff(pp, p)

Compute coefficients for Lagrange interpolation for normalized displacement \a p

- pp(5) : Lagrange interpolations coefficients
- p     : displacement in units of grid spacing

"""
@inline function lagr_5pt_coeff!(pp, p)

    pp[1] = (p*p-1)*p*(p-2)*inv_24
    pp[2] = -(p-1)*p*(p*p-4)*inv_6
    pp[3] = (p*p-1)*(p*p-4)*0.25
    pp[4] = -(p+1)*p*(p*p-4)*inv_6
    pp[5] = (p*p-1)*p*(p+2)*inv_24

end


"""
    lagr_5pt(fm2, fm1, f0, f1, f2, p)

single point 5-pt-lagrange interpolation

- lagr_5pt : interpolated value
- fm2 : known function values at point -2 (relative to where we want to interpolate)
- fm1 : known function values at point -1 (relative to where we want to interpolate)
- f0  : known function values at point 0 (relative to where we want to interpolate)
- f1  : known function values at point 1 (relative to where we want to interpolate))
- f2  : known function values at point 2 (relative to where we want to interpolate)
- p   : displacement in units of grid spacing
"""
function lagr_5pt(lag, fm2, fm1, f0, f1, f2, p)

    lagr_5pt_coeff!(lag.pp, p)

    lag.pp[1] * fm2 + lag.pp[2] * fm1 + lag.pp[3] * f0 + lag.pp[4] * f1 + lag.pp[5] * f2

end

"""
    lagr_5pt_vec(fi, fp, p)

vectorizable 5-pt-lagrange interpolation

- fi : known function values
- fp : interpolated function values
- p  : displacement in units of grid spacing

"""
function lagr_5pt_vec(lag, fi, fp, p)

    lagr_5pt_coeff!(lag.pp, p)

    n = length(fi)

    @simd for i=3:n-2
        fp[i] = ( lag.pp[1] * fi[i-2]
                + lag.pp[2] * fi[i-1]
                + lag.pp[3] * fi[i]  
                + lag.pp[4] * fi[i+1]
                + lag.pp[5] * fi[i+2] )
    end

end


"""
    lagr_7pt_coeff(pp, p)

Compute coefficients for Lagrange interpolation for normalized displacement \a p
    sll_real64, intent(out) :: pp(7) !< Lagrange interpolations coefficients
    sll_real64, intent(in) :: p      !< displacement in units of grid spacing
"""
@inline function lagr_7pt_coeff!(pp, p)

    pp[1] = p*(p-3)*(p^2-4)*(p^2-1)*inv_720
    pp[2] = -p*(p-2)*(p^2-9)*(p^2-1)*inv_120
    pp[3] = p*(p-1)*(p^2-9)*(p^2-4)*inv_48
    pp[4] = -(p^2-9)*(p^2-4)*(p^2-1)*inv_36
    pp[5] = (p+1)*p*(p^2-9)*(p^2-4)*inv_48
    pp[6] = -(p+2)*p*(p^2-9)*(p^2-1)*inv_120
    pp[7] = (p+3)*p*(p^2-4)*(p^2-1)*inv_720

end


"""
    lagr_7pt(fm3, fm2, fm1, f0, f1, f2, f3, p)

single point 7-pt-lagrange interpolation

- lagr_7pt : interpolated value
- fm3 : known function values at point -3 (relative to where we want to interpolate)
- fm2 : known function values at point -2 (relative to where we want to interpolate)
- fm1 : known function values at point -1 (relative to where we want to interpolate)
- f0  : known function values at point 0 (relative to where we want to interpolate)
- f1  : known function values at point 1 (relative to where we want to interpolate))
- f2  : known function values at point 2 (relative to where we want to interpolate)
- f3  : known function values at point 3 (relative to where we want to interpolate)
- p   : displacement in units of grid spacing

"""
function lagr_7pt(lag, fm3, fm2, fm1, f0, f1, f2, f3, p)

    lagr_7pt_coeff!(lag.pp, p)

    ( lag.pp[1] * fm3 
    + lag.pp[2] * fm2 
    + lag.pp[3] * fm1 
    + lag.pp[4] * f0 
    + lag.pp[5] * f1
    + lag.pp[6] * f2
    + lag.pp[7] * f3 )

end 


"""
    lagr_7pt_vec(fi, fp, p)

vectorizable 7-pt-lagrange interpolation
- fi : known function values
- fp : interpolated function values
- p  : displacement in units of grid spacing

"""
function lagr_7pt_vec(lag, fi, fp, p)

    lagr_7pt_coeff!(lag.pp, p)

    n = length(fi)

    @simd for i=4:n-3
        fp[i] = ( lag.pp[1] * fi[i-3] 
                + lag.pp[2] * fi[i-2] 
                + lag.pp[3] * fi[i-1] 
                + lag.pp[4] * fi[i]   
                + lag.pp[5] * fi[i+1] 
                + lag.pp[6] * fi[i+2] 
                + lag.pp[7] * fi[i+3] )
    end    

end


"""
    lagr_9pt_coeff(pp, p)

Compute coefficients for Lagrange interpolation for normalized displacement \a p

"""
@inline function lagr_9pt_coeff!(pp, p)

    pp[1] =  p*(p-4)*(p^2-9)*(p^2-4)*(p^2-1)*inv_40320
    pp[2] = -p*(p-3)*(p^2-16)*(p^2-4)*(p^2-1)*inv_5040
    pp[3] =  p*(p-2)*(p^2-16)*(p^2-9)*(p^2-1)*inv_1440
    pp[4] = -p*(p-1)*(p^2-16)*(p^2-9)*(p^2-4)*inv_720
    pp[5] =  (p^2-16)*(p^2-9)*(p^2-4)*(p^2-1)*inv_576
    pp[6] = -(p+1)*p*(p^2-16)*(p^2-9)*(p^2-4)*inv_720
    pp[7] = (p+2)*p*(p^2-16)*(p^2-9)*(p^2-1)*inv_1440
    pp[8] = -(p+3)*p*(p^2-16)*(p^2-4)*(p^2-1)*inv_5040
    pp[9] = (p+4)*p*(p^2-9)*(p^2-4)*(p^2-1)*inv_40320

end


"""
    lagr_9pt(fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, p)

single point 9-pt-lagrange interpolation

- lagr_9pt : interpolated value
- fm4 : known function values at point -4 (relative to where we want to interpolate)
- fm3 : known function values at point -3 (relative to where we want to interpolate)
- fm2 : known function values at point -2 (relative to where we want to interpolate)
- fm1 : known function values at point -1 (relative to where we want to interpolate)
- f0  : known function values at point 0 (relative to where we want to interpolate)
- f1  : known function values at point 1 (relative to where we want to interpolate))
- f2  : known function values at point 2 (relative to where we want to interpolate)
- f3  : known function values at point 3 (relative to where we want to interpolate)
- f4  : known function values at point 4 (relative to where we want to interpolate)
- p   : displacement in units of grid spacing
"""
function lagr_9pt(lag, fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, p)

    lagr_9pt_coeff!(lag.pp, p)

    ( lag.pp[1] * fm4 
    + lag.pp[2] * fm3 
    + lag.pp[3] * fm2 
    + lag.pp[4] * fm1 
    + lag.pp[5] * f0  
    + lag.pp[6] * f1  
    + lag.pp[7] * f2  
    + lag.pp[8] * f3  
    + lag.pp[9] * f4 )

end


"""
   lagr_9pt_vec(fi, fp, p)

vectorizable 9-pt-lagrange interpolation
- fi : known function values
- fp : interpolated function values
- p  : displacement in units of grid spacing
"""
function lagr_9pt_vec(lag, fi, fp, p)

    lagr_9pt_coeff!(lag.pp, p)

    n = length(fi)

    @simd for i=5:n-4
        fp[i] = ( lag.pp[1] * fi[i-4]
                + lag.pp[2] * fi[i-3]
                + lag.pp[3] * fi[i-2]
                + lag.pp[4] * fi[i-1]
                + lag.pp[5] * fi[i]
                + lag.pp[6] * fi[i+1]
                + lag.pp[7] * fi[i+2]
                + lag.pp[8] * fi[i+3]
                + lag.pp[9] * fi[i+4] )
    end

end


"""
    lagr_11pt_coeff(pp, p)

Compute coefficients for Lagrange interpolation for normalized displacement \a p

- pp[11] : Lagrange interpolations coefficients
- p      : displacement in units of grid spacing

"""
function lagr_11pt_coeff!(pp, p)

    pp[1]  =  p*(p-5)*(p^2-16)*(p^2-9)*(p^2-4)*(p^2-1)*inv_3628800
    pp[2]  = -p*(p-4)*(p^2-25)*(p^2-9)*(p^2-4)*(p^2-1)*inv_362880
    pp[3]  =  p*(p-3)*(p^2-25)*(p^2-16)*(p^2-4)*(p^2-1)*inv_80640
    pp[4]  = -p*(p-2)*(p^2-25)*(p^2-16)*(p^2-9)*(p^2-1)*inv_30240
    pp[5]  =  p*(p-1)*(p^2-25)*(p^2-16)*(p^2-9)*(p^2-4)*inv_17280
    pp[6]  = -(p^2-25)*(p^2-16)*(p^2-9)*(p^2-4)*(p^2-1)*inv_14400
    pp[7]  =  (p+1)*p*(p^2-25)*(p^2-16)*(p^2-9)*(p^2-4)*inv_17280
    pp[8]  = -(p+2)*p*(p^2-25)*(p^2-16)*(p^2-9)*(p^2-1)*inv_30240
    pp[9]  =  (p+3)*p*(p^2-25)*(p^2-16)*(p^2-4)*(p^2-1)*inv_80640
    pp[10] = -(p+4)*p*(p^2-25)*(p^2-9)*(p^2-4)*(p^2-1)*inv_362880
    pp[11] =  (p+5)*p*(p^2-16)*(p^2-9)*(p^2-4)*(p^2-1)*inv_3628800

end 


"""
    lagr_11pt(fm5, fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, f5, p)

single point 11-pt-lagrange interpolation

- lagr_11pt : interpolated value
- fm5 : known function values at point -5 (relative to where we want to interpolate)
- fm4 : known function values at point -4 (relative to where we want to interpolate)
- fm3 : known function values at point -3 (relative to where we want to interpolate)
- fm2 : known function values at point -2 (relative to where we want to interpolate)
- fm1 : known function values at point -1 (relative to where we want to interpolate)
- f0  : known function values at point 0 (relative to where we want to interpolate)
- f1  : known function values at point 1 (relative to where we want to interpolate))
- f2  : known function values at point 2 (relative to where we want to interpolate)
- f3  : known function values at point 3 (relative to where we want to interpolate)
- f4  : known function values at point 4 (relative to where we want to interpolate)
- f5  : known function values at point 5 (relative to where we want to interpolate)
- p   : displacement in units of grid spacing
"""
function lagr_11pt(lag, fm5, fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, f5, p)

    lagr_11pt_coeff!(lag, pp, p)

    ( lag.pp[1]  * fm5 
    + lag.pp[2]  * fm4 
    + lag.pp[3]  * fm3 
    + lag.pp[4]  * fm2 
    + lag.pp[5]  * fm1 
    + lag.pp[6]  * f0  
    + lag.pp[7]  * f1  
    + lag.pp[8]  * f2  
    + lag.pp[9]  * f3  
    + lag.pp[10] * f4  
    + lag.pp[11] * f5  )

end


"""
    lagr_11pt_vec(fi, fp, p)

vectorizable 11-pt-lagrange interpolation

- fi : known function values
- fp : interpolated function values
- p  : displacement in units of grid spacing
"""
function lagr_11pt_vec(lag, fi, fp, p)

    lagr_11pt_coeff!(lag.pp, p)
    n = length(fi)
    @simd for i=6:n-5
      fp[i] = ( lag.pp[1]  * fi[i-5] 
              + lag.pp[2]  * fi[i-4] 
              + lag.pp[3]  * fi[i-3] 
              + lag.pp[4]  * fi[i-2] 
              + lag.pp[5]  * fi[i-1] 
              + lag.pp[6]  * fi[i]   
              + lag.pp[7]  * fi[i+1] 
              + lag.pp[8]  * fi[i+2] 
              + lag.pp[9]  * fi[i+3] 
              + lag.pp[10] * fi[i+4] 
              + lag.pp[11] * fi[i+5] )
    end

end 

"""
    lagrange_interpolation_1d_disp_fixed_no_bc(fi, fp, p, stencil)

Lagrange interpolation, without boundary conditions.  One sided a the outermost points.

- fi : input array of length n
- fp : output array of length n
- p : offset in units of dx (best interpolation result for p close to zero, about [-1,1], but not a requirement)
- stencil : number of points in fi used for interpolation (possible values 3,5)

"""
function lagrange_interpolation_1d_disp_fixed_no_bc(fi, fp, p, stencil)

    lag = Lagrange(stencil)

    if stencil == 5

        fp[1] = lagr_5pt(lag, fi[1], fi[2], fi[3], fi[4], fi[5], p-2)
        fp[2] = lagr_5pt(lag, fi[1], fi[2], fi[3], fi[4], fi[5], p-1)

        lagr_5pt_vec(lag, fi, fp, p)

        fp[end-1] = lagr_5pt(lag, fi[end-4], fi[end-3], fi[end-2], fi[end-1], fi[end], p+1)
        fp[end]   = lagr_5pt(lag, fi[end-4], fi[end-3], fi[end-2], fi[end-1], fi[end], p+2)

    elseif stencil == 3

        fp[1] = lagr_3pt(lag, fi[1], fi[2], fi[3], p-1)
        lagr_3pt_vec(lag, fi, fp, p)
        fp[end] = lagr_3pt(lag, fi[end-2], fi[end-1], fi[end], p+1)

    else

        @error "Lagrange stencil not implemented (3 or 5)."

    end 

end


"""
    interpolate!(fp, fi, p, Lagrange(stencil))

Lagrange interpolation, periodic boundary conditions

- `fi` : input array of length n
- `fp` : output array of length n
- `p` : offset in units of dx (best interpolation when close to zero, about [-1,1], but not a requirement)
- stencil : number of points in fi used for interpolation (currently possible values 3,5,7)
"""
function interpolate!(fp, fi, p, lag :: Lagrange)

    if lag.stencil == 7
        fp[1] = lagr_7pt(lag, fi[end-2], fi[end-1], fi[end], fi[1], fi[2], fi[3], fi[4], p)
        fp[2] = lagr_7pt(lag, fi[end-1], fi[end], fi[1], fi[2], fi[3], fi[4], fi[5], p)
        fp[3] = lagr_7pt(lag, fi[end], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], p)
        lagr_7pt_vec(lag, fi, fp, p)
        fp[end-2] = lagr_7pt(lag, fi[end-5], fi[end-4], fi[end-3], fi[end-2], fi[end-1], fi[end], fi[1], p)
        fp[end-1] = lagr_7pt(lag, fi[end-4], fi[end-3], fi[end-2], fi[end-1], fi[end], fi[1], fi[2], p)
        fp[end] = lagr_7pt(lag, fi[end-3], fi[end-2], fi[end-1], fi[end], fi[1], fi[2], fi[3], p)
    elseif lag.stencil == 5
        fp[1] = lagr_5pt(lag, fi[end-1], fi[end], fi[1], fi[2], fi[3], p)
        fp[2] = lagr_5pt(lag, fi[end], fi[1], fi[2], fi[3], fi[4], p)
        lagr_5pt_vec(lag, fi, fp, p)
        fp[end-1] = lagr_5pt(lag, fi[end-3], fi[end-2], fi[end-1], fi[end], fi[1], p)
        fp[end] = lagr_5pt(lag, fi[end-2], fi[end-1], fi[end], fi[1], fi[2], p)
    elseif lag.stencil == 3
        fp[1] = lagr_3pt(lag, fi[end], fi[1], fi[2], p)
        lagr_3pt_vec(lag, fi, fp, p)
        fp[end] = lagr_3pt(lag, fi[end-1], fi[end], fi[1], p)
    else
        @error "Lagrange stencil not implemented (3, 5 or 7)"
    end

end

"""
    lagrange_interpolation_1d_disp_fixed_periodic(fi, fp, p, stencil)

Lagrange interpolation, periodic boundary conditions

- fi : input array of length n
- fp : output array of length n
- p  : offset in units of dx (best interpolation result for p close to zero, about [-1,1], but not a requirement)
- stencil : number of points in fi used for interpolation (currently possible values 3,5)
"""
function lagrange_interpolation_1d_disp_fixed_periodic(fi, fp, p, stencil)

    lag = Lagrange(stencil)

    if stencil == 7
        fp[1] = lagr_7pt(lag, fi[end-2], fi[end-1], fi[end], fi[1], fi[2], fi[3], fi[4], p)
        fp[2] = lagr_7pt(lag, fi[end-1], fi[end], fi[1], fi[2], fi[3], fi[4], fi[5], p)
        fp[3] = lagr_7pt(lag, fi[end], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], p)
        lagr_7pt_vec(lag, fi, fp, p)
        fp[end-2] = lagr_7pt(lag, fi[end-5], fi[end-4], fi[end-3], fi[end-2], fi[end-1], fi[end], fi[1], p)
        fp[end-1] = lagr_7pt(lag, fi[end-4], fi[end-3], fi[end-2], fi[end-1], fi[end], fi[1], fi[2], p)
        fp[end] = lagr_7pt(lag, fi[end-3], fi[end-2], fi[end-1], fi[end], fi[1], fi[2], fi[3], p)
    elseif stencil == 5
        fp[1] = lagr_5pt(lag, fi[end-1], fi[end], fi[1], fi[2], fi[3], p)
        fp[2] = lagr_5pt(lag, fi[end], fi[1], fi[2], fi[3], fi[4], p)
        lagr_5pt_vec(lag, fi, fp, p)
        fp[end-1] = lagr_5pt(lag, fi[end-3], fi[end-2], fi[end-1], fi[end], fi[1], p)
        fp[end] = lagr_5pt(lag, fi[end-2], fi[end-1], fi[end], fi[1], fi[2], p)
    elseif stencil == 3
        fp[1] = lagr_3pt(lag, fi[end], fi[1], fi[2], p)
        lagr_3pt_vec(lag, fi, fp, p)
        fp[end] = lagr_3pt(lag, fi[end-1], fi[end], fi[1], p)
    else
        @error "Lagrange stencil not implemented (3, 5 or 7)"
    end
end

"""
    lagrange_interpolation_1d_disp_fixed_periodicl(fi, fp, p, stencil)

Lagrange interpolation, periodic boundary conditions, first value repeated at the end
- fi : input array of length n+1
- fp : output array of length n+1
- p  : offset in units of dx (best interpolation result for p close to zero, about [-1,1], but not a requirement)
- stencil : number of points in fi used for interpolation (currently possible values 3,5)
"""
function lagrange_interpolation_1d_disp_fixed_periodicl(fi, fp, p, stencil)

    lag = Lagrange(stencil)
    n = length(fi)-1

    if stencil == 7

        fp[1] = lagr_7pt(lag, fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], p)
        fp[2] = lagr_7pt(lag, fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], p)
        fp[3] = lagr_7pt(lag, fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], p)
        lagr_7pt_vec(lag, fi, fp, p)
        fp[n-1] = lagr_7pt(lag, fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], p)
        fp[n] = lagr_7pt(lag, fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], p)
        fp[end] = fp[1]

    elseif stencil == 5

        fp[1] = lagr_5pt(lag, fi[n-1], fi[n], fi[1], fi[2], fi[3], p)
        fp[2] = lagr_5pt(lag, fi[n], fi[1], fi[2], fi[3], fi[4], p)
        lagr_5pt_vec(lag, fi, fp, p)
        fp[n] = lagr_5pt(lag, fi[n-2], fi[n-1], fi[n], fi[1], fi[2], p)
        fp[n+1] = fp[1]

    elseif stencil == 3

        fp[1] = lagr_3pt(lag, fi[n], fi[1], fi[2], p)
        lagr_3pt_vec(lag, fi, fp, p)
        fp[end] = fp[1]

    else

        @error "Lagrange stencil not implemented."

    end

end

"""
    lagrange_interpolation_1d_fast_centered_periodicl(fi, fp, p, stencil)

Lagrange interpolation centered around the interval of displacement, periodic boundary condtions, first value repeated at the end

"""
function lagrange_interpolation_1d_disp_centered_periodicl(fi, fp, p, stencil)

    lag = Lagrange(stencil)
    n = length(fi)-1
    # compute interval shift
    ip = trunc(Int, p)
    pq = p - ip

    if stencil == 6
        lagr_6pt_vec(lag, fi, fp, pq, ip)
        for i=1:max(0,2-ip)
            fp[i] = lagr_6pt(lag, 
                             fi[mod1(i-2+ip, n)], 
                             fi[mod1(i-1+ip, n)], 
                             fi[mod1(i  +ip, n)], 
                             fi[mod1(i+1+ip, n)], 
                             fi[mod1(i+2+ip, n)], 
                             fi[mod1(i+3+ip, n)], 
                             pq)
        end
        for i=min(n,n-2-ip):n
           fp[i] = lagr_6pt(lag,
                            fi[mod1(i-2+ip, n)], 
                            fi[mod1(i-1+ip, n)],
                            fi[mod1(i  +ip, n)],
                            fi[mod1(i+1+ip, n)],
                            fi[mod1(i+2+ip, n)],
                            fi[mod1(i+3+ip, n)],
                            pq)
        end
        fp[n+1] = fp[1]

    elseif stencil == 4

        lagr_4pt_vec(lag, fi, fp, pq, ip)

        for i=1:max(0,1-ip)
          fp[i] = lagr_4pt(lag,
                           fi[mod1(i-1+ip,n)],
                           fi[mod1(i  +ip,n)],
                           fi[mod1(i+1+ip,n)],
                           fi[mod1(i+2+ip,n)],
                           pq)
        end
        for i=min(n,n-1-ip):n
           fp[i] = lagr_4pt(lag,
                            fi[mod1(i-1+ip, n)], 
                            fi[mod1(i  +ip, n)], 
                            fi[mod1(i+1+ip, n)], 
                            fi[mod1(i+2+ip, n)],
                            pq)
        end
        fp[n+1] = fp[1]

    else

        @error "Lagrange stencil not implemented"

    end

end
