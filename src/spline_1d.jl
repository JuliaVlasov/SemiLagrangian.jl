using OffsetArrays

struct Spline1D

   degree::Int64
   ncells::Int64
   nbasis::Int64
   start::Float64
   stop::Float64
   step::Float64
   bcoef::Vector{Float64}

   function Spline1D(ncells::Int64, degree::Int64, start, stop)

      nbasis = ncells + degree
      step = (stop - start) / ncells
      bcoef = zeros(Float64, ncells + degree)

      new(degree, ncells, nbasis, start, stop, step, bcoef)

   end

end

function get_cell_and_offset(bspl, x)

   if x == bspl.start
      return 1, 0.0
   elseif x == bspl.stop
      return bspl.ncells, 1.0
   else
      offset = (x - bspl.start) / bspl.step
      icell = min(trunc(Int64, offset), bspl.ncells - 1)
      return icell + 1, min(offset - icell, 1.0)
   end

end

"""
    eval_basis!( spl, x, values )

Evaluate value at x of all basis functions with support in local cell
values[j] = B_j(x) for jmin <= j <= jmin+degree
"""
function eval_basis!(spl, x, values)

   jmin, offset = get_cell_and_offset(spl, x)

   values[1] = 1.0
   for j = 1:spl.degree
      xx = -offset
      saved = 0.0
      for r = 0:j-1
         xx = xx + 1.0
         temp = values[r+1] / j
         values[r+1] = saved + xx * temp
         saved = (j - xx) * temp
      end
      values[j+1] = saved
   end

   jmin

end

function eval_basis(spl, x)

   values = zeros(Float64, spl.degree + 1)

   jmin = eval_basis!(spl, x, values)

   jmin, values

end

"""
    eval_deriv!( derivs, spl, x )

Evaluate derivative at x of all basis functions with support in local cell
derivs[j] = B_j'(x) for jmin <= j <= jmin+degree
"""
function eval_deriv!(derivs, spl, x)

   jmin, offset = get_cell_and_offset(spl, x)

   derivs[1] = 1 / spl.step
   for j = 1:spl.degree-1
      xx = -offset
      saved = 0.0
      for r = 0:j-1
         xx = xx + 1.0
         temp = derivs[r+1] / j
         derivs[r+1] = saved + xx * temp
         saved = (j - xx) * temp
      end
      derivs[j+1] = saved
   end

   bjm1 = derivs[1]
   bj = bjm1
   derivs[1] = -bjm1

   for j = 1:spl.degree-1
      bj = derivs[j+1]
      derivs[j+1] = bjm1 - bj
      bjm1 = bj
   end
   derivs[spl.degree+1] = bj

   jmin

end

function eval_basis_and_n_derivs!(derivs, spl::Spline1D, x::Float64, n::Int64)

   ndu = OffsetArray{Float64}(undef, 0:spl.degree, 0:spl.degree)
   a = OffsetArray{Float64}(undef, 0:1, 0:spl.degree)

   icell, offset = get_cell_and_offset(spl, x)

   jmin = icell

   ndu[0, 0] = 1.0
   for j = 1:spl.degree
      xx = -offset
      saved = 0.0
      for r = 0:j-1
         xx = xx + 1.0
         temp = ndu[r, j-1] / j
         ndu[r, j] = saved + xx * temp
         saved = (j - xx) * temp
      end
      ndu[j, j] = saved
   end

   derivs[0, 1:spl.degree+1] .= ndu[0:spl.degree, spl.degree]

   for r = 0:spl.degree
      s1 = 0
      s2 = 1
      a[0, 0] = 1.0
      for k = 1:n
         d = 0.0
         rk = r - k
         pk = spl.degree - k
         if (r >= k)
            a[s2, 0] = a[s1, 0] / (pk + 1)
            d = a[s2, 0] * ndu[rk, pk]
         end
         if (rk > -1)
            j1 = 1
         else
            j1 = -rk
         end
         if (r - 1 <= pk)
            j2 = k - 1
         else
            j2 = spl.degree - r
         end
         for j = j1:j2
            a[s2, j] = (a[s1, j] - a[s1, j-1]) / (pk + 1)
            d = d + a[s2, j] * ndu[rk+j, pk]
         end
         if (r <= pk)
            a[s2, k] = -a[s1, k-1] / (pk + 1)
            d = d + a[s2, k] * ndu[r, pk]
         end
         derivs[k, r+1] = d
         j = s1
         s1 = s2
         s2 = j
      end
   end

   d = spl.degree / spl.step
   for k = 1:n
      derivs[k, :] .= derivs[k, :] .* d
      d = d * (spl.degree - k) / spl.step
   end

end

"""
    eval_value( spl, x )

Evaluate value of 1D spline at location x: y=S(x)
"""
function eval_value(spl, x)

   values = zeros(Float64, spl.degree + 1)

   jmin = eval_basis!(spl, x, values)

   jmax = jmin + spl.degree

   dot(spl.bcoef[jmin:jmax], values)

end

"""
    eval_deriv( spl, x )

Evaluate derivative of 1D spline at location x: y=S'(x)
"""
function eval_deriv(spl, x)

   derivs = zeros(Float64, spl.degree + 1)

   jmin = eval_deriv!(derivs, spl, x)

   jmax = jmin + spl.degree

   dot(spl.bcoef[jmin:jmax], derivs)

end
