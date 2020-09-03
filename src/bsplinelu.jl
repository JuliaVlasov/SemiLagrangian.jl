
# function decLU( A::Matrix{T} ) where{T}
#     n = size(A,1)
#     for k=1:n
#         pivot = A[k,k]
#         for i=k+1:n
#             A[i,k] /= pivot
#         end
#         for i=k+1:n, j=k+1:n
#             # if A[i,k] != 0 && A[k,j] != 0
#             #     println("trace avant k=$k i=$i j=$j A[i,j]=$(A[i,j]) A[i,k]=$(A[i,k]) A[k,j]=$(A[k,j])")
#             # end
#             A[i,j] -= A[i,k]*A[k,j]
#             # if A[i,k] != 0 && A[k,j] != 0
#             #     println("trace après k=$k i=$i j=$j A[i,j]=$(A[i,j])")
#             # end
#         end
#     end
#     return A
# end
# function decLDLt( A::Matrix{T} ) where{T}
#     n = size(A,1)
#     L = zeros(T,n,n)
#     D = zeros(T,n,n)
    
#     for j=1:n
#         L[j,j] = 1
#         D[j,j] = A[j, j] - (j!=1 ? sum(L[j,k]^2 * D[k,k] for k=1:j-1) : 0)
#         for i=j+1:n
#             L[i,j] = (A[i,j] -(j!=1 ? sum(L[i,k]*L[j,k] * D[k,k] for k=1:j-1) : 0 ))/D[j,j]
#         end
#     end
#     return L, D
# end
# function getLU(A::Matrix{T}) where{T}
#     n = size(A,1)
#     L = zeros(T,n,n)
#     U = zeros(T,n,n)
#     for i=1:n
#         L[i,i] = 1
#         for j=1:i-1
#             L[i,j] = A[i,j]
#         end
#         for j=i:n
#             U[i,j] = A[i,j]
#         end
#     end
#     return L, U
# end

# function topl(n, t, iscirc=true)
#     res=zeros(Rational{BigInt},n,n)
#     kl, ku = get_kl_ku(size(t,1))
#     for i=1:n
#         for (j,v) in enumerate(t)
#             ind = i+j-kl-1
#             if 1 <= ind <= n
#                 res[i, ind] = v
#             elseif iscirc
#                 if ind < 1
#                     ind += n
#                 else
#                     ind -= n
#                 end
#                 res[i,ind] = v
#             end
#         end
#     end
#     return res
# end
# # gradconj : solve A*x = b for A symetric and positive
# function gradconj(A::Matrix{T}, b::Vector{T}, tol) where T
#     n = size(A,1)
#     x = zeros(T,n)
#     r = b - A*x
#     p = copy(r)
#     r2 = sum(r .^ 2)
#     for k=1:1000
#         ap = A*p
#         alpha = r2 / sum( p .* ap )
#         x .+= alpha*p
#         r .-= alpha*ap
#         nm = norm(r, Inf)
#         println("k=$k prec=$nm")
#         if nm < tol
#             break
#         end
#         newr2 = sum(r .^ 2)
#         beta = newr2/r2
#         p .= r + beta*p
#         r2 = newr2
#     end
#     return x
# end
# # gradconj_np : solve A*x = b for A symetric and not necessary positive
# function gradconj_np(A::AbstractMatrix{T}, b::Vector{T}, tolabs, tolrel) where T
#     n = size(A,1)
#     x = zeros(T,n)
#     r = A*b - A*(A*x)
#     p = copy(r)
#     r2 = sum(r .^ 2)
#     for k=1:1000
#         ap = A*(A*p)
#         alpha = r2 / sum( p .* ap )
#         x .+= alpha*p
#         r .-= alpha*ap
#         nm = norm(r)
#         nmrel = nm/norm(x)
#         println("k=$k prec=$nm precrel=$nmrel")
#         if nm < tolabs && nmrel < tolrel
#             break
#         end
#         newr2 = sum(r .^ 2)
#         beta = newr2/r2
#         p .= r + beta*p
#         r2 = newr2
#     end
#     return x
# end




function decLULu(iscirc, band, lastcols, lastrows)
    wd = size(band, 1)
    ku = div(wd, 2)
    kl = wd-1-ku
    szb = size(band, 2)
    n = iscirc ? size(lastrows, 2) : szb
    begrow = n-ku
    begcol = n-kl
    for k=1:(iscirc ? begrow : n)
        pivot = band[ku+1, k]
        for i=ku+2:wd
            band[i,k] /= pivot
        end
        if iscirc
            for i=1:ku
                lastrows[i, k] /= pivot
            end
        end
#        println("trace1 k=$k")
        for i=1:kl, j=1:ku
            if k+j <= szb
#                println("i=$i j=$j avant band=$(band[ku+1+i-j, k+j]) b1=$(band[ku+1+i, k]) b2=$(band[ku+1+j, k+j])")
                # A[k+i,k+j] -= A[k+i, k]*A[k, k+j]
                band[ku+1+i-j, k+j] -= band[ku+1+i, k] * band[ku+1-j, k+j]
#                println("i=$i j=$j apres band=$(band[ku+1+i-j, k+j])")
            end
        end
        if iscirc
            bkl = min(kl, begrow-k)
            for i=1:kl, j=1:bkl
                # A[k+j, begcol+i] = A[k+j, k]*A[k, begcol+i]
                lastcols[k+j, i] -= band[ku+1+j, k]*lastcols[k, i]
            end
            bku = min(ku, szb-k)
            for i=1:bku, j=1:ku
                # A[begrow+j, k+i] = A[begrow+j, k]*A[k, k+i]
                lastrows[j, k+i] -= lastrows[j, k]* band[ku+1-i, k+i]
            end
            for i=1:ku, j=1:kl
                # A[begrow+i, begcol+j] = A[begrow+i, k]*A[k, begcol+j]
                lastrows[i,begcol+j] -= lastrows[i, k]*lastcols[k,j]
            end
        end
    end
    if iscirc
        for k=begrow+1:n
            i_k = k-begrow
            pivot = lastrows[i_k, k]
            for i=i_k+1:ku
                lastrows[i,k] /= pivot
            end
            for i=i_k+1:ku, j=k+1:n
                # A[i+begrow, j] -= A[i+begrow, k]*A[k, j]
                lastrows[i, j] -= lastrows[i, k]*lastrows[k-begrow, j]
            end
        end
    end
end


   

function get_kl_ku(order)
    ku = div(order,2)
    kl = order-1-ku
    return kl, ku
end


struct LuSpline{T}
    band::Matrix{T}
    ku::Int64
    kl::Int64
    iscirc::Bool
    isLU::Bool
    lastcols::Union{Matrix{T},Missing} # missing when iscirc=false
    lastrows::Union{Matrix{T},Missing} # missing when iscirc=false
    function LuSpline(n, t::Vector{T}; iscirc=true, isLU=true) where{T}
        wd = size(t,1) # band width
        kl, ku = get_kl_ku(wd)
        szb = iscirc ? n-kl : n
        band = zeros(T, wd, szb)
        for i=1:wd
            jbeg = i <= ku + 1 ? ku - i + 2 : 1
            if iscirc
                jend = i >= kl + 2 ? szb - i + kl + 1 : szb
            else
                jend = i >= ku + 2 ? n - i + ku + 1 : n
            end
            for j=jbeg:jend
                band[i,j] = t[i]
            end
        end
        if iscirc
            lastrows = zeros(T, ku, n)
            lastcols = zeros(T, n-ku, kl)
            for i=1:ku
                for ind=1:wd
                    j=n-wd+i+ind
                    lastrows[i,(j-1)%n+1] = t[ind]
                end
            end
            for i=1:kl
                for j=1:kl+1-i
                    lastcols[j,j+i-1] = t[i]
                    lastcols[n-kl-ku+i+j-1,j] = t[i]
                 end
            end
        else
            lastcols = lastrows = missing
        end
        if isLU
            decLULu(iscirc, band, lastcols, lastrows)
        end
        return new{T}(band, ku, kl, iscirc, isLU, lastcols, lastrows)
    end
    function LuSpline( A::Matrix{T}, ku, kl; iscirc=true, isLU=false) where {T}
        n = size(A,1)
        if iscirc
            lastrows = copy(A[end-ku+1:end,:])
            lastcols = copy(A[1:end-ku,end-kl+1:end])
        else
            lastcols = lastrows = missing
        end
        szb = iscirc ? n-kl : n
        wd = kl+ku+1
        band = zeros(T,wd,szb)
        for j=1:szb
            for i=j-ku:j+kl
                if 1 <= i <= n - (iscirc ? ku : 0)
                    band[ku+i+1-j,j] = A[i,j]
                end
            end
        end
        return new{T}(band, ku, kl, iscirc, isLU, lastcols, lastrows)
    end
end
function ==(la::LuSpline{T}, lb::LuSpline{T}) where{T}
    return (la.ku == lb.ku && la.kl == la.kl && la.iscirc == lb.iscirc 
            && la.isLU == lb.isLU && la.band == lb.band 
            && (!la.iscirc || (la.lastrows == lb.lastrows && la.lastcols == lb.lastcols)))
end
function sol( A::Matrix{T}, Y::Vector{T}) where{T}
    L, U = getLU(A)
    n = size(A,1)
    Y1 = zeros(T,n)
    for i=1:n
        Y1[i] = Y[i] - sum(Y1[1:i-1] .* A[i, 1:i-1])
    end
 #   @assert Y1 == (L^(-1))*Y "Erreur 1" 
 #   @assert isapprox(Y1, (L^(-1))*Y, atol=1e-60) "Erreur 1" 
    X = zeros(T,n)
    for i=n:-1:1
        X[i] = (Y1[i] - sum(X[i+1:n] .* A[i, i+1:n]))/A[i,i]
    end
#    @assert X == (U^(-1))*Y1 "Erreur 2" 
#    @assert isapprox(X,(U^(-1))*Y1,atol=1e-60) "Erreur 2" 
#    @assert X == ((L*U)^(-1))*Y "Erreur3"
#    @assert isapprox(X, ((L*U)^(-1))*Y, atol=1e-60) "Erreur3"
    return X, Y1
end

function sol(spA::LuSpline{T}, b::Vector{T}) where{T}
    szb = size(spA.band,2)
    n = spA.iscirc ? size(spA.lastrows, 2) : szb
    begrow = n-spA.ku
    begcol = n-spA.kl
    Y = zeros(T, n)
    Y = copy(b)
    endmat = spA.iscirc ? begrow : n
    endmat2 = spA.iscirc ? begcol : n
    for i=2:endmat
        fin = i-1
        deb = max( 1, i-spA.kl)
        Y[i] -= sum([ Y[j]*spA.band[spA.ku+1+i-j, j] for j=deb:fin])
    end
    if spA.iscirc
        for i=begrow+1:n
            Y[i] -= sum( Y[1:i-1] .* spA.lastrows[i-begrow,1:i-1])
        end
    end
    X = zeros(T,n)
    if spA.iscirc
        for i=n:-1:begrow+1
            X[i] =(Y[i] - sum(X[i+1:n] .* spA.lastrows[i-begrow,i+1:n]))/spA.lastrows[i-begrow,i]
        end
    end
    for i=endmat:-1:1
        deb = i+1
        fin = min(i+spA.ku, endmat2)
        if deb <= fin
            s = sum(spA.band[ spA.ku+1+i-j,j]   * X[j] for j=deb:fin)
        else
            s = 0
        end
        if spA.iscirc
            s += sum(spA.lastcols[i,1:spA.kl] .* X[end-spA.kl+1:end])
        end
        X[i] = (Y[i]-s)/spA.band[spA.ku+1,i]
    end
    return X, Y
end
get_n(sp::LuSpline)=sp.iscirc ? size(sp.lastrows, 2) : size(sp.band, 2)
get_order(sp::LuSpline)=sp.ku+sp.kl+1




abstract type InterpolationType end
abstract type B_Spline{T,iscirc} <: InterpolationType end
struct B_SplineLU{T,iscirc} <: B_Spline{T,iscirc}
    ls::LuSpline{T}
    bspline::Spline
    function B_SplineLU( order, n, eltfortype::T; iscirc=true) where{T}
        bspline = getbspline(order, 0)
        ls = LuSpline(n,convert.(T,bspline.(1:order)), iscirc=iscirc, isLU=true)
        return new{T, iscirc}(ls, bspline)
    end
    # B_Spline(o, n, t::DataType ; kwargs... )=B_Spline(o, n, one(t) ; kwargs... )
end

sol(bsp::B_SplineLU{T}, b::Vector{T}) where{T}=sol(bsp.ls, b)[1]
get_n(bsp::B_SplineLU{T}) where{T}=get_n(bsp.ls)
get_order(bsp::B_SplineLU{T}) where{T}=get_order(bsp.ls)
get_bspline(bsp::B_SplineLU{T}) where{T}=bsp.bspline
get_type(bsp::B_SplineLU{T, iscirc}) where{T,iscirc}="B_SplineLU{$T, $iscirc}"

function interpolate!( adv, fp, fi, dec, 
    bsp::B_Spline{T, iscirc}
) where {T, iscirc}
    println("begin interpolate with $(get_type(bsp))")
    res = sol(bsp, fi)
    n = get_n(bsp)
    order = get_order(bsp)
    kl, ku = get_kl_ku(order)

    decint = convert(Int, floor(dec))
    decfloat = dec-decint
    if decfloat > 0.5
        decfloat -= 1.
        decint += 1
    end

#    println("n=$n size(res)=$(size(res,1)) size(fi)=$(size(fi,1)) size(fp)=$(size(fp,1))")
    precal = get_bspline(bsp).((1:order) .- decfloat)
    if iscirc
        decall = n-kl+decint-2
        for i=1:n
            fp[i] = 0
            dec = decall + i
            for j=1:order
                # fp[i] += res[(i+n-bsp.ls.kl-2+j)%n+1]*bsp.bspline(j-1+dec)
                fp[i] += res[(dec+j)%n+1]*precal[j]
            end
            # diff = fp[i] -fi[i]
            # println("i2=$i diff=$diff")
        end
    else
        decall = -kl+decint-1
        for i=1:n
            # deb = max(1, bsp.ls.kl+2-i)
            # fin = min(order, n-i + bsp.ls.ku+1)
            fp[i] = 0
            dec = decall+i
            for j=1:order
                ind = dec+j
                if 1 <= ind <= n
                    fp[i] += res[ind]*precal[j]
                end
            end
        end
    end
    println("end interpolate with $(get_type(bsp))")
end
