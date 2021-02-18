import Base: ==

function decLULu(iscirc, band, lastcols, lastrows)
    wd = size(band, 1)
    kl, ku = get_kl_ku(wd)
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
                lastrows[i, begcol+j] -= lastrows[i, k]*lastcols[k, j]
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

"""
    struct LuSpline{T}
    LuSpline(n, t::Vector{T}; iscirc=true, isLU=true) where{T}

Structure of a LU decomposition of circular banded matrix,
a LU decomposition can be stored in a Matrix which is equal to L + U - I.
For a circular Banded matrix all non zero coefficients are in the band and in the last columns and lines


# Implementation
- band::Matrix{T} : matrix of size (kl+ku+1, n-kl)
- ku : size of the band upper the diagonal
- kl : size of the band lower the diagonal
- iscirc : true if and only if original matrix is circular
- isLU : true if LU decomposition has been perform
- lastcols : only in circular case, Matrix of size (n-ku, kl) that represents teh last columns of matrix
- lastrows : only in circular case, Matrix of size (n, ku) that represents the last rows of matrix

# Arguments
- 'n' : size of the matrix
- 't::Vector{T}` : vector of all values, the size is order+1, where order is the order of the spline.

"""
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
"""
    LuSpline( A::Matrix{T}, ku, kl; iscirc=true, isLU=false) where {T}

Contructor from a matrix, for test only

"""
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
function sol(spA::LuSpline{T}, b::AbstractVector{T}) where{T}
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

"""
    B_SplineLU{T,edge, order, N} <: B_Spline{T, edge, order}

Type containing spline coefficients for b-spline interpolation

# Type parameters
- `T` : the type of data that is interpolate
- `edge::TypeEdge=CircEdge` : true if function is circular
- `order::Int`: order of lagrange interpolation
- `N` : type of integer, in fact Int64 or BigInt that is used in the Spline object

# Implementation :
- `ls::LuSpline{T}` : the LU matrix
- `bspline::SplineInt{N}` : object that contain the bspline

# Arguments : 
- `n` : size of the matrix
- `order` : the order of interpolation
- `[T::DataType=Float64]` : The type values to interpolate 

"""

struct B_SplineLU{T, edge, order} <: B_Spline{T, edge, order}
    ls::LuSpline{T}
    tabfct::Vector{Polynomial{T}}
    function B_SplineLU( order::Int, n::Int, T::DataType=Float64)
        (order%2 == 0) && throw(ArgumentError("order=$order B_SplineLU for even  order is not implemented n=$n")) 
        bspline = getbspline(order, 0)
        tabfct_rat = map(x -> bspline[order-x](Polynomial([order-x,1])), 0:order)
        # bspline = SplineInt(order)
        # N = typeof(bspline.fact_order)
        # tabpol = map(x -> bspline[order-x](Polynomial([order-x,1])), 0:order)
        ls = LuSpline(n,convert.(T, bspline.(1:order)), iscirc=true, isLU=true)
        return new{T, CircEdge, order}(ls, convert.(Polynomial{T}, tabfct_rat))
    end
    B_SplineLU(o::Int, n::Int, elt::T; kwargs...) where {T<:Number}=B_SplineLU(o, n, T; kwargs...)
end


sol(bsp::B_SplineLU{T}, b::AbstractVector{T}) where {T<:Number}=sol(bsp.ls, b)[1]
get_n(bsp::B_SplineLU{T}) where{T}=get_n(bsp.ls)


