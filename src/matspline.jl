
function decLU( A::Matrix{T} ) where{T}
    n = size(A,1)
    for k=1:n
        pivot = A[k,k]
        for i=k+1:n
            A[i,k] /= pivot
        end
        for i=k+1:n, j=k+1:n
            A[i,j] -= A[i,k]*A[k,j]
        end
    end
    return A
end

function getLU(A::Matrix{T}) where{T}
    n = size(A,1)
    L = zeros(T,n,n)
    U = zeros(T,n,n)
    for i=1:n
        L[i,i] = 1
        for j=1:i-1
            L[i,j] = A[i,j]
        end
        for j=i:n
            U[i,j] = A[i,j]
        end
    end
    return L, U
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
    return X
end


function topl(n, t, iscirc=true)
    res=zeros(Rational{BigInt},n,n)
    dec = div(size(t,1), 2)
    for i=1:n
        for (j,v) in enumerate(t)
            ind = i+j-dec-1
            if 1 <= ind <= n
                res[i, ind] = v
            elseif iscirc
                if ind < 1
                    ind += n
                else
                    ind -= n
                end
                res[i,ind] = v
            end
        end
    end
    return res
end

   



struct LuSpline{T}
    band::Matrix{T}
    ku::Int64
    kl::Int64
    iscirc::Bool
    lastcols::Matrix{T} # only when iscirc=true
    lastrows::Matrix{T} # only when iscirc=true
    function LuSpline(n, t::Vector{T}, iscirc=true) where{T}
        wd = size(t,1) # band width
        ku = div(wd-1,2)
        kl = wd-1-ku
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
        return new{T}(band, ku, kl, iscirc, lastcols, lastrows)
    #   decLULu(iscirc, band, lastcols, lastrows)
    end
    function LuSpline( A::Matrix{T}, ku, kl, iscirc) where {T}
        n = size(A,1)
        if iscirc
            lastrows = copy(A[end-ku+1:end,:])
            lastcols = copy(A[1:end-ku,end-kl+1:end])
        end
        szb = iscirc ? n-kl : n
        wd = kl+ku+1
        band = zeros(T,wd,szb)
        for j=1:szb
            for i=j-ku:j+ku
                if 1 <= i <= n-kl
                    band[ku+i+1-j,j] = A[i,j]
                end
            end
        end
        return new{T}(band, ku, kl, iscirc, lastcols, lastrows)
    end
end
function ==(la::LuSpline{T}, lb::LuSpline{T}) where{T}
    return (la.ku == lb.ku && la.kl == la.kl && la.iscirc == lb.iscirc 
            && la.band == lb.band && la.lastrows == lb.lastrows 
            && la.lastcols == lb.lastcols)
end
