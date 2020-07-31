
function decLU( A::Matrix{T} ) where{T}
    n = size(A,1)
    for k=1:n
        pivot = A[k,k]
        for i=k+1:n
            A[i,k] /= pivot
        end
        for i=k+1:n, j=k+1:n
            # if A[i,k] != 0 && A[k,j] != 0
            #     println("trace avant k=$k i=$i j=$j A[i,j]=$(A[i,j]) A[i,k]=$(A[i,k]) A[k,j]=$(A[k,j])")
            # end
            A[i,j] -= A[i,k]*A[k,j]
            # if A[i,k] != 0 && A[k,j] != 0
            #     println("trace après k=$k i=$i j=$j A[i,j]=$(A[i,j])")
            # end
        end
    end
    return A
end
function decLULu(iscirc, band, lastcols, lastrows)
    wd = size(band, 1)
    ku = div(wd, 2)
    kl = wd-1-ku
    szb = size(band, 2)
    for k=1:szb
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
                band[ku+1+i-j, k+j] -= band[ku+1+i, k] * band[ku+1-j, k+j]
#                println("i=$i j=$j apres band=$(band[ku+1+i-j, k+j])")
            end
        end 
#        println("trace2")

# not finish
    end
 #   lastrows[k-szb, k]
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
    ku = div(size(t,1), 2)
    kl = size(t,1)-1-ku
    for i=1:n
        for (j,v) in enumerate(t)
            ind = i+j-kl-1
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
    isLU::Bool
    lastcols::Union{Matrix{T},Missing} # only when iscirc=true
    lastrows::Union{Matrix{T},Missing} # only when iscirc=true
    function LuSpline(n, t::Vector{T}; iscirc=true, isLU=true) where{T}
        wd = size(t,1) # band width
        ku = div(wd,2)
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
