



struct GeoConst{T,N}
    adv::Advection{T,N}
    coefrsqk::Array{T,N}
    function GeoConst(adv::Advection{T,N}) where {T,N}
        sz = sizeall(adv)
        coefrsqk = zeros(T,sz)
        for ind in CartesianIndices(sz)
            k = sqrt(sum([adv.t_mesh[i].points[ind.I[i]]^2 for i=1:N]))
            coeffrsqk[ind] = k2 == 0 ? 0 : 1/k
        end
        return new{T,N}(adv,coefrsqk)
    end
end

function caldec(geoc::GeoConst{T,N}, advd::AdvectionData{T,N})
    
