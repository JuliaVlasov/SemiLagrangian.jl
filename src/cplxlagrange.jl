using Polynomials

function cplxlagrange(tab::Array{Complex{T}}, k::CartesianIndex) where {T}
    p = Polynomial([one(T) + zero(T) * im, 0])
    for ind in CartesianIndices(tab)
        if tab[k] != tab[ind]
            p *= Polynomial([-tab[ind], one(Complex{T})]) / (tab[k] - tab[ind])
        end
    end
    return p
end

function getpoly(tab::Array{Complex{T},2}) where {T}
    p = Polynomial([zero(Complex{T}), zero(Complex{T})])
    for ind in CartesianIndices(tab)
        p += (T(ind.I[1]) + im * T(ind.I[2])) * cplxlagrange(tab, ind)
    end
    return p
end

function getpolyvec(tab::Array{Complex{T},2}) where {T}
    p = Polynomial([zero(Complex{T}), zero(Complex{T})])
    order = size(tab, 1) - 1
    origin = div(order, 2)
    vec = [tab[:, origin]; tab[origin, :]]
    for ind in CartesianIndices(vec)
        v = if ind.I[1] <= order + 1
            T(ind.I[1]) + im * T(origin)
        else
            T(origin) + im * T(ind.I[1])
        end
        p += v * cplxlagrange(vec, ind)
    end
    return p
end
