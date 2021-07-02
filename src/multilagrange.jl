

using DataStructures

function simplelagrange(v::Vector{T}, indice::Int,dim::Int,  N::Int) where{T}
    p =  MPoly(OrderedDict([0,0] =>T(1)), [:x,:y])
    base = [[1,0],[0,1]][indice]
    for i=1:length(v)
        if i != indice 
            p *= MPoly(OrderedDict(base => T(1), [0,0] => -v[i]), [:x,:y])/(v[indice]-v[i])
        end
    end
    return p
end



function multilagrange(tab::Array{T,2}) where{T}
    N = 2
    sz = size(tab)
    var = [:x,:y]
    order = sz[1]-1
    origin = div(order,2)
    p =  MPoly(OrderedDict([0,0] =>T(0)), [:x,:y])
    tabpol = Array{MPoly{T},2}(undef,sz...)
    for ind in CartesianIndices(sz)
        p1 = simplelagrange(tab[:,ind.I[2]], ind.I[1], 1, N)
        p2 = simplelagrange(tab[ind.I[1], :], ind.I[2], 2, N)
        tabpol[ind] = p1*p2
    end



    




