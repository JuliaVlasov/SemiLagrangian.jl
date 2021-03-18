@enum EdgeType CircEdge=1 InsideEdge=2

"""
    AbstractInterpolation{T, edge, order, nd}

Abstract supertype for all interpolation type

# Type parameters
- `T` : type of number on witch interpolation works
- `edge` : type of edge treatment
- `order` : order of interpolation
- `nd` : number of dimensions of the interpolation

# Implementation constraint
- `tabfct::Vector` : this attribut must be on the implementation, it is a table of function of size order+1
"""
abstract type AbstractInterpolation{T, edge, order, nd} end
struct NullInterpoLation{T} <: AbstractInterpolation{T, CircEdge, 1, 0} end

"""
    get_order(_::AbstractInterpolation{T, edge, order}) where{T, edge, order}
Return the order of interpolation implementation       
"""
get_order(_::AbstractInterpolation{T, edge, order}) where{T, edge, order}=order

"""
    get_order(_::AbstractInterpolation{T, edge, order}) where{T, edge, order}
Return the numer of dims of interpolation implementation       
"""
get_ndims(_::AbstractInterpolation{T, edge, order, nd}) where{T, edge, order, nd}=nd

"""
    sol(_::AbstractInterpolation, line::AbstractVector)

Interface method to transform the treated line, by default this method does nothing

# Arguments :
- `_::AbstractInterpolation` : interpolation implementation
- `line::AbstractVector` : line to transform

# Return :
The transformed line

"""
sol(_::AbstractInterpolation, b::AbstractArray)=b

isbspline(_::AbstractInterpolation)=false

Base.show(io::IO, interp::AbstractInterpolation)=print(io, typeof(interp))




@inline function get_precal(interp::AbstractInterpolation{T}, decf::T) where{T}
    return [T(fct(decf)) for fct in interp.tabfct]
end

@inline function get_precal!(v::Vector{T}, bsp::AbstractInterpolation{T}, decf::T) where{T}
    v .= get_precal(bsp,decf)
end

# modulo for "begin to one" array
modone(ind, n)=(n+ind-1)%n+1
gettabmod(lg)=modone.(1:10lg, lg) # 
function get_allprecal(interp::AbstractInterpolation{T, InsideEdge,order}, decint::Int, decfloat::T) where {T,order}
    origin = -div(order, 2)
    indbeg = origin+decint
    indend = indbeg+order
    return [get_precal(interp, decfloat+i) for i=indbeg:indend]
end

"""
    interpolate!( fp::AbstractVector{T}, 
        fi::AbstractVector{T},
        decint::Int, 
        precal::Vector{T}, 
        interp::AbstractInterpolation{T, CircEdge, order},
        tabmod=gettabmod(length(fi)) ) where {T, order}

apply an offset to the function fi interpolate by interp struct, the result is in fp vector,
decint and precal are precompute with get_precal method, the TypeEdge is CircEdge

# Arguments
- `fp::AbstractVector` : output vector
- `fi::AbstractVector` : input vector
- `decint` : offset in units of dx
- `precal::Vector` : vector of length order+1 precompute with get_precal(interp, dec) (dec is the offset)
- `interp::AbstractInterpolation{T, CircEdge, order}` : interpolation implementation, note that TypeEdge is CircEdge
- `tabmod=gettabmod(length(fi))` : precompute for "begin at one" modulo

# Returns :
- No return
"""
function interpolate!( 
    fp::AbstractVector{T}, 
    fi::AbstractVector{T},
    decint::Int, 
    precal::Vector{T}, 
    interp::AbstractInterpolation{T, CircEdge, order, 1},
    tabmod=gettabmod(length(fi))
) where {T, order}
    res = sol(interp,fi)
    origin = -div(order,2)
    lg = length(fi)
    for i=1:lg
        indbeg=i+origin+decint+5lg
        indend=indbeg+order
        fp[i] = sum(res[tabmod[indbeg:indend]] .* precal)
    end
    missing  
end
"""
    interpolate!( 
    fp::AbstractVector{T}, fi::AbstractVector{T}, decint::Int, 
    allprecal::Vector{Vector{T}}, 
    interp::AbstractInterpolation{T, InsideEdge, order},
    tabmod=gettabmod(length(fi))
    ) where {T, order}

apply an offset to the function fi interpolate by interp struct, the result is in fp vector,
decint and precal are precompute with get_precal method, the TypeEdge is InsideEdge, it is a marginal case

# Arguments
- `fp::AbstractVector` : output vector
- `fi::AbstractVector` : input vector
- `decint` : offset in units of dx
- `allprecal::Vector{Vector{T}}` : vector of vector of length order+1 precompute with get_precal(interp, dec) (dec is the offset)
- `interp::AbstractInterpolation{T, InsideEdge, order}` : interpolation implementation, note that TypeEdge is CircEdge
- `tabmod=gettabmod(length(fi))` : precompute for "begin at one" modulo

# Returns :
- No return
"""
function interpolate!( 
    fp::AbstractVector{T}, fi::AbstractVector{T}, decint::Int, 
    allprecal::Vector{Vector{T}}, 
    interp::AbstractInterpolation{T, InsideEdge, order},
    tabmod=gettabmod(length(fi))
) where {T, order}
    res = sol(interp,fi)
    origin = -div(order,2)
    lg = length(fi)
    lgp = length(allprecal)
    borne1=-decint-origin
    borne2=lg-decint+origin-1
    for i=1:borne1
        indbeg=1
        indend=order+1
        ind=i
#        @show 1, i, ind, indbeg, indend, decint, lgp
        fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
    end
    for i=borne1+1:borne2
        indbeg = i-borne1
        indend = indbeg+order
        ind = borne1+1
#        @show 2, i, ind, indbeg, indend, decint, lgp
        fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
    end
    for i=borne2+1:lg
        indbeg = lg-order
        indend = lg
        ind = lgp-(lg-i)
#        @show 3, i, ind, indbeg, indend, decint, lgp
        fp[i] = sum(res[indbeg:indend] .* allprecal[ind])
    end
    missing  
end

"""
    interpolate!( fp, fi, dec, interp)

apply the offset dec to the function fi interpolate by interp struct, the result is in fp Vector

# Arguments
- `fp` : output vector of length n
- `fi` : input vector of length n
- `dec` : offset in units of dx
- `interp::AbstractInterpolation` : interpolation implementation

# Returns :
- No return
"""
function interpolate!( fp, fi, dec, interp::AbstractInterpolation{T, edge, order}) where {T, edge, order}
    decint = convert(Int, floor(dec)) 
    decfloat = dec - decint
    if edge == CircEdge
        return interpolate!(fp, fi, decint, get_precal(interp, decfloat), interp )
    else
        return interpolate!(fp, fi, decint, get_allprecal(interp, decint, decfloat), interp )
    end
    missing
end

# function interpolate!(
#     fp::AbstractArray{T,2}, 
#     fi::AbstractArray{T,2}, 
#     dec::NTuple{2, Array{T, 2}}, 
#     interp::AbstractInterpolation{T, edge, order, 2},
#     tabmod=gettabmod.(size(fi))
# ) where {T, edge, order}

#     res = sol(interp, fi)

#     tabfct = interp.tabfct
#     decint1 = Int.(floor.(dec[1]))
#     decint2 = Int.(floor.(dec[2]))
#     decfl1 = dec[1] - decint1
#     decfl2 = dec[2] - decint2

#     origin = -div(order,2)
#     dec1=origin+5size(fp,1)
#     dec2=origin+5size(fp,2)

#     for i=1:size(fp,1), j=1:size(fp,2)
#         deb_i = i + decint1[i,j]+dec1
#         end_i = deb_i+order
#         deb_j = j+decint2[i,j]+dec2
#         end_j = deb_j + order
# #            @show size(tabdec1), size(tabdec2)
# #    tab = tabdec1[:, tabmod[1][i+decint1[i]+size(fp,1)]] .* transpose(tabdec2[:, tabmod[2][j+decint2[j]+size(fp,2)]])
# # tab = tabdec1[:, tabmod[2][j+decint1[j]+size(fp,2)]] .* transpose(tabdec2[:, tabmod[1][i+decint1[i]+size(fp,1)]])
#         fl_i = decfl1[i,j]
#         fl_j = decfl2[i,j]
# #        tab = [f(fl_i) for f in tabfct] .* transpose([f(fl_j) for f in tabfct])
#         tab = dotprod(([f(fl_i) for f in tabfct], [f(fl_j) for f in tabfct]))
        
# #            @show size(tab), size(fi[tabmod[1][deb_i:end_i], tabmod[2][deb_j:end_j]])
#         fp[i,j] = sum( tab .* fi[tabmod[1][deb_i:end_i], tabmod[2][deb_j:end_j]])
#     end
# end
# OK pour n'importe quel nd
function interpolate!(
    fp::AbstractArray{T,nd}, 
    fi::AbstractArray{T,nd}, 
    dec::NTuple{nd,Any}, 
    interp::AbstractInterpolation{T, edge, order,nd}
) where {T, edge, order, nd}

    sz = size(fp)
    tabmod = gettabmod.(sz)

    res = sol(interp, fi)

 
    tabfct = interp.tabfct

    origin = -div(order,2)
    
    decall = 5 .* sz .+ origin

    decfl = zeros(T, nd)
    oldfl = zeros(T, nd)
    deb_i = zeros(Int, nd)
    end_i = zeros(Int, nd)
    tab = dotprod(ntuple(x -> [f(decfl[x]) for f in tabfct], nd))
    for ind in CartesianIndices(sz)
        for i=1:nd
            dfl = dec[i](ind)
            dint = Int(floor(dfl))
            dfl -= dint
            deb = ind.I[i] + dint+ decall[i]
            end_i[i] = deb + order
            deb_i[i] = deb
            decfl[i] = dfl
        end
        if decfl != oldfl 
            tab = dotprod(ntuple(x -> [f(decfl[x]) for f in tabfct], nd))
            oldfl .= decfl
        end
#        c = ntuple(x -> deb_i[x]:end_i[x], nd)
        c = ntuple(x -> tabmod[x][deb_i[x]:end_i[x]], nd)
        # if nd == 3
        #     @show ind, size(tab), size(fi[c...])
        # end 
        fp[ind] = sum(tab .* fi[c...])
    end
end
