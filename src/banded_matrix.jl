import LinearAlgebra.LAPACK: gbtrf!, gbtrs!

export BandedMatrix, set_element!, factorize!, solve!

mutable struct BandedMatrix{T}

    n::Int64
    kl::Int64
    ku::Int64
    ipiv::AbstractVector
    q::AbstractMatrix

    function BandedMatrix(x::T, n, kl, ku) where {T<:Union{AbstractFloat,Complex{AbstractFloat}}}

        ipiv = zeros(Int64, n)
        q = zeros(T, (2 * kl + ku + 1, n))

        new{T}(n, kl, ku, ipiv, q)

    end
    BandedMatrix(T1::DataType, n, kl, ku)=BandedMatrix(zero(T1), n, kl, ku)

    BandedMatrix(n, kl, ku)=BandedMatrix(Float64, n, kl, ku)

end

function set_element!(matrix, i, j, a_ij)

    matrix.q[matrix.kl+matrix.ku+1+i-j, j] = a_ij

end

function factorize!(matrix)

 #   println(" avant gbtrf matrix=$matrix")

    matrix.q, matrix.ipiv = gbtrfgen!(matrix.kl, matrix.ku, matrix.n, matrix.q)

#    println(" après gbtrf matrix=$matrix")

end

function solve!(matrix, b)

    gbtrsgen!( matrix.kl, matrix.ku, matrix.n, matrix.q, matrix.ipiv, b)

end
