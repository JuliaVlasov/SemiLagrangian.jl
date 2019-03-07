import LinearAlgebra.LAPACK: gbtrf!, gbtrs!

export BandedMatrix, set_element!, factorize!, solve!

mutable struct BandedMatrix

    n    :: Int64
    kl   :: Int64
    ku   :: Int64
    ipiv :: AbstractVector
    q    :: AbstractMatrix

    function BandedMatrix( n, kl, ku )

        ipiv = zeros(Int64,n)
        q = zeros(Float64, (2*kl+ku+1,n))

        new( n, kl, ku, ipiv, q)

    end

end 

function set_element!( matrix, i, j, a_ij )

    matrix.q[matrix.kl+matrix.ku+1+i-j,j] = a_ij

end

function factorize!( matrix )

    matrix.q, matrix.ipiv = gbtrf!( matrix.kl, matrix.ku, matrix.n, matrix.q )

end

function solve!( matrix, b )

    gbtrs!( 'N', matrix.kl, matrix.ku, matrix.n, matrix.q, matrix.ipiv, b )

end 
