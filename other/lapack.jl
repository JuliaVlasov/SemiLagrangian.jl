
function pttrfgen!(D::Vector{T}, E::Vector{T}) where {T<:Union{AbstractFloat,Complex{AbstractFloat},Rational,Complex{Rational}}}
    n1 = size(D,1)
    n2 = size(E,1)
    @assert n1 == n2+1 "n1=$n1 and n2=$n2 n1 must be greater than n2 by one"
    for i=1:n2
        E[i] = E[i]/D[i]
        D[i+1] -= E[i]^2*D[i]
    end
    return (D, E)
end
function pttrsgen!(D::Vector{T}, E::Vector{T}, B::Array{T}) where {T<:Union{AbstractFloat,Complex{AbstractFloat},Rational,Complex{Rational}}}
    n1 = size(D,1)
    n2 = size(E,1)
    n3 = size(B,1)
    @assert n1 == n2+1 "n1=$n1 and n2=$n2, n1 must be greater than n2 by one"
    @assert n1 == n3 "n1=$n1 and n3=$n3, n1 must be equal to n3"
    for i=1:n2
        B[i+1,:] .-= E[i]*B[i,:]
    end
    B ./= D
    for i=n2:-1:1
        B[i,:] -= E[i]*B[i+1,:]
    end
    return B
end

# *> \brief \b DGBTF2 computes the LU factorization of a general band matrix using the unblocked version of the algorithm.
#  *
#  *  =========== DOCUMENTATION ===========
#  *
#  * Online html documentation available at
#  *            http://www.netlib.org/lapack/explore-html/
#  *
#  *> \htmlonly
#  *> Download DGBTF2 + dependencies
#  *> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbtf2.f">
#  *> [TGZ]</a>
#  *> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbtf2.f">
#  *> [ZIP]</a>
#  *> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbtf2.f">
#  *> [TXT]</a>
#  *> \endhtmlonly
#  *
#  *  Definition:
#  *  ===========
#  *
#  *       SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
#  *
#  *       .. Scalar Arguments ..
#  *       INTEGER            INFO, KL, KU, LDAB, M, N
#  *       ..
#  *       .. Array Arguments ..
#  *       INTEGER            IPIV( * )
#  *       DOUBLE PRECISION   AB( LDAB, * )
#  *       ..
#  *
#  *
#  *> \par Purpose:
#  *  =============
#  *>
#  *> \verbatim
#  *>
#  *> DGBTF2 computes an LU factorization of a real m-by-n band matrix A
#  *> using partial pivoting with row interchanges.
#  *>
#  *> This is the unblocked version of the algorithm, calling Level 2 BLAS.
#  *> \endverbatim
#  *
#  *  Arguments:
#  *  ==========
#  *
#  *> \param[in] M
#  *> \verbatim
#  *>          M is INTEGER
#  *>          The number of rows of the matrix A.  M >= 0.
#  *> \endverbatim
#  *>
#  *> \param[in] N
#  *> \verbatim
#  *>          N is INTEGER
#  *>          The number of columns of the matrix A.  N >= 0.
#  *> \endverbatim
#  *>
#  *> \param[in] KL
#  *> \verbatim
#  *>          KL is INTEGER
#  *>          The number of subdiagonals within the band of A.  KL >= 0.
#  *> \endverbatim
#  *>
#  *> \param[in] KU
#  *> \verbatim
#  *>          KU is INTEGER
#  *>          The number of superdiagonals within the band of A.  KU >= 0.
#  *> \endverbatim
#  *>
#  *> \param[in,out] AB
#  *> \verbatim
#  *>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
#  *>          On entry, the matrix A in band storage, in rows KL+1 to
#  *>          2*KL+KU+1; rows 1 to KL of the array need not be set.
#  *>          The j-th column of A is stored in the j-th column of the
#  *>          array AB as follows:
#  *>          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
#  *>
#  *>          On exit, details of the factorization: U is stored as an
#  *>          upper triangular band matrix with KL+KU superdiagonals in
#  *>          rows 1 to KL+KU+1, and the multipliers used during the
#  *>          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
#  *>          See below for further details.
#  *> \endverbatim
#  *>
#  *> \param[in] LDAB
#  *> \verbatim
#  *>          LDAB is INTEGER
#  *>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
#  *> \endverbatim
#  *>
#  *> \param[out] IPIV
#  *> \verbatim
#  *>          IPIV is INTEGER array, dimension (min(M,N))
#  *>          The pivot indices; for 1 <= i <= min(M,N), row i of the
#  *>          matrix was interchanged with row IPIV(i).
#  *> \endverbatim
#  *>
#  *> \param[out] INFO
#  *> \verbatim
#  *>          INFO is INTEGER
#  *>          = 0: successful exit
#  *>          < 0: if INFO = -i, the i-th argument had an illegal value
#  *>          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
#  *>               has been completed, but the factor U is exactly
#  *>               singular, and division by zero will occur if it is used
#  *>               to solve a system of equations.
#  *> \endverbatim
#  *
#  *  Authors:
#  *  ========
#  *
#  *> \author Univ. of Tennessee
#  *> \author Univ. of California Berkeley
#  *> \author Univ. of Colorado Denver
#  *> \author NAG Ltd.
#  *
#  *> \date December 2016
#  *
#  *> \ingroup doubleGBcomputational
#  *
#  *> \par Further Details:
#  *  =====================
#  *>
#  *> \verbatim
#  *>
#  *>  The band storage scheme is illustrated by the following example, when
#  *>  M = N = 6, KL = 2, KU = 1:
#  *>
#  *>  On entry:                       On exit:
#  *>
#  *>      *    *    *    +    +    +       *    *    *   u14  u25  u36
#  *>      *    *    +    +    +    +       *    *   u13  u24  u35  u46
#  *>      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
#  *>     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
#  *>     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
#  *>     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
#  *>
#  *>  Array elements marked * are not used by the routine; elements marked
#  *>  + need not be set on entry, but are required by the routine to store
#  *>  elements of U, because of fill-in resulting from the row
#  *>  interchanges.
#  *> \endverbatim
#  *>
#  *  =====================================================================
#        SUBROUTINE dgbtf2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
#  *
#  *  -- LAPACK computational routine (version 3.7.0) --
#  *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
#  *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
#  *     December 2016
#  *
#  *     .. Scalar Arguments ..
#        INTEGER            INFO, KL, KU, LDAB, M, N
#  *     ..
#  *     .. Array Arguments ..
#        INTEGER            IPIV( * )
#        DOUBLE PRECISION   AB( LDAB, * )
#  *     ..
#  *
#  *  =====================================================================
#  *
#  *     .. Parameters ..
#        DOUBLE PRECISION   ONE, ZERO
#        parameter( one = 1.0d+0, zero = 0.0d+0 )
#  *     ..
#  *     .. Local Scalars ..
#        INTEGER            I, J, JP, JU, KM, KV
#  *     ..
#  *     .. External Functions ..
#        INTEGER            IDAMAX
#        EXTERNAL           idamax
#  *     ..
#  *     .. External Subroutines ..
#        EXTERNAL           dger, dscal, dswap, xerbla
#  *     ..
#  *     .. Intrinsic Functions ..
#        INTRINSIC          max, min
#  *     ..
#  *     .. Executable Statements ..
#  *
function gbtrfgen!(kl, ku, m, 
    AB::Matrix{T}
) where {T<:Union{AbstractFloat,Complex{AbstractFloat},Rational,Complex{Rational}}}

    ipiv = collect(1:m)

    info = 0

#  *     KV is the number of superdiagonals in the factor U, allowing for
#  *     fill-in.
#  *
#        kv = ku + kl
    kv = ku + kl

#  *
#  *     Gaussian elimination with partial pivoting
#  *
#  *     Set fill-in elements in columns KU+2 to KV to zero.
#  *
#        DO 20 j = ku + 2, min( kv, n )
#           DO 10 i = kv - j + 2, kl
#              ab( i, j ) = zero
#     10    CONTINUE
#     20 CONTINUE
    borne=min(kv,m)
    for j = ku+2:borne
        for i = kv -j+2:kl
            AB[i,j] = 0
        end
    end
#  *
#  *     JU is the index of the last column affected by the current stage
#  *     of the factorization.
#  *
#        ju = 1
    ju = 1
#  *
#        DO 40 j = 1, min( m, n )
    for j=1:m
                
#  *
#  *        Set fill-in elements in column J+KV to zero.
#  *
#           IF( j+kv.LE.n ) THEN
#              DO 30 i = 1, kl
#                 ab( i, j+kv ) = zero
#     30       CONTINUE
#           END IF
        if j+kv <= m
            for i=1:kl
                AB[i,j+kv] = 0
            end
        end

#  *
#  *        Find pivot and test for singularity. KM is the number of
#  *        subdiagonal elements in the current column.
#  *
#           km = min( kl, m-j )
#           jp = idamax( km+1, ab( kv+1, j ), 1 )
#           ipiv( j ) = jp + j - 1
        km = min(kl, m-j)

        amax = abs(AB[kv+1, j])
        jp = 1
        for k=2:km+1
            v=abs(AB[kv+k, j])
            if v > amax
                jp = k
                amax = v
            end
        end
        ipiv[j] = j + jp - 1

#           IF( ab( kv+jp, j ).NE.zero ) THEN
#              ju = max( ju, min( j+ku+jp-1, n ) )

        if AB[kv+jp, j] != 0 
            ju = max( ju, min(j+ku+jp-1, m) )
#  *
#  *           Apply interchange to columns J to JU.
#  *
#              IF( jp.NE.1 )
#       $         CALL dswap( ju-j+1, ab( kv+jp, j ), ldab-1,
#       $                     ab( kv+1, j ), ldab-1 )
            if jp != 1
                for k=j:ju
                    i1 = kv + jp + j - k
                    i2 = kv + 1 + j - k
                    AB[i1, k], AB[i2, k] = AB[i2,k], AB[i1,k]
                end
            end
#  *
#              IF( km.GT.0 ) THEN
#  *
#  *              Compute multipliers.
#  *
#                 CALL dscal( km, one / ab( kv+1, j ), ab( kv+2, j ), 1 )
            if km > 0
                coef = AB[kv+1,j]
                for k=1:km
                    AB[kv+1+k,j] /= coef
                end
            end

#  *
#  *              Update trailing submatrix within the band.
#  *
#                 IF( ju.GT.j )
#       $            CALL dger( km, ju-j, -one, ab( kv+2, j ), 1,
#       $                       ab( kv, j+1 ), ldab-1, ab( kv+1, j+1 ),
#       $                       ldab-1 )
# SUBROUTINE dger(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
# DO 20 j = 1,n
# IF (y(jy).NE.zero) THEN
#     temp = alpha*y(jy)
#     DO 10 i = 1,m
#         a(i,j) = a(i,j) + x(i)*temp
# 10             CONTINUE
# END IF
# jy = jy + incy
# 20     CONTINUE
            for j2 = 1:(ju-j)
                coef = AB[kv+1-j2,j+j2]
                for i2 = 1:km
                    AB[kv+1+i2-j2,j+j2] -= coef*AB[kv+1+i2,j]
                end
            end
#              END IF
#           ELSE
#  *
#  *           If pivot is zero, set INFO to the index of the pivot
#  *           unless a zero pivot has already been found.
#  *
#              IF( info.EQ.0 )
#       $         info = j
        else
            if info == 0
                info = j
            end
        end
#           END IF
#     40 CONTINUE
#        RETURN
#  *
#  *     End of DGBTF2
#  *
#        END
    end
    return (AB, ipiv)
end
# *
#  *  -- LAPACK computational routine (version 3.7.0) --
#  *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
#  *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
#  *     December 2016
#  *
#  *     .. Scalar Arguments ..
#        CHARACTER          TRANS
#        INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
#  *     ..
#  *     .. Array Arguments ..
#        INTEGER            IPIV( * )
#        DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
#  *     ..
#  *
#  *  =====================================================================
#  *
#  *     .. Parameters ..
#        DOUBLE PRECISION   ONE
#        parameter( one = 1.0d+0 )
#  *     ..
#  *     .. Local Scalars ..
#        LOGICAL            LNOTI, NOTRAN
#        INTEGER            I, J, KD, L, LM
#  *     ..
#  *     .. External Functions ..
#        LOGICAL            LSAME
#        EXTERNAL           lsame
#  *     ..
#  *     .. External Subroutines ..
#        EXTERNAL           dgemv, dger, dswap, dtbsv, xerbla
#  *     ..
#  *     .. Intrinsic Functions ..
#        INTRINSIC          max, min
#  *     ..
#  *     .. Executable Statements ..
#  *
#  *     Test the input parameters.
#  *
function gbtrsgen!(
#    trans, 
    kl, 
    ku, 
    m, 
    AB::Matrix{T}, 
    ipiv, 
    B::Array{T}
)    where {T<:Union{AbstractFloat,Complex{AbstractFloat},Rational,Complex{Rational}}}
    info = 0

#        info = 0
#        notran = lsame( trans, 'N' )
#        IF( .NOT.notran .AND. .NOT.lsame( trans, 'T' ) .AND. .NOT.
#       $    lsame( trans, 'C' ) ) THEN
#           info = -1
#        ELSE IF( n.LT.0 ) THEN
#           info = -2
#        ELSE IF( kl.LT.0 ) THEN
#           info = -3
#        ELSE IF( ku.LT.0 ) THEN
#           info = -4
#        ELSE IF( nrhs.LT.0 ) THEN
#           info = -5
#        ELSE IF( ldab.LT.( 2*kl+ku+1 ) ) THEN
#           info = -7
#        ELSE IF( ldb.LT.max( 1, n ) ) THEN
#           info = -10
#        END IF
#        IF( info.NE.0 ) THEN
#           CALL xerbla( 'DGBTRS', -info )
#           RETURN
#        END IF
#  *
#  *     Quick return if possible
#  *
#        IF( n.EQ.0 .OR. nrhs.EQ.0 )
#       $   RETURN
#  *
#        kd = ku + kl + 1
#        lnoti = kl.GT.0
    kd = ku + kl + 1
#  *
#        IF( notran ) THEN
#  *
#  *        Solve  A*X = B.
#  *
#  *        Solve L*X = B, overwriting B with X.
#  *
#  *        L is represented as a product of permutations and unit lower
#  *        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
#  *        where each transformation L(i) is a rank-one modification of
#  *        the identity matrix.
#  *
#           IF( lnoti ) THEN
#              DO 10 j = 1, n - 1
#                 lm = min( kl, n-j )
#                 l = ipiv( j )
#                 IF( l.NE.j )
#       $            CALL dswap( nrhs, b( l, 1 ), ldb, b( j, 1 ), ldb )
#                 CALL dger( lm, nrhs, -one, ab( kd+1, j ), 1, b( j, 1 ),
#       $                    ldb, b( j+1, 1 ), ldb )
#     10       CONTINUE
#           END IF
# DO 20 j = 1,n
# IF (y(jy).NE.zero) THEN
#     temp = alpha*y(jy)
#     DO 10 i = 1,m
#         a(i,j) = a(i,j) + x(i)*temp
# 10             CONTINUE
# END IF
# jy = jy + incy
# 20     CONTINUE
    if kl > 0
        for j=1:m-1
            lm = min(kl, m-j)
            l = ipiv[j]
            if l != j
                B[l,:], B[j,:] = B[j,:], B[l,:]
            end
            for i=1:lm
                B[j+i,:] -= AB[kd+i,j]*B[j,:]
            end
        end
    end
#  *
#           DO 20 i = 1, nrhs
#  *
#  *           Solve U*X = B, overwriting B with X.
#  *
#              CALL dtbsv( 'Upper', 'No transpose', 'Non-unit', n, kl+ku,
#       $                  ab, ldab, b( 1, i ), 1 )
#     20    CONTINUE
#  *
# DO 20 j = n,1,-1
# IF (x(j).NE.zero) THEN
#     l = kplus1 - j
#     IF (nounit) x(j) = x(j)/a(kplus1,j)
#     temp = x(j)
#     DO 10 i = j - 1,max(1,j-k),-1
#         x(i) = x(i) - temp*a(l+i,j)
# 10                     CONTINUE
# END IF
# 20             CONTINUE
    k = kd-1
    for i1=1:size(B,2)
        for j=m:-1:1
            l = kd - j
            B[j,i1] /= AB[kd,j]
            for i = j-1:-1:max(1,j-k)
                B[i,i1] -= AB[l+i,j]*B[j,i1]
            end
        end
    end

    return B
end



#        ELSE
#  *
#  *        Solve A**T*X = B.
#  *
#           DO 30 i = 1, nrhs
#  *
#  *           Solve U**T*X = B, overwriting B with X.
#  *
#              CALL dtbsv( 'Upper', 'Transpose', 'Non-unit', n, kl+ku, ab,
#       $                  ldab, b( 1, i ), 1 )
#     30    CONTINUE
#  *
#  *        Solve L**T*X = B, overwriting B with X.
#  *
#           IF( lnoti ) THEN
#              DO 40 j = n - 1, 1, -1
#                 lm = min( kl, n-j )
#                 CALL dgemv( 'Transpose', lm, nrhs, -one, b( j+1, 1 ),
#       $                     ldb, ab( kd+1, j ), 1, one, b( j, 1 ), ldb )
#                 l = ipiv( j )
#                 IF( l.NE.j )
#       $            CALL dswap( nrhs, b( l, 1 ), ldb, b( j, 1 ), ldb )
#     40       CONTINUE
#           END IF
#        END IF
#        RETURN
#  *
#  *     End of DGBTRS
