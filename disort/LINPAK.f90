! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! $Rev: 42 $ $Date: 2014-11-0712:42:45 -0500 (Fri, 07 Nov 2014) $
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Call tree:
!
!    SGBCO
!       SASUM
!       SDOT
!       SAXPY
!       SGBFA
!           ISAMAX
!           SAXPY
!           SSCAL
!       SSCAL
!   SGBSL
!       SDOT
!       SAXPY
!   SGECO
!       SASUM
!       SDOT
!       SAXPY
!       SGEFA
!           ISAMAX
!           SAXPY
!           SSCAL
!       SSCAL
!   SGESL
!       SDOT
!       SAXPY
!   SSWAP
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module m_tenstr_disort_linpack
  contains

      SUBROUTINE SGBCO( ABD, LDA, N, ML, MU, IPVT, RCOND, Z )

!         Factors a real band matrix by Gaussian elimination
!         and estimates the condition of the matrix.
!
!         Revision date:  8/1/82
!         Author:  Moler, C. B. (U. of New Mexico)
!
!     If  RCOND  is not needed, SGBFA is slightly faster.
!     To solve  A*X = B , follow SBGCO by SGBSL.
!
!     input:
!
!        ABD     REAL(LDA, N)
!                contains the matrix in band storage.  The columns
!                of the matrix are stored in the columns of  ABD  and
!                the diagonals of the matrix are stored in rows
!                ML+1 through 2*ML+MU+1 of  ABD .
!                See the comments below for details.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!                LDA must be .GE. 2*ML + MU + 1 .
!
!        N       INTEGER
!                the order of the original matrix.
!
!        ML      INTEGER
!                number of diagonals below the main diagonal.
!                0 .LE. ML .LT. N .
!
!        MU      INTEGER
!                number of diagonals above the main diagonal.
!                0 .LE. MU .LT. N .
!                more efficient if  ML .LE. MU .
!
!     on return
!
!        ABD     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        RCOND   REAL
!                an estimate of the reciprocal condition of  A .
!                For the system  A*X = B , relative perturbations
!                in  A  and  B  of size  epsilon  may cause
!                relative perturbations in  X  of size  epsilon/RCOND .
!                If  RCOND  is so small that the logical expression
!                           1.0 + RCOND .EQ. 1.0
!                is true, then  A  may be singular to working
!                precision.  In particular,  RCOND  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        Z       REAL(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is close to a singular matrix, then  Z  is
!                an approximate null vector in the sense that
!                norm(a*z) = rcond*norm(a)*norm(z) .
!
!     Band storage
!
!           If  A  is a band matrix, the following program segment
!           will set up the input.
!
!                   ML = (band width below the diagonal)
!                   MU = (band width above the diagonal)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX(1, J-MU)
!                      I2 = MIN(N, J+ML)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
!           In addition, the first  ML  rows in  ABD  are used for
!           elements generated during the triangularization.
!           The total number of rows needed in  ABD  is  2*ML+MU+1 .
!           The  ML+MU by ML+MU  upper left triangle and the
!           ML by ML  lower right triangle are not referenced.
!
!     Example:  if the original matrix is
!
!           1112 13  0  0  0
!           2122 2324  0  0
!            032 3334 35  0
!            0  043 4445 46
!            0  0  054 5556
!            0  0  0  065 66
!
!      then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABD should contain
!
!            *  *  *  +  +  +  , * = not used
!            *  * 1324 3546  , + = used for pivoting
!            * 1223 3445 56
!           1122 3344 5566
!           2132 4354 65  *
!
! --------------------------------------------------------------------


!     .. Scalar Arguments ..

      INTEGER   LDA, ML, MU, N
      REAL      RCOND
!     ..
!     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      ABD( LDA, * ), Z( * )
!     ..
!     .. Local Scalars ..

      INTEGER   INFO, IS, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM
      REAL      ANORM, EK, S, SM, T, WK, WKM, YNORM
!     ..
!     .. External Functions ..

      !REAL      SASUM, SDOT
      !EXTERNAL  SASUM, SDOT
!     ..
!     .. External Subroutines ..

      !EXTERNAL  SAXPY, SGBFA, SSCAL
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC ABS, MAX, MIN, SIGN
!     ..


!                       ** compute 1-norm of A
      ANORM  = 0.0E0
      L  = ML + 1
      IS = L + MU

      DO 10 J = 1, N

         ANORM  = MAX( ANORM, SASUM( L,ABD( IS,J ),1 ) )

         IF( IS.GT.ML + 1 ) IS = IS - 1

         IF( J.LE.MU ) L  = L + 1

         IF( J.GE.N - ML ) L  = L - 1

   10 CONTINUE
!                                               ** factor

      CALL SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

!     RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))) .
!     estimate = norm(Z)/norm(Y) where  A*Z = Y  and  trans(A)*Y = E.
!     trans(A) is the transpose of A.  The components of E  are
!     chosen to cause maximum local growth in the elements of W  where
!     trans(U)*W = E.  The vectors are frequently rescaled to avoid
!     overflow.

!                     ** solve trans(U)*W = E
      EK = 1.0E0

      DO 20 J = 1, N
         Z( J ) = 0.0E0
   20 CONTINUE


      M  = ML + MU + 1
      JU = 0

      DO 50 K = 1, N

         IF( Z( K ).NE.0.0E0 ) EK = SIGN( EK, -Z( K ) )

         IF( ABS( EK - Z( K ) ).GT.ABS( ABD( M,K ) ) ) THEN

            S  = ABS( ABD( M,K ) ) / ABS( EK - Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            EK = S*EK

         END IF

         WK   = EK - Z( K )
         WKM  = -EK - Z( K )
         S    = ABS( WK )
         SM   = ABS( WKM )

         IF( ABD( M,K ).NE.0.0E0 ) THEN

            WK   = WK / ABD( M, K )
            WKM  = WKM / ABD( M, K )

         ELSE

            WK   = 1.0E0
            WKM  = 1.0E0

         END IF

         KP1  = K + 1
         JU   = MIN( MAX( JU,MU + IPVT( K ) ), N )
         MM   = M

         IF( KP1.LE.JU ) THEN

            DO 30 J = KP1, JU
               MM     = MM - 1
               SM     = SM + ABS( Z( J ) + WKM*ABD( MM,J ) )
               Z( J ) = Z( J ) + WK*ABD( MM, J )
               S      = S + ABS( Z( J ) )
   30       CONTINUE

            IF( S.LT.SM ) THEN

               T  = WKM - WK
               WK = WKM
               MM = M

               DO 40 J = KP1, JU
                  MM = MM - 1
                  Z( J ) = Z( J ) + T*ABD( MM, J )
   40          CONTINUE

            END IF

         END IF

         Z( K ) = WK

   50 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

!                         ** solve trans(L)*Y = W
      DO 60 KB = 1, N
         K  = N + 1 - KB
         LM = MIN( ML, N - K )

         IF( K.LT.N )&
     &       Z( K ) = Z( K ) + SDOT( LM, ABD( M+1, K ), 1, Z( K+1 ), 1 )

         IF( ABS( Z( K ) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

         END IF

         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T

   60 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

      YNORM  = 1.0E0
!                         ** solve L*V = Y
      DO 70 K = 1, N

         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T
         LM     = MIN( ML, N - K )

         IF( K.LT.N )&
     &       CALL SAXPY( LM, T, ABD( M+1, K ), 1, Z( K+1 ), 1 )

         IF( ABS( Z(K) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z(K) )

            CALL SSCAL( N, S, Z, 1 )

            YNORM  = S*YNORM

         END IF

   70 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

      YNORM  = S*YNORM

!                           ** solve  U*Z = W
      DO 80 KB = 1, N

         K  = N + 1 - KB

         IF( ABS( Z( K ) ).GT.ABS( ABD( M,K ) ) ) THEN

            S  = ABS( ABD( M,K ) ) / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            YNORM  = S*YNORM

         END IF

         IF( ABD( M,K ).NE.0.0E0 ) Z( K ) = Z( K ) / ABD( M, K )
         IF( ABD( M,K ).EQ.0.0E0 ) Z( K ) = 1.0E0

         LM = MIN( K, M ) - 1
         LA = M - LM
         LZ = K - LM
         T  = -Z( K )

         CALL SAXPY( LM, T, ABD( LA,K ), 1, Z( LZ ), 1 )

   80 CONTINUE
!                              ** make znorm = 1.0

      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

      YNORM  = S*YNORM
      IF( ANORM.NE.0.0E0 ) RCOND  = YNORM / ANORM
      IF( ANORM.EQ.0.0E0 ) RCOND  = 0.0E0

      END

      SUBROUTINE SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

!         Factors a real band matrix by elimination.

!         Revision date:  8/1/82
!         Author:  Moler, C. B. (U. of New Mexico)

!     SGBFA is usually called by SBGCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.

!     Input:  same as SGBCO

!     On return:

!        ABD,IPVT    same as SGBCO

!        INFO    INTEGER
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0.  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that SGBSL will divide by zero if
!                     called.  Use  RCOND  in SBGCO for a reliable
!                     indication of singularity.

!     (see SGBCO for description of band storage mode)

! ----------------------------------------------------------------


!     .. Scalar Arguments ..

      INTEGER   INFO, LDA, ML, MU, N
!     ..
!     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      ABD( LDA, * )
!     ..
!     .. Local Scalars ..

      INTEGER   I, I0, J, J0, J1, JU, JZ, K, KP1, L, LM, M, MM, NM1
      REAL      T
!     ..
!     .. External Functions ..

      !INTEGER   ISAMAX
      !EXTERNAL  ISAMAX
!     ..
!     .. External Subroutines ..

      !EXTERNAL  SAXPY, SSCAL
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC MAX, MIN
!     ..


      M    = ML + MU + 1
      INFO = 0
!                        ** zero initial fill-in columns
      J0 = MU + 2
      J1 = MIN( N, M ) - 1

      DO 20 JZ = J0, J1

         I0 = M + 1 - JZ

         DO 10 I = I0, ML
            ABD( I, JZ ) = 0.0E0
   10    CONTINUE

   20 CONTINUE

      JZ = J1
      JU = 0
!                       ** Gaussian elimination with partial pivoting
      NM1  = N - 1

      DO 50 K = 1, NM1

         KP1 = K + 1
!                                  ** zero next fill-in column
         JZ = JZ + 1

         IF( JZ.LE.N ) THEN

            DO 30 I = 1, ML
               ABD( I, JZ ) = 0.0E0
   30       CONTINUE

         END IF
!                                  ** find L = pivot index
         LM  = MIN( ML, N - K )
         L   = ISAMAX( LM + 1, ABD( M, K ), 1 ) + M - 1
         IPVT( K ) = L + K - M

         IF( ABD( L,K ).EQ.0.0E0 ) THEN
!                                      ** zero pivot implies this column
!                                      ** already triangularized
            INFO = K

         ELSE
!                                ** interchange if necessary
            IF( L.NE.M ) THEN

               T           = ABD( L, K )
               ABD( L, K ) = ABD( M, K )
               ABD( M, K ) = T
            END IF
!                                      ** compute multipliers
            T  = - 1.0E0 / ABD( M, K )

            CALL SSCAL( LM, T, ABD( M + 1,K ), 1 )

!                               ** row elimination with column indexing

            JU = MIN( MAX( JU,MU + IPVT( K ) ), N )
            MM = M

            DO 40 J = KP1, JU

               L  = L - 1
               MM = MM - 1
               T  = ABD( L, J )

               IF( L.NE.MM ) THEN

                  ABD( L, J ) = ABD( MM, J )
                  ABD( MM, J ) = T

               END IF

               CALL SAXPY( LM, T, ABD( M+1, K ), 1, ABD( MM+1, J ), 1)

   40       CONTINUE

         END IF

   50 CONTINUE


      IPVT( N ) = N
      IF( ABD( M,N ).EQ.0.0E0 ) INFO = N

      END

      SUBROUTINE SGBSL( ABD, LDA, N, ML, MU, IPVT, B, JOB )

!         Solves the real band system
!            A * X = B  or  transpose(A) * X = B
!         using the factors computed by SBGCO or SGBFA.

!         Revision date:  8/1/82
!         Author:  Moler, C. B. (U. of New Mexico)

!     Input:

!        ABD     REAL(LDA, N)
!                the output from SBGCO or SGBFA.

!        LDA     INTEGER
!                the leading dimension of the array  ABD .

!        N       INTEGER
!                the order of the original matrix.

!        ML      INTEGER
!                number of diagonals below the main diagonal.

!        MU      INTEGER
!                number of diagonals above the main diagonal.

!        IPVT    INTEGER(N)
!                the pivot vector from SBGCO or SGBFA.

!        B       REAL(N)
!                the right hand side vector.

!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  transpose(A)*X = B

!     On return

!        B       the solution vector  X

!     Error condition

!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically, this indicates singularity,
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if SBGCO has set RCOND .GT. 0.0
!        or SGBFA has set INFO .EQ. 0 .

!     To compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
!        10 continue

! --------------------------------------------------------

!     .. Scalar Arguments ..

      INTEGER   JOB, LDA, ML, MU, N
!     ..
!     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      ABD( LDA, * ), B( * )
!     ..
!     .. Local Scalars ..

      INTEGER   K, KB, L, LA, LB, LM, M, NM1
      REAL      T
!     ..
!     .. External Functions ..

      !REAL      SDOT
      !EXTERNAL  SDOT
!     ..
!     .. External Subroutines ..

      !EXTERNAL  SAXPY
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC MIN
!     ..


      M   = MU + ML + 1
      NM1 = N - 1

      IF( JOB.EQ.0 ) THEN
!                           ** solve  A * X = B

!                               ** first solve L*Y = B
         IF( ML.NE.0 ) THEN

            DO 10 K = 1, NM1

               LM = MIN( ML, N - K )
               L  = IPVT( K )
               T  = B( L )

               IF( L.NE.K ) THEN

                  B( L ) = B( K )
                  B( K ) = T

               END IF

               CALL SAXPY( LM, T, ABD( M + 1,K ), 1, B( K + 1 ), 1 )

   10       CONTINUE

         END IF

!                           ** now solve  U*X = Y
         DO 20 KB = 1, N

            K      = N + 1 - KB
            B( K ) = B( K ) / ABD( M, K )
            LM     = MIN( K, M ) - 1
            LA     = M - LM
            LB     = K - LM
            T      = -B( K )

            CALL SAXPY( LM, T, ABD( LA,K ), 1, B( LB ), 1 )

   20    CONTINUE


      ELSE
!                          ** solve  trans(A) * X = B

!                                  ** first solve  trans(U)*Y = B
         DO 30 K = 1, N

            LM     = MIN( K, M ) - 1
            LA     = M - LM
            LB     = K - LM
            T      = SDOT( LM, ABD( LA,K ), 1, B( LB ), 1 )
            B( K ) = ( B( K ) - T ) / ABD( M, K )

   30    CONTINUE

!                                  ** now solve trans(L)*X = Y
         IF( ML.NE.0 ) THEN

            DO 40 KB = 1, NM1

               K      = N - KB
               LM     = MIN( ML, N - K )
               B( K ) = B( K ) + SDOT( LM, ABD( M+1, K ), 1,&
     &                                 B( K+1 ), 1 )
               L      = IPVT( K )

               IF( L.NE.K ) THEN

                  T    = B( L )
                  B( L ) = B( K )
                  B( K ) = T

               END IF

   40       CONTINUE

         END IF

      END IF

      END

      SUBROUTINE SGECO( A, LDA, N, IPVT, RCOND, Z )

!         Factors a real matrix by Gaussian elimination
!         and estimates the condition of the matrix.

!         Revision date:  8/1/82
!         Author:  Moler, C. B. (U. of New Mexico)

!         If  RCOND  is not needed, SGEFA is slightly faster.
!         To solve  A*X = B , follow SGECO by SGESL.

!     On entry

!        A       REAL(LDA, N)
!                the matrix to be factored.

!        LDA     INTEGER
!                the leading dimension of the array  A .

!        N       INTEGER
!                the order of the matrix  A .

!     On return

!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U , where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.

!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.

!        RCOND   REAL
!                an estimate of the reciprocal condition of  A .
!                For the system  A*X = B , relative perturbations
!                in  A  and  B  of size  epsilon  may cause
!                relative perturbations in  X  of size  epsilon/RCOND .
!                If  RCOND  is so small that the logical expression
!                           1.0 + RCOND .EQ. 1.0
!                is true, then  A  may be singular to working
!                precision.  In particular,  RCOND  is zero  if
!                exact singularity is detected or the estimate
!                underflows.

!        Z       REAL(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is close to a singular matrix, then  Z  is
!                an approximate null vector in the sense that
!                norm(A*Z) = RCOND*norm(A)*norm(Z) .

! ------------------------------------------------------------------

!     .. Scalar Arguments ..

      INTEGER   LDA, N
      REAL      RCOND
!     ..
!     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      A( LDA, * ), Z( * )
!     ..
!     .. Local Scalars ..

      INTEGER   INFO, J, K, KB, KP1, L
      REAL      ANORM, EK, S, SM, T, WK, WKM, YNORM
!     ..
!     .. External Functions ..

      !REAL      SASUM, SDOT
      !EXTERNAL  SASUM, SDOT
!     ..
!     .. External Subroutines ..

      !EXTERNAL  SAXPY, SGEFA, SSCAL
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC ABS, MAX, SIGN
!     ..


!                        ** compute 1-norm of A
      ANORM  = 0.0E0
      DO 10 J = 1, N
         ANORM  = MAX( ANORM, SASUM( N,A( 1,J ),1 ) )
   10 CONTINUE
!                                      ** factor

      CALL SGEFA( A, LDA, N, IPVT, INFO )

!     RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))) .
!     estimate = norm(Z)/norm(Y) where  A*Z = Y  and  trans(A)*Y = E .
!     trans(A) is the transpose of A.  The components of E  are
!     chosen to cause maximum local growth in the elements of W  where
!     trans(U)*W = E.  The vectors are frequently rescaled to avoid
!     overflow.

!                        ** solve trans(U)*W = E
      EK = 1.0E0

      DO 20 J = 1, N
         Z( J ) = 0.0E0
   20 CONTINUE


      DO 50 K = 1, N

         IF( Z( K ).NE.0.0E0 ) EK = SIGN( EK, -Z( K ) )

         IF( ABS( EK - Z( K ) ).GT.ABS( A( K,K ) ) ) THEN

            S  = ABS( A( K,K ) ) / ABS( EK - Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            EK = S*EK

         END IF

         WK   = EK - Z( K )
         WKM  = -EK - Z( K )
         S    = ABS( WK )
         SM   = ABS( WKM )

         IF( A( K,K ).NE.0.0E0 ) THEN

            WK   = WK / A( K, K )
            WKM  = WKM / A( K, K )

         ELSE

            WK   = 1.0E0
            WKM  = 1.0E0

         END IF

         KP1  = K + 1

         IF( KP1.LE.N ) THEN

            DO 30 J = KP1, N
               SM     = SM + ABS( Z( J ) + WKM*A( K,J ) )
               Z( J ) = Z( J ) + WK*A( K, J )
               S      = S + ABS( Z( J ) )
   30       CONTINUE

            IF( S.LT.SM ) THEN

               T  = WKM - WK
               WK = WKM

               DO 40 J = KP1, N
                  Z( J ) = Z( J ) + T*A( K, J )
   40          CONTINUE

            END IF

         END IF

         Z( K ) = WK

   50 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )
!                                ** solve trans(L)*Y = W
      DO 60 KB = 1, N
         K  = N + 1 - KB

         IF( K.LT.N )&
     &       Z( K ) = Z( K ) + SDOT( N - K, A( K+1, K ), 1, Z( K+1 ), 1)

         IF( ABS( Z( K ) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

         END IF

         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T
   60 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )
!                                 ** solve L*V = Y
      YNORM  = 1.0E0

      DO 70 K = 1, N
         L      = IPVT( K )
         T      = Z( L )
         Z( L ) = Z( K )
         Z( K ) = T

         IF( K.LT.N ) CALL SAXPY( N - K, T, A( K + 1,K ), 1, Z( K + 1 ),&
     &                            1 )

         IF( ABS( Z( K ) ).GT.1.0E0 ) THEN

            S  = 1.0E0 / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            YNORM  = S*YNORM
         END IF

   70 CONTINUE


      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )
!                                  ** solve  U*Z = V
      YNORM  = S*YNORM

      DO 80 KB = 1, N

         K  = N + 1 - KB

         IF( ABS( Z( K ) ).GT.ABS( A( K,K ) ) ) THEN

            S  = ABS( A( K,K ) ) / ABS( Z( K ) )

            CALL SSCAL( N, S, Z, 1 )

            YNORM  = S*YNORM

         END IF

         IF( A( K,K ).NE.0.0E0 ) Z( K ) = Z( K ) / A( K, K )

         IF( A( K,K ).EQ.0.0E0 ) Z( K ) = 1.0E0

         T  = -Z( K )

         CALL SAXPY( K - 1, T, A( 1,K ), 1, Z( 1 ), 1 )

   80 CONTINUE
!                                   ** make znorm = 1.0
      S  = 1.0E0 / SASUM( N, Z, 1 )

      CALL SSCAL( N, S, Z, 1 )

      YNORM  = S*YNORM

      IF( ANORM.NE.0.0E0 ) RCOND = YNORM / ANORM
      IF( ANORM.EQ.0.0E0 ) RCOND = 0.0E0

      END

      SUBROUTINE SGEFA( A, LDA, N, IPVT, INFO )

!         Factors a real matrix by Gaussian elimination.

!         Revision date:  8/1/82
!         Author:  Moler, C. B. (U. of New Mexico)

!     SGEFA is usually called by SGECO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (time for SGECO) = (1 + 9/N) * (time for SGEFA) .

!     Input:  same as SGECO

!     On return:

!        A,IPVT  same as SGECO

!        INFO    INTEGER
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0.  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that SGESL or SGEDI will divide by zero
!                     if called.  Use  RCOND  in SGECO for a reliable
!                     indication of singularity.

! ---------------------------------------------------------------------

!     .. Scalar Arguments ..

      INTEGER   INFO, LDA, N
!     ..
!     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      A( LDA, * )
!     ..
!     .. Local Scalars ..

      INTEGER   J, K, KP1, L, NM1
      REAL      T
!     ..
!     .. External Functions ..

      !INTEGER   ISAMAX
      !EXTERNAL  ISAMAX
!     ..
!     .. External Subroutines ..

      !EXTERNAL  SAXPY, SSCAL
!     ..


!                      ** Gaussian elimination with partial pivoting
      INFO = 0
      NM1  = N - 1

      DO 20 K = 1, NM1

         KP1  = K + 1
!                                            ** find L = pivot index

         L  = ISAMAX( N - K + 1, A( K,K ), 1 ) + K - 1
         IPVT( K ) = L

         IF( A( L,K ).EQ.0.0E0 ) THEN
!                                     ** zero pivot implies this column
!                                     ** already triangularized
            INFO = K

         ELSE
!                                     ** interchange if necessary
            IF( L.NE.K ) THEN

               T         = A( L, K )
               A( L, K ) = A( K, K )
               A( K, K ) = T

            END IF
!                                     ** compute multipliers
            T  = -1.0E0 / A( K, K )

            CALL SSCAL( N - K, T, A( K + 1,K ), 1 )

!                              ** row elimination with column indexing
            DO 10 J = KP1, N

               T  = A( L, J )

               IF( L.NE.K ) THEN

                  A( L, J ) = A( K, J )
                  A( K, J ) = T

               END IF

               CALL SAXPY( N-K, T, A( K+1, K ), 1, A( K+1, J ), 1 )

   10       CONTINUE

         END IF

   20 CONTINUE


      IPVT( N ) = N
      IF( A( N,N ) .EQ. 0.0E0 ) INFO = N

      END

      SUBROUTINE SGESL( A, LDA, N, IPVT, B, JOB )

!         Solves the real system
!            A * X = B  or  transpose(A) * X = B
!         using the factors computed by SGECO or SGEFA.

!         Revision date:  8/1/82
!         Author:  Moler, C. B. (U. of New Mexico)

!     On entry

!        A       REAL(LDA, N)
!                the output from SGECO or SGEFA.

!        LDA     INTEGER
!                the leading dimension of the array  A

!        N       INTEGER
!                the order of the matrix  A

!        IPVT    INTEGER(N)
!                the pivot vector from SGECO or SGEFA.

!        B       REAL(N)
!                the right hand side vector.

!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  transpose(A)*X = B

!     On return

!        B       the solution vector  X

!     Error condition

!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically, this indicates singularity,
!        but it is often caused by improper arguments or improper
!        setting of LDA.  It will not occur if the subroutines are
!        called correctly and if SGECO has set RCOND .GT. 0.0
!        or SGEFA has set INFO .EQ. 0 .

!     To compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call sgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call sgesl(a,lda,n,ipvt,c(1,j),0)
!        10 continue

! ---------------------------------------------------------------------

!     .. Scalar Arguments ..

      INTEGER   JOB, LDA, N
!     ..
!     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      A( LDA, * ), B( * )
!     ..
!     .. Local Scalars ..

      INTEGER   K, KB, L, NM1
      REAL      T
!     ..
!     .. External Functions ..

      !REAL      SDOT
      !EXTERNAL  SDOT
!     ..
!     .. External Subroutines ..

      !EXTERNAL  SAXPY
!     ..


      NM1  = N - 1

      IF( JOB.EQ.0 ) THEN
!                                 ** solve  A * X = B

!                                     ** first solve  L*Y = B
         DO 10 K = 1, NM1

            L  = IPVT( K )
            T  = B( L )

            IF( L.NE.K ) THEN

               B( L ) = B( K )
               B( K ) = T

            END IF

            CALL SAXPY( N - K, T, A( K+1, K ), 1, B( K+1 ), 1 )

   10    CONTINUE
!                                    ** now solve  U*X = Y
         DO 20 KB = 1, N

            K      = N + 1 - KB
            B( K ) = B( K ) / A( K, K )
            T      = - B( K )

            CALL SAXPY( K-1, T, A( 1, K ), 1, B(1), 1 )

   20    CONTINUE


      ELSE
!                         ** solve  trans(A) * X = B

!                                    ** first solve  trans(U)*Y = B
         DO 30 K = 1, N

            T      = SDOT( K - 1, A( 1,K ), 1, B( 1 ), 1 )
            B( K ) = ( B( K ) - T ) / A( K, K )

   30    CONTINUE

!                                    ** now solve  trans(l)*x = y
         DO 40 KB = 1, NM1

            K      = N - KB
            B( K ) = B( K ) + SDOT( N - K, A( K+1, K ), 1, B( K+1 ), 1)
            L      = IPVT( K )

            IF( L.NE.K ) THEN

               T      = B( L )
               B( L ) = B( K )
               B( K ) = T

            END IF

   40    CONTINUE

      END IF

      END

      REAL FUNCTION SASUM( N, SX, INCX )

!  INPUT--    N  Number of elements in vector to be summed
!            SX  Sing-prec array, length 1+(N-1)*INCX, containing vector
!          INCX  Spacing of vector elements in SX

!  OUTPUT-- SASUM   Sum from 0 to N-1 of  ABS(SX(1+I*INCX))
! ----------------------------------------------------------

!     .. Scalar Arguments ..

      INTEGER   INCX, N
!     ..
!     .. Array Arguments ..

      REAL      SX( * )
!     ..
!     .. Local Scalars ..

      INTEGER   I, M
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC ABS, MOD
!     ..

      SASUM  = 0.0

      IF( N.LE.0 ) RETURN

      IF( INCX.NE.1 ) THEN
!                                          ** non-unit increments
         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            SASUM  = SASUM + ABS( SX( I ) )
   10    CONTINUE

      ELSE
!                                          ** unit increments
         M  = MOD( N, 6 )

         IF( M.NE.0 ) THEN
!                             ** clean-up loop so remaining vector
!                             ** length is a multiple of 6.
            DO 20 I = 1, M
               SASUM  = SASUM + ABS( SX( I ) )
   20       CONTINUE

         END IF
!                              ** unroll loop for speed
         DO 30 I = M + 1, N, 6
            SASUM  = SASUM + ABS( SX( I ) ) + ABS( SX( I + 1 ) ) +&
     &               ABS( SX( I + 2 ) ) + ABS( SX( I + 3 ) ) +&
     &               ABS( SX( I + 4 ) ) + ABS( SX( I + 5 ) )
   30    CONTINUE

      END IF

      END

      SUBROUTINE SAXPY( N, SA, SX, INCX, SY, INCY )

!          Y = A*X + Y  (X, Y = vectors, A = scalar)

!  INPUT--
!        N  Number of elements in input vectors X and Y
!       SA  Single precision scalar multiplier A
!       SX  Sing-prec array containing vector X
!     INCX  Spacing of elements of vector X in SX
!       SY  Sing-prec array containing vector Y
!     INCY  Spacing of elements of vector Y in SY

! OUTPUT--
!       SY   For I = 0 to N-1, overwrite  SY(LY+I*INCY) with
!                 SA*SX(LX+I*INCX) + SY(LY+I*INCY),
!            where LX = 1          if INCX .GE. 0,
!                     = (-INCX)*N  if INCX .LT. 0
!            and LY is defined analogously using INCY.
! ------------------------------------------------------------

!     .. Scalar Arguments ..

      INTEGER   INCX, INCY, N
      REAL      SA
!     ..
!     .. Array Arguments ..

      REAL      SX( * ), SY( * )
!     ..
!     .. Local Scalars ..

      INTEGER   I, IX, IY, M
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC MOD
!     ..


      IF( N.LE.0 .OR. SA.EQ.0.0 ) RETURN

      IF( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN

         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            SY( I ) = SY( I ) + SA*SX( I )
   10    CONTINUE

      ELSE IF( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN

!                                        ** equal, unit increments
         M  = MOD( N, 4 )

         IF( M.NE.0 ) THEN
!                            ** clean-up loop so remaining vector length
!                            ** is a multiple of 4.
            DO 20 I = 1, M
               SY( I ) = SY( I ) + SA*SX( I )
   20       CONTINUE

         END IF
!                              ** unroll loop for speed
         DO 30 I = M + 1, N, 4
            SY( I ) = SY( I ) + SA*SX( I )
            SY( I + 1 ) = SY( I + 1 ) + SA*SX( I + 1 )
            SY( I + 2 ) = SY( I + 2 ) + SA*SX( I + 2 )
            SY( I + 3 ) = SY( I + 3 ) + SA*SX( I + 3 )
   30    CONTINUE


      ELSE
!               ** nonequal or nonpositive increments.
         IX = 1
         IY = 1
         IF( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )

         DO 40 I = 1, N
            SY( IY ) = SY( IY ) + SA*SX( IX )
            IX = IX + INCX
            IY = IY + INCY
   40    CONTINUE

      END IF

      END

      REAL FUNCTION SDOT( N, SX, INCX, SY, INCY )

!        Single-prec dot product of vectors  X  and  Y

!  INPUT--
!        N  Number of elements in input vectors X and Y
!       SX  Sing-prec array containing vector X
!     INCX  Spacing of elements of vector X in SX
!       SY  Sing-prec array containing vector Y
!     INCY  Spacing of elements of vector Y in SY

! OUTPUT--
!     SDOT   Sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
!            where  LX = 1          if INCX .GE. 0,
!                      = (-INCX)*N  if INCX .LT. 0,
!            and LY is defined analogously using INCY.
! ------------------------------------------------------------------

!     .. Scalar Arguments ..

      INTEGER   INCX, INCY, N
!     ..
!     .. Array Arguments ..

      REAL      SX( * ), SY( * )
!     ..
!     .. Local Scalars ..

      INTEGER   I, IX, IY, M
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC MOD
!     ..


      SDOT = 0.0

      IF( N.LE.0 ) RETURN

      IF( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN

         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            SDOT = SDOT + SX( I )*SY( I )
   10    CONTINUE


      ELSE IF( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN

!                                        ** equal, unit increments
         M  = MOD( N, 5 )

         IF( M.NE.0 ) THEN
!                            ** clean-up loop so remaining vector length
!                            ** is a multiple of 4.
            DO 20 I = 1, M
               SDOT = SDOT + SX( I )*SY( I )
   20       CONTINUE

         END IF
!                              ** unroll loop for speed
         DO 30 I = M + 1, N, 5
            SDOT = SDOT + SX( I )*SY( I ) + SX( I + 1 )*SY( I + 1 ) +&
     &               SX( I + 2 )*SY( I + 2 ) + SX( I + 3 )*SY( I + 3 ) +&
     &               SX( I + 4 )*SY( I + 4 )
   30    CONTINUE

      ELSE
!               ** nonequal or nonpositive increments.
         IX = 1
         IY = 1

         IF( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )

         DO 40 I = 1, N
            SDOT = SDOT + SX( IX )*SY( IY )
            IX   = IX + INCX
            IY   = IY + INCY
   40    CONTINUE

      END IF

      END

      SUBROUTINE SSCAL( N, SA, SX, INCX )

!         Multiply vector SX by scalar SA

!  INPUT--  N  Number of elements in vector
!          SA  Single precision scale factor
!          SX  Sing-prec array, length 1+(N-1)*INCX, containing vector
!        INCX  Spacing of vector elements in SX

! OUTPUT-- SX  Replace  SX(1+I*INCX)  with  SA * SX(1+I*INCX)
!                for I = 0 to N-1
! ---------------------------------------------------------------------

!     .. Scalar Arguments ..

      INTEGER   INCX, N
      REAL      SA
!     ..
!     .. Array Arguments ..

      REAL      SX( * )
!     ..
!     .. Local Scalars ..

      INTEGER   I, M
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC MOD
!     ..


      IF( N.LE.0 ) RETURN

      IF( INCX.NE.1 ) THEN

         DO 10 I = 1, 1 + ( N - 1 )*INCX, INCX
            SX( I ) = SA*SX( I )
   10    CONTINUE


      ELSE

         M  = MOD( N, 5 )

         IF( M.NE.0 ) THEN
!                           ** clean-up loop so remaining vector length
!                           ** is a multiple of 5.
            DO 20 I = 1, M
               SX( I ) = SA*SX( I )
   20       CONTINUE

         END IF
!                             ** unroll loop for speed
         DO 30 I = M + 1, N, 5
            SX( I ) = SA*SX( I )
            SX( I + 1 ) = SA*SX( I + 1 )
            SX( I + 2 ) = SA*SX( I + 2 )
            SX( I + 3 ) = SA*SX( I + 3 )
            SX( I + 4 ) = SA*SX( I + 4 )
   30    CONTINUE

      END IF

      END

      SUBROUTINE SSWAP( N, SX, INCX, SY, INCY )

!          Interchange s.p vectors  X  and  Y, as follows:

!     For I = 0 to N-1, interchange  SX(LX+I*INCX) and SY(LY+I*INCY),
!     where LX = 1          if INCX .GE. 0,
!              = (-INCX)*N  if INCX .LT. 0
!     and LY is defined analogously using INCY.


!  INPUT--
!        N  Number of elements in input vectors X and Y
!       SX  Sing-prec array containing vector X
!     INCX  Spacing of elements of vector X in SX
!       SY  Sing-prec array containing vector Y
!     INCY  Spacing of elements of vector Y in SY

! OUTPUT--
!       SX  Input vector SY (unchanged if N .LE. 0)
!       SY  Input vector SX (unchanged IF N .LE. 0)
! --------------------------------------------------------------

!     .. Scalar Arguments ..

      INTEGER   INCX, INCY, N
!     ..
!     .. Array Arguments ..

      REAL      SX( * ), SY( * )
!     ..
!     .. Local Scalars ..

      INTEGER   I, IX, IY, M
      REAL      STEMP1, STEMP2, STEMP3
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC MOD
!     ..


      IF( N.LE.0 ) RETURN

      IF( INCX.EQ.INCY .AND. INCX.GT.1 ) THEN

         DO 10 I = 1, 1 + ( N-1 )*INCX, INCX
            STEMP1 = SX( I )
            SX( I ) = SY( I )
            SY( I ) = STEMP1
   10    CONTINUE


      ELSE IF( INCX.EQ.INCY .AND. INCX.EQ.1 ) THEN

!                                        ** equal, unit increments
         M  = MOD( N, 3 )

         IF( M.NE.0 ) THEN
!                            ** clean-up loop so remaining vector length
!                            ** is a multiple of 3.
            DO 20 I = 1, M
               STEMP1 = SX( I )
               SX( I ) = SY( I )
               SY( I ) = STEMP1
   20       CONTINUE

         END IF
!                              ** unroll loop for speed
         DO 30 I = M + 1, N, 3
            STEMP1 = SX( I )
            STEMP2 = SX( I + 1 )
            STEMP3 = SX( I + 2 )
            SX( I ) = SY( I )
            SX( I + 1 ) = SY( I + 1 )
            SX( I + 2 ) = SY( I + 2 )
            SY( I ) = STEMP1
            SY( I + 1 ) = STEMP2
            SY( I + 2 ) = STEMP3
   30    CONTINUE


      ELSE
!               ** nonequal or nonpositive increments.
         IX = 1
         IY = 1

         IF( INCX.LT.0 ) IX = 1 + ( N - 1 )*( -INCX )
         IF( INCY.LT.0 ) IY = 1 + ( N - 1 )*( -INCY )

         DO 40 I = 1, N
            STEMP1 = SX( IX )
            SX( IX ) = SY( IY )
            SY( IY ) = STEMP1
            IX   = IX + INCX
            IY   = IY + INCY
   40    CONTINUE

      END IF

      END

      INTEGER FUNCTION ISAMAX( N, SX, INCX )

! INPUT--  N     Number of elements in vector of interest
!          SX    Sing-prec array, length 1+(N-1)*INCX, containing vector
!          INCX  Spacing of vector elements in SX

! OUTPUT-- ISAMAX   First I, I = 1 to N, to maximize
!                         ABS(SX(1+(I-1)*INCX))
! ---------------------------------------------------------------------

!     .. Scalar Arguments ..

      INTEGER   INCX, N
!     ..
!     .. Array Arguments ..

      REAL      SX( * )
!     ..
!     .. Local Scalars ..

      INTEGER   I, II
      REAL      SMAX, XMAG
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC ABS
!     ..

      ISAMAX = 0
      
      IF( N.LE.0 ) THEN

         ISAMAX = 0

      ELSE IF( N.EQ.1 ) THEN

         ISAMAX = 1

      ELSE

         SMAX = 0.0
         II   = 1

         DO 10 I = 1, 1 + ( N-1 )*INCX, INCX

            XMAG = ABS( SX( I ) )

            IF( SMAX.LT.XMAG ) THEN

               SMAX   = XMAG
               ISAMAX = II

            END IF

            II = II + 1

   10    CONTINUE

      END IF

      END FUNCTION ISAMAX

end module
