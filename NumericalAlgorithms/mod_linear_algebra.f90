    module mod_linear_algebra
    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY : IEEE_FMA 
    
!** FMACHP is the machine constant                                      
!** (i.e., the smallest positive machine number for which  1+FMACHP > 1 
    DOUBLEPRECISION, PARAMETER :: FMACHP = EPSILON(1.0D0) 
    
    INTERFACE GAUSS
        MODULE PROCEDURE :: GAUSS_VEC, GAUSS_MAT
    END INTERFACE GAUSS
    
    contains
    
    FUNCTION GAUSS_VEC(A, Y) RESULT(X)
    DOUBLEPRECISION, INTENT(in) :: A(:,:), Y(:)
    DOUBLEPRECISION X(SIZE(A,2))
    
    DOUBLEPRECISION LR(SIZE(A,1),SIZE(A,2)), D(SIZE(A,1))
    INTEGER N, M, MARK, IPIVOT(SIZE(A,2))
    
        N = SIZE(A, 1)
        M = SIZE(A, 2)
        
        LR = A
        
        CALL GAUSSA( M, LR, N, Y, X, MARK, D, IPIVOT )
    
    END FUNCTION
    
    FUNCTION GAUSS_MAT(A, Y) RESULT(X)
    DOUBLEPRECISION, INTENT(in) :: A(:,:), Y(:,:)
    DOUBLEPRECISION X(SIZE(A,2), SIZE(Y,2))
    
    DOUBLEPRECISION LR(SIZE(A,1),SIZE(A,2)), D(SIZE(A,1))
    INTEGER N, M, K, MARK, IPIVOT(SIZE(A,2))
    
        N = SIZE(A, 1)
        M = SIZE(A, 2)
        K = SIZE(Y, 2)
        
        LR = A
        
        CALL GAUSRS( M, LR, N, K, Y, X, MARK, D, IPIVOT )
    
    END FUNCTION
      
![KA{P 4}{Direct Methods for Solving Linear Systems}                    
![       {Direct Methods for Solving Systems of Linear                  
![        Equations}*)                                                  
![           {Gau"s Algorithm with Column Pivot Search}*)               
      SUBROUTINE GAUSSA(N, A, LDA, Y, X, MARK, D, IPIVOT) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Solving a linear system of equations  A * X = Y  by applying  *      
!  the Gauss-elimination method with scaling and column pivot    *      
!  search.                                                       *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N        : order of the linear system.                        *      
!  A        : 2-dimensional array A(1:LDA,1:N); the matrix A is  *      
!             the system matrix of the equations, (A = A(ORG)).  *      
!  LDA      : leading dimension of A as defined in the calling   *      
!             program.                                           *      
!  Y        : N-vector Y(1:N); the right hand side of the system *      
!             of equations.                                      *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  A        : 2-dimensional array A(1:LDA,1:N), containing the   *      
!             factors L and R with P * A(ORG) = L * R.           *      
!             P = permutation matrix, L = unit lower triangular  *      
!             matrix, and R = upper triangular matrix.           *      
!  X        : N-vector X(1:N); the solution vector of the system *      
!             of equations.                                      *      
!  MARK     : = 1, even number of row permutations.              *      
!             =-1, odd number of row permutations.               *      
!             = 0, input array A is numerically singular.        *      
!             The determinant of A can be computed as :          *      
!                DET(A(ORG)) = MARK * A(1,1) * ... * A(N,N).     *      
!  D        : N-vector D(1:N); the reciprocals of the row sum    *      
!             norms of A(ORG) that serve as scaling factors:     *      
!             D(I) = 1./(ABS(A(I,1)) + ... + ABS(A(I,N)))  for   *      
!             I = 1, ..., N.                                     *      
!  IPIVOT   : N-vector IPIVOT(1:N); it indicates the row         *      
!             permutations for the scaled column pivot search    *      
!             and thereby defines the permutation matrix P.      *      
!             If e.g. IPIVOT(2) = 7, then the 7th row of A(ORG)  *      
!             is permuted to become the 2nd row of P * A(ORG).   *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: GAUSSP, GAUSSS                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  authors   : Gisela Engeln-Muellges, Guido Dubois              *      
!  date      : 04.25.88                                          *      
!  source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER LDA, N, MARK, IPIVOT (1:N)
      DOUBLEPRECISION A, Y, X, D
      DIMENSION A (1:LDA, 1:N), Y (1:N), X (1:N), D (1:N) 
!                                                                       
!  Factor the matrix A by using SUBROUTINE GAUSSP.                      
!                                                                       
      CALL GAUSSP (N, A, LDA, IPIVOT, MARK, D) 
!                                                                       
!  Updating and backsubstitution via SUBROUTINE GAUSSS                  
!  in order to find the solution of the system of equations.            
!                                                                       
      IF (MARK.NE.0) CALL GAUSSS (N, A, LDA, IPIVOT, Y, X) 
      RETURN 
      END SUBROUTINE GAUSSA    
      
      SUBROUTINE GAUSRS (N, A, LDA, M, RS, XL, MARK, D, IPIVOT) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Solving a linear systems of equations  A * XL = RS  for M     *      
!  right hand sides using the Gauss-elimination method with      *      
!  scaling and column pivot search .                             *      
!  If the system has the form                                    *      
!         A * XL = I  , where I = identity matrix and A, XL, I   *      
!  are all (NxN) matrices, then the solution XL is the matrix    *      
!  inverse of A.                                                 *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETER:                                              *      
!  ================                                              *      
!  N        : order of the system of equations.                  *      
!  A        : 2-dimensional array A(1:LDA,1:N), containing the   *      
!             LDAxN matrix A common to all M systems of equations*      
!             (A = A(ORG)).                                      *      
!  LDA      : leading dimension of A, RS and XL, as defined in   *      
!             the calling program.                               *      
!  M        : number of right hand sides and hence the number of *      
!             solution vectors.                                  *      
!  RS       : 2-dimensional array RS(1:LDA,1:M), that is formed  *      
!             with the M right hand sides as columns.            *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  A        : 2-dimensional array A(1:LDA,1:N), containing the   *      
!             factors L and R with  P * A(ORG) = L * R. Here     *      
!             P = permutation matrix, L = unit lower triangular  *      
!             matrix and R = upper triangular matrix.            *      
!  XL       : 2-dimensional array XL(1:LDA,1:M) that contains    *      
!             the M solution vectors as columns for each of the  *      
!             M systems of equations.                            *      
!  MARK     : = 1, even number of row permutations.              *      
!             =-1, odd number of row permutations.               *      
!             = 0, input matrix A is numerically singular.       *      
!             The determinant of A is :                          *      
!                DET(A(ORG)) = MARK * A(1,1) * ... * A(N,N).     *      
!  D        : N-vector D(1:N); the reciprocals of the row sum    *      
!             norms of A(ORG), used for scaling:                 *      
!             D(I) = 1./(ABS(A(I,1)) + ... + ABS(A(I,N)))  for   *      
!             I = 1, ..., N.                                     *      
!  IPIVOT   : N-vector IPIVOT(1:N); it indicates the row per-    *      
!             mutations and thus defines the permutation matrix  *      
!             P. If e.g. IPIVOT(2) = 7, then the 7th row in      *      
!             of A(ORG) will become the 2nd row of P * A(ORG).   *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: GAUSSP, GAUSSS                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  authors   : Gisela Engeln-Muellges, Guido Dubois              *      
!  date      : 04.25.88                                          *      
!  source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER LDA, N, M, K, MARK
      DOUBLEPRECISION A, RS, XL, D
      DIMENSION A (1:LDA, 1:N), RS (1:LDA, 1:M), XL (1:LDA, 1:M)
      DIMENSION D (1:N)                                                           
      INTEGER IPIVOT (1:N) 
!                                                                       
!  Factoring the matrix A by applying SUBROUTINE GAUSSP.                
!                                                                       
      CALL GAUSSP (N, A, LDA, IPIVOT, MARK, D) 
!                                                                       
!  Updating and bachsubstitution using SUBROUTINE GAUSSS in order to    
!  calculate the solution vectors for the M systems of equations.       
!                                                                       
      IF (MARK.NE.0) THEN 
         DO 10 K = 1, M 
            CALL GAUSSS (N, A, LDA, IPIVOT, RS (1, K), XL (1, K) ) 
   10    END DO 
      ENDIF 
      RETURN 
      END SUBROUTINE GAUSRS                         
      
!                                                                       
!                                                                       
      SUBROUTINE GAUSSP (N, A, LDA, IPIVOT, MARK, D) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Factoring the matrix A into the product of two matrices L and *      
!  R so that  P * A = L * R, where P = permutation matrix,       *      
!  L = unit lower triangular matrix and R = upper triangular     *      
!  matrix by applying the Gauss-elimination method with          *      
!  scaling and column pivot search.                              *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N        : order of the system of equations.                  *      
!  A        : 2-dimensional array A(1:LDA,1:N); the system matrix*      
!             of the system of equations, (A = A(ORG)).          *      
!  LDA      : leading dimension of A as defined in the calling   *      
!             program.                                           *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  A        : 2-dimensional array A(1:LDA,1:N), containing the   *      
!             factors L and R with  P * A(ORG) = L * R  for      *      
!             P a permutation matrix. The upper triangular R     *      
!             is stored in the upper triangle of A. The unit     *      
!             lower triangular matrix L, except for the diagonal *      
!             ones, is stored in the lower triangle of A.        *      
!  IPIVOT   : N-vector IPIVOT(1:N); it indicates the row         *      
!             permutations of the scaled column pivot search     *      
!             algorithm and thus defines the permutation matrix  *      
!             P. If e.g. IPIVOT(2) = 7, then the 7th row of      *      
!             A(ORG) has become the 2nd row of P * A(ORG).       *      
!  MARK     : = 1, even number of row permutations.              *      
!             =-1, odd number of row permutations.               *      
!             = 0, system matrix A is numerically singular.      *      
!             The determinant of A is :                          *      
!                DET(A(ORG)) = MARK * A(1,1) * ... * A(N,N).     *      
!  D        : N-vector D(1:N); the reciprocals of the row sum    *      
!             norms of A(ORG) that serve as scaling factors:     *      
!             D(I) = 1./(ABS(A(I,1)) + ... + ABS(A(I,N)))  for   *      
!             I = 1, ..., N.                                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required:                                         *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  authors   : Gisela Engeln-Muellges, Guido Dubois              *      
!  date      : 04.25.88                                          *      
!  source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER N, J, I, K, MARK, LDA, IPVT
      INTEGER IPIVOT (1:N)
      DOUBLEPRECISION A, D, RELERR, SUM, FAK, DUMMY, PIVOT
      DIMENSION A (1:LDA, 1:N), D (1:N)
!                                                                       
!  Local storage of error parameter RELERR in case that the             
!  SUBROUTINE is called repeatedly.                                     
!                                                                       
      MARK = 1 
!                                                                       
!  Calculation of the machine constant and initializing the relative    
!  error.                                                               
!                                                                       
      
      RELERR = 8.00D0 * FMACHP
!                                                                       
!  Calculation of row sum norms of A and initializing                   
!  the PIVOT vector.                                                    
!                                                                       
      DO 20 I = 1, N 
         IPIVOT (I) = I 
         SUM = ABS (A (I, 1) ) 
         DO 30 K = 2, N 
            SUM = SUM + ABS (A (I, K) ) 
   30    END DO 
         IF (SUM.EQ.0.0D0) THEN 
            MARK = 0 
            RETURN 
         ELSE 
            D (I) = 1.0D0 / SUM 
         ENDIF 
   20 END DO 
      IF (N.EQ.1) RETURN 
!                                                                       
!  Triangular factorization.                                            
!                                                                       
      DO 40 I = 1, N - 1 
!                                                                       
!  Determine the pivot row.                                             
!                                                                       
         PIVOT = ABS (A (I, I) ) * D (I) 
         IPVT = I 
         DO 50 J = I + 1, N 
            DUMMY = ABS (A (J, I) ) * D (J) 
            IF (DUMMY.GT.PIVOT) THEN 
               PIVOT = DUMMY 
               IPVT = J 
            ENDIF 
   50    END DO 
         IF (PIVOT.LT.RELERR) THEN 
            MARK = 0 
            RETURN 
         ELSE 
            IF (IPVT.NE.I) THEN 
!                                                                       
!  Interchange the I-th and the IPVT-th row of A.                       
!                                                                       
               MARK = - MARK 
               J = IPIVOT (I) 
               IPIVOT (I) = IPIVOT (IPVT) 
               IPIVOT (IPVT) = J 
               DUMMY = D (I) 
               D (I) = D (IPVT) 
               D (IPVT) = DUMMY 
               DO 60 J = 1, N 
                  DUMMY = A (I, J) 
                  A (I, J) = A (IPVT, J) 
                  A (IPVT, J) = DUMMY 
   60          END DO 
            ENDIF 
!                                                                       
!  Perform the elimination step.                                        
!                                                                       
            DO 70 J = I + 1, N 
               A (J, I) = A (J, I) / A (I, I) 
               FAK = A (J, I) 
               DO 80 K = I + 1, N 
                  A (J, K) = A (J, K) - FAK * A (I, K) 
   80          END DO 
   70       END DO 
         ENDIF 
   40 END DO 
      IF (ABS (A (N, N) ) .LT.RELERR) MARK = 0 
      RETURN 
      END SUBROUTINE GAUSSP                         
!                                                                       
!                                                                       
      SUBROUTINE GAUSSS (N, A, LDA, IPIVOT, Y, X) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Calculating the solution X of a linear system of equations    *      
!  A * X = Y, where A has been factored via Gauss-elimination    *      
!  in SUBROUTINE GAUSSP.                                         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N        : order of the system of equations.                  *      
!  A        : 2-dimensional array A(1:LDA,1:N) containing the    *      
!             factors L and  R  with  P * A(ORG) = L * R  for    *      
!             P a permutation matrix. This array is the output   *      
!             matrix of SUBROUTINE GAUSSP.                       *      
!  LDA      : leading dimension of A as defined in the calling   *      
!             program.                                           *      
!  IPIVOT   : N-vector IPIVOT(1:N); it indicates the row         *      
!             interchanges in P * A relative to A(ORG). It is an *      
!             output of SUBROUTINE GAUSSP.                       *      
!  Y        : N-vector Y(1:N); the right hand side of the system *      
!             of equations.                                      *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  X        : N-vector X(1:N); the solution vector for the       *      
!             system of equations.                               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  authors   : Gisela Engeln-Muellges, Guido Dubois              *      
!  date      : 04.25.88                                          *      
!  source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z)       
      INTEGER J, K, I, N, IPVT, LDA
      DOUBLEPRECISION SUM, A, X, Y
      INTEGER IPIVOT (1:N)
      DIMENSION A (1:LDA, 1:N), Y (1:N), X (1:N) 
      IF (N.EQ.1) THEN 
         X (1) = Y (1) / A (1, 1) 
         RETURN 
      ENDIF 
!                                                                       
!  Updating the right hand side.                                        
!                                                                       
      IPVT = IPIVOT (1) 
      X (1) = Y (IPVT) 
      DO 10 I = 2, N 
         SUM = 0.0D0 
         DO 20 J = 1, I - 1 
            SUM = SUM + A (I, J) * X (J) 
   20    END DO 
         IPVT = IPIVOT (I) 
         X (I) = Y (IPVT) - SUM 
   10 END DO 
!                                                                       
!  Compute the solution vector X by backsubstitution.                   
!                                                                       
      X (N) = X (N) / A (N, N) 
      DO 50 I = N - 1, 1, - 1 
         SUM = 0.0D0 
         DO 40 K = N, I + 1, - 1 
            SUM = SUM + A (I, K) * X (K) 
   40    END DO 
         X (I) = (X (I) - SUM) / A (I, I) 
   50 END DO 
      RETURN 
      END SUBROUTINE GAUSSS                         

    subroutine test_linear_algebra()
    use mod_show_matrix
    implicit none
    integer i, j, s, k
    doubleprecision e, d, t1, t2, r
    doubleprecision, allocatable :: A(:,:), X(:), Y(:), EX(:)
    doubleprecision, allocatable :: B(:,:), P(:,:), EP(:,:)
    integer, parameter :: order = 3
    
        print *, "TESTING mod_linear_algebra"
        print *, " - Solving linear systems of equations using gaussian elimination."
        
        print *, ''
        print *, 'VECTOR-MATRIX'
        print '(1x,a12,1x,a12,1x,a12,1x,a12,1x,a12)', 'SIZE', 'REPEATS', 'RESIDUAL', 'TIME', 'RATE (MFP)'
        
        do i=1, 8
            s = 2**(i+1)
            k = 1024000/(s*s)
            
            if( allocated(A) )   deallocate(A )
            if( allocated(X) )   deallocate(X )
            if( allocated(Y) )   deallocate(Y )
            if( allocated(EX) )  deallocate(EX)          
            
            allocate(A(s,s))
            allocate(Y(s))
            
            
            call RANDOM_NUMBER(A)
            call RANDOM_NUMBER(Y)
            
            call CPU_TIME(t1)
            do j = 1, k
                X = gauss(A, Y)
            end do
            call CPU_TIME(t2)
            
            d = dble( t2 - t1 )
            EX = Y - matmul(A, X)    
            e = maxval( abs( ex ) )
            r = (s**order)*k/(d*1e6)
                        
            print '(1x,i12,1x,i12,1x,g12.4,1x,g12.4,1x,g12.4)', s, k, e, d, r
        end do
    
        print *, ''
        print *, 'MATRIX-MATRIX'
        print '(1x,a12,1x,a12,1x,a12,1x,a12,1x,a12)', 'SIZE', 'REPEATS', 'RESIDUAL', 'TIME', 'RATE (MFP)'
        
        do i=1, 8
            s = 2**(i+1)
            k = 1024000/(s*s)
            
            if( allocated(A) )   deallocate(A )
            if( allocated(B) )   deallocate(B )
            if( allocated(P) )   deallocate(P )
            if( allocated(EP) )  deallocate(EP)
            
            allocate(A(s,s))
            allocate(P(s,s))
                        
            call RANDOM_NUMBER(A)            
            call RANDOM_NUMBER(P)
            
            call CPU_TIME(t1)
            do j = 1, k
                B = gauss(A, P)
            end do
            call CPU_TIME(t2)
            
            d = dble( t2 - t1 )
            EP = P - matmul(A, B)
            e = maxval( abs( ep ) )
            r = (s**order)*k/(d*1e6)
            
            print '(1x,i12,1x,i12,1x,g12.4,1x,g12.4,1x,g12.4)', s, k, e, d, r
        end do
      
        end subroutine
    
    end module