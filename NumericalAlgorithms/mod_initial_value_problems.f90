    module mod_initial_value_problems
    use mod_linear_algebra
    
!                                                                       
!** Using the machine constant FMACHP, we determine EPS1 in order       
!** to avoid excessively small final steps near XE, EPS2 to check       
!** for zero.
!** EPS1 bounds the admissable initial step size HK below.              
!** EPS2 is used to determine a level for negligable quantities.        
!                                                                       
    DOUBLEPRECISION, PARAMETER :: EPS1 = (2*FMACHP)**0.75D0 
    DOUBLEPRECISION, PARAMETER :: EPS2 = 100.0D0 * (2* FMACHP) 
!                                                                       
!** HS as the optimal step size for approximating the      
!** Jacobi matrix. (This is done only once at the start.)               
!                                                                       
    DOUBLEPRECISION, PARAMETER :: HS = 10.0D0 * SQRT (2*FMACHP) 
    DOUBLEPRECISION, PARAMETER :: TOL = 26.0D0 * FMACHP
    
    INTEGER, PARAMETER :: NDGL = 20
    
      INTERFACE
        SUBROUTINE DERIVATIVE(X,Y,YS)
        IMPORT
        DOUBLEPRECISION X, Y(:), YS(:)
        END SUBROUTINE
      END INTERFACE    
    
    contains
    
      SUBROUTINE IVP (XK, HK, YK, N, DES, XE, EPSABS, EPSREL,           &
      INDEX, NMAX, NUSED, IERR)                                                
!                                                                       
!*****************************************************************      
!                                                                *      
!  For an approximation YK of the solution Y at XK of a system   *      
!  of ordinary differential equations of 1st order               *      
!                    Y' = F(X,Y),                                *      
!  this program computes an approximation for the solution Y     *      
!  at XE.                                                        *      
!  We use step size control in such a way that the error of the  *      
!  computed approximation falls either absolutely or relatively  *      
!  within the given error bounds EPSABS or EPSREL.               *      
!  In case of the Prince-Dormand embedding formula we also check *      
!  for stiffness of the system.                                  *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
! XK     - initial value of the independent variable X           *      
! HK     - proposed step size for the next step                  *      
! YK     - vector YK(1:N); value of the solution of the          *      
!          differential equation at XK                           *      
! N      - number of differential equations ( 1 <= N <= 20 )     *      
! DES    - right-hand side of the differential equation which    *      
!          must be provided as a SUBROUTINE in the form:         *      
!             SUBROUTINE DES (X, Y, F, N)                        *      
!          (starting with:  DOUBLE PRECISION Y(N), F(N), etc.).  *      
!          Here F represents the value of the differential       *      
!          equation's right-hand side at (X,Y). (DES has to be   *      
!          defined as EXTERNAL in the calling program)           *      
! XE     - location where the solution is desired; XE may not    *      
!          be chosen smaller than XK.                            *      
! EPSABS - error bound for the absolute accuracy of the desired  *      
!          solution. EPSABS has to be >= 0; if EPSABS = 0, only  *      
!          the relative accuracy is considered.                  *      
! EPSREL - error bound for the relative accuracy of the desired  *      
!          solution. EPSREL has to be >= 0; if EPSREL = 0,       *      
!          only the absolute accuracy is considered.             *      
! INDEX  - chooses the embedding formula with step size control: *      
!             = 0: RUNGE-KUTTA method  2nd/3rd order             *      
!             =-1: Prince-Dormand embedding formula of 4th/5th   *      
!                  order (checking for stiffness of the system   *      
!                  of the DEs here, see output IERR)             *      
!            else: formula of ENGLAND  4th/5th order             *      
! NMAX   - upper limit for the number of function evaluations    *      
!          allowed for the right-hand side F.                    *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
! XK     - location that was reached during last integration.    *      
!          If IERR = 0 normally XK is equal to XE                *      
! HK     - local step size used last ( it should remain un-      *      
!          changed for the next step )                           *      
! YK     - approximate value of the solution Y at the new        *      
!          location XK                                           *      
! NUSED  - number of function evaluations actually used          *      
! IERR   - error parameter:                                      *      
!          = 0: everything o.k.                                  *      
!               If Prince-Dormand was specified, the sytsem was  *      
!               not found to be stiff.                           *      
!          = 1: both error bounds  EPS...  are too small         *      
!               (relative to the machine constant)               *      
!          = 2: XE <= XK       (wrt machine constant)            *      
!          = 3: step size  HK <= 0   (wrt machine constant)      *      
!          = 4: N > 20  or  N <= 0                               *      
!          = 5: NUSED > NMAX:  the number of allowed functional  *      
!               evaluations is insufficient to determine an      *      
!               adequate approximate solution with the required  *      
!               accuracy; on termination, XK and HK contain      *      
!               the current values                               *      
!          =-1: The computations terminated ok, but the Prince-  *      
!               Dormand formula has detected possible stiffness. *      
!          =-2: The computations have terminated ok, but the     *      
!               Prince-Dormand formula has recognized the system *      
!               as stiff using two criteria: we recommend to use *      
!               a method suited for stiff DEs instead.           *      
!          =-3: NUSED > NMAX:  The number of allowed functional  *      
!               evaluations does not suffice, see IERR = 5.      *      
!               Moreover, the Prince-Dormand formula seems to    *      
!               indicate that the system is stiff; we recommend  *      
!               to use a suitable stiff DE solver.               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: RUKU23, ENGL45, PRDO45                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Klaus Niederdrenk                               *      
!  Date        : 11.01.1985 / 4.1.1995                           *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      
      !INTERFACE
      !  SUBROUTINE DES(X,Y,N,YS)
      !  IMPORT
      !  INTEGER N
      !  DOUBLEPRECISION X, Y(N), YS(N)
      !  END SUBROUTINE
      !END INTERFACE
      
      PROCEDURE (DERIVATIVE) :: DES
      
      INTEGER N, NUSED
      DOUBLEPRECISION YK (N), HK, XE, DIFF, S, XK, EPSABS, EPSREL
      DOUBLEPRECISION Y (20), YT (20), SG, XEND, YMAX, DUMMY
      LOGICAL IEND, STEIF1, STEIF2, STEIFA 
      INTEGER IANZ, IERR, INDEX, I, NMAX
!                                                                       
!** Preassign local variables                                           
!                                                                       
      SG = SIGN (1.0D0, XE) 
      XEND = (1.0D0 - SG * EPS2) * XE 
      IERR = 0 
      NUSED = 0 
      IEND = .FALSE. 
      STEIF1 = .FALSE. 
      STEIF2 = .FALSE. 
      STEIFA = .FALSE. 
      IANZ = 0 
!                                                                       
!** Check input parameters                                              
!                                                                       
      YMAX = MAXVAL(ABS(YK)) !DVNORM (YK, Y00, N) 
      IF (EPSABS.LE.EPS2 * YMAX.AND.EPSREL.LE.EPS2) THEN 
         IERR = 1 
      ELSEIF (XEND.LT.XK) THEN 
         IERR = 2 
      ELSEIF (HK.LT.EPS2 * ABS (XK) ) THEN 
         IERR = 3 
      ELSEIF (N.LE.0.OR.N.GT.20) THEN 
         IERR = 4 
      ENDIF 
      IF (IERR.NE.0) RETURN 
!                                                                       
!**********  C O N T R O L L I N G   a l g o r i t h  **********        
!                                                                       
      IF (XK + HK.GT.XEND) THEN 
         HK = XE-XK 
         DUMMY = HK 
         IEND = .TRUE. 
      ENDIF 
!                                                                       
!** Integrate on the interval *[*XK, XE*]* in suitable steps            
!                                                                       
   50 CONTINUE 
!                                                                       
!** Call the desired one step method                                    
!                                                                       
      IF (INDEX.EQ.0) THEN 
         CALL RUKU23 (XK, HK, YK, N, DES, Y, YT) 
         NUSED = NUSED+3 
      ELSEIF (INDEX.EQ. - 1) THEN 
         CALL PRDO45 (XK, HK, YK, N, DES, Y, YT, STEIF1, STEIFA) 
         NUSED = NUSED+7 
         IF (STEIFA) THEN 
            IANZ = IANZ + 1 
            IF (IANZ.GE.3) STEIF2 = .TRUE. 
         ELSE 
            IANZ = 0 
         ENDIF 
      ELSE 
         CALL ENGL45 (XK, HK, YK, N, DES, Y, YT) 
         NUSED = NUSED+6 
      ENDIF 
!                                                                       
      DIFF = MAXVAL(ABS(Y-YT)) !DVNORM (Y, YT, N) 
!                                                                       
      IF (DIFF.LT.EPS2) THEN 
         S = 2.0D0 
      ELSE 
         YMAX = MAXVAL(ABS(YT)) ! DVNORM (YT, Y00, N) 
         S = SQRT (HK * (EPSABS + EPSREL * YMAX) / DIFF) 
         IF (INDEX.NE.0) S = SQRT (S) 
      ENDIF 
!                                                                       
      IF (S.GT.1.0D0) THEN 
!                                                                       
!** accept the performed integration with step size HK                  
!                                                                       
         DO 60 I = 1, N 
            YK (I) = YT (I) 
   60    END DO 
!                                                                       
         XK = XK + HK 
!                                                                       
!** if then endpoint XE has been reached or if more than the            
!** allowable function evaluations were used: go back                   
!                                                                       
   70    IF (IEND) THEN 
            HK = DUMMY 
            IF (INDEX.EQ. - 1) THEN 
               IF (STEIF1.OR.STEIF2) IERR = - 1 
               IF (STEIF1.AND.STEIF2) IERR = - 2 
            ENDIF 
            RETURN 
         ELSEIF (NUSED.GT.NMAX) THEN 
            IERR = 5 
            IF (INDEX.EQ. - 1.AND. (STEIF1.OR.STEIF2) ) IERR = - 3 
            RETURN 
         ENDIF 
!                                                                       
!** increase the step size for the next step maximally by a factor of tw
!                                                                       
         HK = HK * MIN (2.0D0, 0.98D0 * S) 
!                                                                       
         IF ( (XK + HK) .GE.XEND) THEN 
            DUMMY = HK 
            HK = XE-XK 
            IEND = .TRUE. 
!                                                                       
!** if very close to XE: go back                                        
!                                                                       
            IF (HK.LT.EPS1 * ABS (XE) ) GOTO 70 
         ENDIF 
      ELSE 
!                                                                       
!** the previous step is unaccaptable, the step size HK must be decrease
!** at most it must be halved                                           
!                                                                       
         HK = HK * MAX (0.5D0, 0.98D0 * S) 
         IEND = .FALSE. 
      ENDIF 
!                                                                       
      GOTO 50 
!                                                                       
      END SUBROUTINE IVP                            
!                                                                       
!                                                                       
      SUBROUTINE RUKU23 (X, H, Y, N, DES, Y2, Y3) 
!                                                                       
!*****************************************************************      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Starting with an approximation Y at X, this program computes  *      
!  approximations Y2 and Y3 at X + H using the RUNGE-KUTTA       *      
!  embedding formula of 2nd and 3rd order for a system of        *      
!  ordinary differential equations of 1st order                  *      
!                  Y' = F(X,Y).                                  *      
!  The system of N ordinary differential equations of 1st order  *      
!  must be provided by a SUBROUTINE DES.                         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  X     - initial value for the independent variable X          *      
!  H     - step size                                             *      
!  Y     - vector Y(1:N); value for the solution of the          *      
!          differential equation at X                            *      
!  N     - number of differential equations  ( 1 <= N <= 20 )    *      
!  DES   - right hand side of the differential equation. It has  *      
!          provided by the user as a SUBROUTINE in the form:     *      
!             SUBROUTINE DES (X, Y, N, F)                        *      
!          (starting with: DOUBLE PRECISION Y(N), F(N), etc. ).  *      
!          Here F denotes the right hand side of the differential*      
!          equation at (X,Y). ( In the calling program DES has to*      
!          declared as EXTERNAL)                                 *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  Y2    - vector Y2(1:N); approximate solution using the 2nd    *      
!          order method for solving the differential equation    *      
!          at X+H                                                *      
!  Y3    - vector Y3(1:N); approximate solution using the 3rd    *      
!          order method for solving the differential equation    *      
!          at X+H                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : none                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Klaus Niederdrenk                                  *      
!  date     : 11.01.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      
      !INTERFACE
      !  SUBROUTINE DES(X,Y,N,YS)
      !  IMPORT
      !  INTEGER N
      !  DOUBLEPRECISION X, Y(N), YS(N)
      !  END SUBROUTINE
      !END INTERFACE
      
      PROCEDURE (DERIVATIVE) :: DES
      
      INTEGER N, I
      DOUBLEPRECISION X, H, Y (N), Y2 (N), Y3 (N), DUMMY (20) 
      DOUBLEPRECISION K1 (20), K2 (20), K3 (20) 
!                                                                       
      CALL DES (X, Y, K1) 
      DO 10 I = 1, N 
         DUMMY (I) = Y (I) + H * K1 (I) 
   10 END DO 
      CALL DES (X + H, DUMMY, K2) 
      DO 20 I = 1, N 
         DUMMY (I) = Y (I) + 0.25D0 * H * (K1 (I) + K2 (I) ) 
   20 END DO 
      CALL DES (X + 0.5D0 * H, DUMMY, K3) 
!                                                                       
      DO 100 I = 1, N 
         Y2 (I) = Y (I) + 0.5D0 * H * (K1 (I) + K2 (I) ) 
         Y3 (I) = Y (I) + H / 6.0D0 * (K1 (I) + K2 (I) + 4.0D0 * K3 (I) &
         )                                                              
  100 END DO 
      RETURN 
      END SUBROUTINE RUKU23                         
!                                                                       
!                                                                       
      SUBROUTINE ENGL45 (X, H, Y, N, DES, Y4, Y5) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Starting from an approximation Y at X, this program uses the  *      
!  ENGLAND embedding formula of 4th and 5th order to find appro- *      
!  ximations Y4 and Y5 at X + H for the solution of the system of*      
!  differential equations of 1st order                           *      
!                  Y' = F(X,Y) .                                 *      
!  The system contains N ordinary differential equations of 1st  *      
!  order, with the right hand side F(X,Y) provided by the user   *      
!  in a SUBROUTINE DES.                                          *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  X     - initial value for the independent variable X          *      
!  H     - step size                                             *      
!  Y     - vector Y(1:N); value for the solution of the          *      
!          differential equation at X                            *      
!  N     - number of differential equations  ( 1 <= N <= 20 )    *      
!  DES   - right hand side of the differential equation. It has  *      
!          provided by the user as a SUBROUTINE in the form:     *      
!             SUBROUTINE DES (X, Y, N, F)                        *      
!          (starting with: DOUBLE PRECISION Y(N), F(N), etc. ).  *      
!          Here F denotes the right hand side of the differential*      
!          equation at (X,Y). ( In the calling program DES has to*      
!          declared as EXTERNAL)                                 *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  Y4    - vector Y4(1:N); approximate solution using the 4th    *      
!          order method for solving the differential equation    *      
!          at X+H                                                *      
!  Y5    - vector Y5(1:N); approximate solution using the 5th    *      
!          order method for solving the differential equation    *      
!          at X+H                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : none                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Klaus Niederdrenk                                  *      
!  date     : 11.01.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      
      !INTERFACE
      !  SUBROUTINE DES(X,Y,N,YS)
      !  IMPORT
      !  INTEGER N
      !  DOUBLEPRECISION X, Y(N), YS(N)
      !  END SUBROUTINE
      !END INTERFACE
      
      PROCEDURE (DERIVATIVE) :: DES

      INTEGER N, I
      DOUBLEPRECISION X, H, Y (N), Y4 (N), Y5 (N) 
      DOUBLEPRECISION DUMMY (20) 
      DOUBLEPRECISION K1 (20), K2 (20), K3 (20), K4 (20), K5 (20),      &
      K6 (20)                                                           
!                                                                       
      CALL DES (X, Y, K1) 
      DO 10 I = 1, N 
         DUMMY (I) = Y (I) + 0.5D0 * H * K1 (I) 
   10 END DO 
      CALL DES (X + 0.5D0 * H, DUMMY, K2) 
!                                                                       
      DO 20 I = 1, N 
         DUMMY (I) = Y (I) + 0.25D0 * H * (K1 (I) + K2 (I) ) 
   20 END DO 
      CALL DES (X + 0.5D0 * H, DUMMY, K3) 
!                                                                       
      DO 30 I = 1, N 
         DUMMY (I) = Y (I) + H * ( - K2 (I) + 2.0D0 * K3 (I) ) 
   30 END DO 
      CALL DES (X + H, DUMMY, K4) 
!                                                                       
      DO 40 I = 1, N 
         DUMMY (I) = Y (I) + H / 27.0D0 * (7.0D0 * K1 (I) + 10.0D0 * K2 &
         (I) + K4 (I) )                                                 
   40 END DO 
      CALL DES (X + 2.0D0 / 3.0D0 * H, DUMMY, K5) 
!                                                                       
      DO 50 I = 1, N 
         DUMMY (I) = Y (I) + 0.16D-02 * H * (28.0D0 * K1 (I) - 125.0D0 *&
         K2 (I) + 546.0D0 * K3 (I) + 54.0D0 * K4 (I) - 378.0D0 * K5 (I) &
         )                                                              
   50 END DO 
      CALL DES (X + 0.2D0 * H, DUMMY, K6) 
!                                                                       
      DO 100 I = 1, N 
         Y4 (I) = Y (I) + H / 6.0D0 * (K1 (I) + 4.0D0 * K3 (I) + K4 (I) &
         )                                                              
         Y5 (I) = Y (I) + H / 336.0D0 * (14.0D0 * K1 (I) + 35.0D0 * K4 (&
         I) + 162.0D0 * K5 (I) + 125.0D0 * K6 (I) )                     
  100 END DO 
      RETURN 
      END SUBROUTINE ENGL45                         
!                                                                       
!                                                                       
      SUBROUTINE PRDO45 (X, H, Y, N, DES, Y4, Y5, ST1, ST2) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Starting from an approximation Y at X, this program uses the  *      
!  Prince-Dormand embedding formula of 4th and 5th order to find *      
!  apprximations Y4 and Y5 at X + H for the solution of the      *      
!  system of differential equations of 1st order                 *      
!                  Y' = F(X,Y) .                                 *      
!  The system contains N ordinary differential equations of 1st  *      
!  order, with the right hand side F(X,Y) provided by the user   *      
!  in a SUBROUTINE DES.                                          *      
!  This program tests for stiffness in two ways.                 *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  X     - initial value for the independent variable X          *      
!  H     - step size                                             *      
!  Y     - vector Y(1:N); value for the solution of the          *      
!          differential equation at X                            *      
!  N     - number of differential equations  ( 1 <= N <= 20 )    *      
!  DES   - right hand side of the differential equation. It has  *      
!          provided by the user as a SUBROUTINE in the form:     *      
!             SUBROUTINE DES (X, Y, N, F)                        *      
!          (starting with: DOUBLE PRECISION Y(N), F(N), etc. ).  *      
!          Here F denotes the right hand side of the differential*      
!          equation at (X,Y). ( In the calling program DES has to*      
!          be declared as EXTERNAL)                              *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  Y4    - vector Y4(1:N); approximate solution using the 4th    *      
!          order method for solving the differential equation    *      
!          at X+H                                                *      
!  Y5    - vector Y5(1:N); approximate solution using the 5th    *      
!          order method for solving the differential equation    *      
!          at X+H                                                *      
!  ST1   - logical variable: .TRUE. , if the test for stiffness  *      
!          (using the approximate dominat eigenvalue) was        *      
!          positive ; else the value of ST! remains unaltered.   *      
!  ST2   - logical variable: .TRUE. , if the test of stiffness   *      
!          due to  Shampine  and  Hiebert  is positive; else the *      
!          value of ST2 is set to  .FALSE. .                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required :                                        *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Klaus Niederdrenk                                  *      
!  date     : 1.4.1995                                           *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      !INTERFACE
      !  SUBROUTINE DES(X,Y,N,YS)
      !  IMPORT
      !  INTEGER N
      !  DOUBLEPRECISION X, Y(N), YS(N)
      !  END SUBROUTINE
      !END INTERFACE
      
      PROCEDURE (DERIVATIVE) :: DES

      INTEGER N, I
      DOUBLEPRECISION X, H, Y (N), Y4 (N), Y5 (N) 
      DOUBLEPRECISION DUMMY (20) 
      DOUBLEPRECISION K1 (20), K2 (20), K3 (20), K4 (20), K5 (20),      &
      K6 (20), K7 (20), G6 (20), G7 (20)                                
      LOGICAL ST1, ST2 
!                                                                       
      CALL DES (X, Y, K1) 
      DO 10 I = 1, N 
         DUMMY (I) = Y (I) + 0.2D0 * H * K1 (I) 
   10 END DO 
      CALL DES (X + 0.2D0 * H, DUMMY, K2) 
!                                                                       
      DO 20 I = 1, N 
         DUMMY (I) = Y (I) + 0.075D0 * H * (K1 (I) + 3.0D0 * K2 (I) ) 
   20 END DO 
      CALL DES (X + 0.3D0 * H, DUMMY, K3) 
!                                                                       
      DO 30 I = 1, N 
         DUMMY (I) = Y (I) + H / 45.0D0 * (44.0D0 * K1 (I) - 168.0D0 *  &
         K2 (I) + 160.0D0 * K3 (I) )                                    
   30 END DO 
      CALL DES (X + 0.8D0 * H, DUMMY, K4) 
!                                                                       
      DO 40 I = 1, N 
         DUMMY (I) = Y (I) + H / 6561.0D0 * (19372.0D0 * K1 (I) -       &
         76080.0D0 * K2 (I) + 64448.0D0 * K3 (I) - 1908.0D0 * K4 (I) )  
   40 END DO 
      CALL DES (X + 8.0D0 / 9.0D0 * H, DUMMY, K5) 
!                                                                       
      DO 50 I = 1, N 
         G6 (I) = Y (I) + H / 167904.0D0 * (477901.0D0 * K1 (I) -       &
         1806240.0D0 * K2 (I) + 1495424.0D0 * K3 (I) + 46746.0D0 * K4 ( &
         I) - 45927.0D0 * K5 (I) )                                      
   50 END DO 
      CALL DES (X + H, G6, K6) 
!                                                                       
      DO 60 I = 1, N 
         G7 (I) = Y (I) + H / 142464.0D0 * (12985.0D0 * K1 (I) +        &
         64000.0D0 * K3 (I) + 92750.0D0 * K4 (I) - 45927.0D0 * K5 (I)   &
         + 18656.0D0 * K6 (I) )                                         
   60 END DO 
      CALL DES (X + H, G7, K7) 
!                                                                       
      DO 100 I = 1, N 
         Y5 (I) = G7 (I) 
         Y4 (I) = Y (I) + H / 21369600.0D0 * (1921409.0D0 * K1 (I)      &
         + 9690880.0D0 * K3 (I) + 13122270.0D0 * K4 (I) - 5802111.0D0 * &
         K5 (I) + 1902912.0D0 * K6 (I) + 534240.0D0 * K7 (I) )          
100 END DO 
!                                                                       
!**  Test for stiffness via dominant eigenvalue (approximately)         
!                                                                       
      !IF (DVNORM (K7, K6, N) .GT. 3.3D0 * DVNORM (G7, G6, N) ) THEN
      IF( MAXVAL(ABS(K7-K6)) > 3.30D0 * MAXVAL(ABS(G7-G6)) ) THEN
          ST1 = .TRUE.
      END IF
!                                                                       
!**  One step of testing stiffness according to Shampine and Hiebert    
!                                                                       
      DO 110 I = 1, N 
         G6 (I) = H * (2.2D0 * K2 (I) + 0.13D0 * K4 (I) + 0.144D0 * K5 (&
         I) )                                                           
         G7 (I) = H * (2.134D0 * K1 (I) + 0.24D0 * K3 (I) + 0.1D0 * K6 (&
         I) )                                                           
  110 END DO 
      !IF (DVNORM (G6, G7, N) .LT.DVNORM (Y4, Y5, N) ) THEN 
      IF (MAXVAL(ABS(G6 - G7)) <  MAXVAL(ABS(Y4 - Y5)) ) THEN 
         ST2 = .TRUE. 
      ELSE 
         ST2 = .FALSE. 
      ENDIF 
      RETURN 
    END SUBROUTINE PRDO45                         

    
    SUBROUTINE FRKFSY (A, DA, N, Y, DES, H, HMX, ABSERR, RELERR, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  A system of ordinary differential equations of 1st order is   *      
!  integrated by applying the RUNGE-KUTTA-FEHLBERG method        *      
!  [ of order O(H**5) ] with estimates for the local error and   *      
!  step size control.                                            *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  A      : starting value for the integration interval          *      
!  DA     : length of the integration interval;                  *      
!           DA may be < 0.0 if we want to integrate to the left. *      
!  N      : number of equations; N < 11.                         *      
!  Y      : vector Y(1:N); the initial values at A               *      
!  DES    : SUBROUTINE, that describes the system of differential*      
!           equations, given in the following form:              *      
!                      SUBROUTINE  DES(X,Y,YS)                   *      
!                      X : independent variable                  *      
!                      Y : vector of dependent variables         *      
!                      YS: vector YS(I)=DY(I)/DX of derivatives  *      
!                          at X, I=1,...,N                       *      
!                      Y and YS are dimensioned as DOUBLE        *      
!                      PRECISION Y(1), YS(1), however, they may  *      
!                      be to be used as a vector of length N.    *      
!           example :  SUBROUTINE DES(X,Y,YS)                    *      
!                      DOUBLE PRECISION Y(1),YS(1)               *      
!                      YS(1)=Y(1)                                *      
!                      YS(2)=-Y(2)                               *      
!                      RETURN                                    *      
!                      END                                       *      
!  H      : initial step size; if H is provided unrealistically, *      
!           H is modified internally; H may be negative if       *      
!           DA < 0.0.                                            *      
!  HMX    : upper bound for the step size magnitude used during  *      
!           calculation. HMX > 0.0                               *      
!  ABSERR :] bounds for the acceptable local error, relative to  *      
!  RELERR :] the current step size. If the following holds for   *      
!         :] each component of the computed solution Y(I)        *      
!                  ABS ( estimate of the local error) .LE.       *      
!                      ABS(H)*(RELERR*ABS(Y(I))+ABSERR),         *      
!           then the solution is accepted in the current step.   *      
!           If ABSERR = 0.0, we test for the relative error;     *      
!           If RELERR = 0.0, we test for the absolute error.     *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  A      : last x value for which a solution was successfully   *      
!           determined. Normally the following will hold:        *      
!           A on output = A on input + DA.                       *      
!  Y      : computed solution vector at  A on output             *      
!  H      : optimal step size, which was used for the last step. *      
!  IERR   : = 1, everything o.k.; solution found at A + DA.      *      
!           = 2, after 3000 calls of SUBROUTINE DES we stop with-*      
!                out having reached the endpoint A+DA. If com-   *      
!                putations are to be continued, call FRKFSY again*      
!                with unchanged parameters.                      *      
!           = 3, false input data; i.e.                          *      
!                ABSERR.LT.0.0     or    RELERR.LT.0.0     or    *      
!                ABSERR + RELERR = 0.0  or  HMX.LE.0.0: Return.  *      
!           = 4, the optimal step size cannot be achieved for the*      
!                computer.     RETURN                            *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : MACHPD                                 *      
!                                                                *      
!                                                                *      
!  sources : SHAMPINE/ALLEN, see [SHAM73].                       *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Richard Reuter                                     *      
!  date     : 02.09.1983                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      
      !INTERFACE
      !  SUBROUTINE DES(X,Y,YS)
      !  IMPORT
      !  DOUBLEPRECISION X, Y(:), YS(:)
      !  END SUBROUTINE
      !END INTERFACE
      
      PROCEDURE (DERIVATIVE) :: DES
    
      INTEGER N, IERR, I, LFD, IAD
      DOUBLEPRECISION A, DA, Y (N), H, HMX, ABSERR, RELERR, B, HMAX, HF, X
      DOUBLEPRECISION YT (10), T (10), R (10), K1 (10), K2 (10) 
      DOUBLEPRECISION K3 (10), K4 (10), K5 (10), K6 (10) 
      DOUBLEPRECISION QUOT, TR
      
!                                                                       
!     check the input data                                              
!                                                                       
      IERR = 3 
      IF (RELERR.LT.0.0D0.OR.ABSERR.LT.0.0D0.OR.RELERR +                &
      ABSERR.EQ.0.0D0.OR.HMX.LE.0.0D0) RETURN                           
!                                                                       
      IERR = 4 
      B = A + DA 
      IF (ABS (DA) .LE. TOL * MAX (ABS (A), ABS (B) ) ) &
      RETURN                                                            
!                                                                       
      HMAX = MIN (HMX, ABS (DA) ) 
      IF (ABS (H) .LE. TOL * ABS (A) ) H = HMAX 
!                                                                       
!     Initialize counter for calls of SUBROUTINE DES                    
!                                                                       
      LFD = 0 
      IAD = 0 
!                                                                       
!     H is bounded by HMAX and is chosen so that the                    
!     endpoint B is reached, if possible.                               
!                                                                       
    3 H = DSIGN (MIN (ABS (H), HMAX), DA) 
      IF (ABS (B - A) .LE.1.25D0 * ABS (H) ) THEN 
         HF = H 
!                                                                       
!        if IAD=1 and H=B-A acceptable, we stop after                   
!        the next integration step.                                     
!                                                                       
         IAD = 1 
         H = B - A 
      ENDIF 
!                                                                       
!     an integration step is executed                                   
!                                                                       
      CALL DES (A, Y, K1) 
      LFD = LFD+1 
    5 CONTINUE 
      X = 0.25D0 * H 
      DO 6 I = 1, N 
         YT (I) = Y (I) + X * K1 (I) 
    6 END DO 
      X = A + X 
      CALL DES (X, YT, K2) 
      DO 7 I = 1, N 
         YT (I) = Y (I) + H * (K1 (I) * (3.0D0 / 32.0D0) + K2 (I)       &
         * (9.0D0 / 32.0D0) )                                           
    7 END DO 
      X = A + H * (3.0D0 / 8.0D0) 
      CALL DES (X, YT, K3) 
      DO 8 I = 1, N 
         YT (I) = Y (I) + H * (K1 (I) * (1932.0D0 / 2197.0D0) - K2 (I)  &
         * (7200.0D0 / 2197.0D0) + K3 (I) * (7296.0D0 / 2197.0D0) )     
    8 END DO 
      X = A + H * (12.0D0 / 13.0D0) 
      CALL DES (X, YT, K4) 
      DO 9 I = 1, N 
         YT (I) = Y (I) + H * (K1 (I) * (439.0D0 / 216.0D0) - 8.0D0 *   &
         K2 (I) + K3 (I) * (3680.0D0 / 513.0D0) - K4 (I) * (845.0D0 /   &
         4104.0D0) )                                                    
    9 END DO 
      X = A + H 
      CALL DES (X, YT, K5) 
      DO 10 I = 1, N 
         YT (I) = Y (I) + H * ( - K1 (I) * (8.0D0 / 27.0D0) + 2.0D0 *   &
         K2 (I) - K3 (I) * (3544.0D0 / 2565.0D0) + K4 (I) * (1859.0D0 / &
         4104.0D0) - K5 (I) * (11.0D0 / 40.0D0) )                       
   10 END DO 
      X = A + 0.5D0 * H 
      CALL DES (X, YT, K6) 
      DO 11 I = 1, N 
         T (I) = K1 (I) * (25.0D0 / 216.0D0) + K3 (I) * (1408.0D0 /     &
         2565.0D0) + K4 (I) * (2197.0D0 / 4104.0D0) - K5 (I) * 0.20D0   
         YT (I) = Y (I) + H * T (I) 
   11 END DO 
!                                                                       
!     YT(I) now represents the latest result of this pass.              
!     Determine R(I), the estimate of the local                         
!     error, relative to the current step size.                         
!                                                                       
      DO 12 I = 1, N 
         R (I) = K1 (I) / 360.0D0 - K3 (I) * (128.0D0 / 4275.0D0)       &
         - K4 (I) * (2197.0D0 / 75240.0D0) + K5 (I) / 50.0D0 + K6 (I)   &
         * (2.0D0 / 55.0D0)                                             
   12 END DO 
!                                                                       
!     Check accuracy                                                    
!                                                                       
      QUOT = 0.0D0 
      DO 13 I = 1, N 
         TR = ABS (R (I) ) / (RELERR * ABS (YT (I) ) + ABSERR) 
         QUOT = MAX (QUOT, TR) 
   13 END DO 
!                                                                       
!     If  QOUOT.LE.1.0   ==> integration step is accepted               
!                                                                       
      IF (QUOT.LE.1.0D0) THEN 
!                                                                       
!        result is accepted                                             
!                                                                       
         DO 14 I = 1, N 
            Y (I) = YT (I) 
   14    END DO 
         A = A + H 
!                                                                       
!        if A=B ,  RETURN                                               
!                                                                       
         IF (IAD.EQ.1) THEN 
            IERR = 1 
            H = HF 
            RETURN 
         ENDIF 
!                                                                       
!        prepare next step                                              
!                                                                       
         QUOT = MAX (QUOT, 6.5536D-4) 
      ENDIF 
      QUOT = MIN (QUOT, 4096.0D0) 
      H = 0.8D0 * H / SQRT (SQRT (QUOT) ) 
!                                                                       
!     We just achieved that H was increased by at most a factor of 5,   
!     or alternatively, that it was decreased by a factor of 10         
!     at most                                                           
!                                                                       
      IF (ABS (H) .LE. TOL * ABS (A) ) THEN 
         IERR = 4 
         RETURN 
      ENDIF 
      LFD = LFD+5 
      IF (LFD.GE.2995) THEN 
         IERR = 2 
         RETURN 
      ENDIF 
      IF (QUOT.LE.1.0D0) THEN 
!                                                                       
!        the step was successful. Continue with another step.           
!                                                                       
         GOTO 3 
      ELSE 
!                                                                       
!        the step is repeated for a smaller H.                          
!                                                                       
         IAD = 0 
         GOTO 5 
      ENDIF 
      END SUBROUTINE FRKFSY                         
    
    
      SUBROUTINE GEAR4 (XK, HK, YK, N, DES, XE, EPSABS, EPSREL,         &
      NMAX, NUSED, IERR)                                                      
!                                                                       
!*****************************************************************      
!                                                                *      
!  Starting from an approximation YK for the solution Y of a     *      
!  system of ordinary differential equations of first order      *      
!                    Y' = F(X,Y)                                 *      
!  at XK, this program computes an approximate solution YE at XE.*      
!  Here we compute internally with a step size control, that     *      
!  ensures the error of the computed solution to be less than the*      
!  given absolute or relative error bounds EPSABS and EPSREL.    *      
!  These bounds must be specified small enough for good results. *      
!  The method used is the multistep method of Gear of fourth     *      
!  order which is highly capable of solving stiff DEs. (Stiff DEs*      
!  are those DE systems which have solution components of very   *      
!  disparate growths.)                                           *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  XK    - starting value for X                                  *      
!  HK    - proposed step size for first step                     *      
!  YK    - vector YK(1:N); Y value of the solution to the DE     *      
!          at XK                                                 *      
!  N     - number of DEs  ( 1 <= N <= 20 )                       *      
!  DES   - right hand side of the DE, given as a subroutine:     *      
!             SUBROUTINE DES (X, Y, N, F)                        *      
!          ( starting with: DOUBLE PRECISION Y(N), F(N), X ).    *      
!          Here F is the value of the right hand side at (X,Y).  *      
!          ( DES must be declared as EXTERNAL in the calling     *      
!            program.)                                           *      
!  XE    - X value for desired solution; XE > XK.                *      
! EPSABS - error bound for absolute error; EPSABS >= 0; if       *      
!          EPSABS = 0 the algorithm maintains only the relative  *      
!          accuracy.                                             *      
! EPSREL - error bound for relative error; EPSREL >= 0; if       *      
!          EPSREL = 0 the algorithm maintains only the absolute  *      
!          accuracy.                                             *      
!  NMAX  - maximal number of evaluations of the right hand side  *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  XK    - final X value of the integration. If IERR = 0,        *      
!          usually XK = XE (within machine precision).           *      
!  HK    - terminal step size used; should be used for subsequent*      
!          integrations                                          *      
!  YK    - approximate value for the solution at XK              *      
!  NUSED - number of actual evaluations of the right hand side   *      
!  IERR  - error parameter:                                      *      
!          = 0: all is o.k.                                      *      
!          = 1: both error bounds  EPS...  too small             *      
!                              (relative to the machine constant)*      
!          = 2: XE <= XK       (relative to the machine constant)*      
!          = 3: step size  HK <= 0  (rel. to machine precision)  *      
!          = 4: N > 20   or   N <= 0                             *      
!          = 5: NUSED > NMAX: Number of allowed function         *      
!               evaluations was exceeded; try to restart with    *      
!               XK, YK and HK.                                   *      
!          = 6: The Jacobi matrix is singular; XK, YK, HK contain*      
!               the values reached.                              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines used:  IVP, GAUSSP, GAUSSS                        *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Klaus Niederdrenk                               *      
!  Date        : 1.22.1996                                       *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      
      !INTERFACE
      !  SUBROUTINE DES(X,Y,N,YS)
      !  IMPORT
      !  INTEGER N
      !  DOUBLEPRECISION X, Y(N), YS(N)
      !  END SUBROUTINE
      !END INTERFACE
      
      PROCEDURE (DERIVATIVE) :: DES
          
      INTEGER N, NMAX, NUSED, IERR, I, K, NANL, IPIVOT, ITER, MARK
      DOUBLEPRECISION XK, YK (1:N), SG, XE, HK, EPSABS, EPSREL, HK1
      DOUBLEPRECISION XEND, YMAX, DUMMY, HELP, CON, HALT, ZJP1, YKP1
      DOUBLEPRECISION XKA, XKE, HKA, F, ZJ, FS, FSG, D, DIFF, EPS
      DOUBLEPRECISION QUOT1, QUOT2, QUOT3, QUOT4, Q
      DIMENSION ZJ (0:4, NDGL), ZJP1 (0:4, NDGL), F(NDGL) 
      DIMENSION FS (NDGL, NDGL), HELP (NDGL)
      DIMENSION YKP1 (NDGL), CON(NDGL)
      DIMENSION D (NDGL), IPIVOT(NDGL), FSG(NDGL, NDGL) 
      LOGICAL IEND
!                                                                       
!** Initialize                                                          
!                                                                       
      SG = DSIGN (1.0D0, XE) 
      XEND = (1.0D0 - SG * EPS2) * XE 
      IERR = 0 
      NUSED = 0 
      IEND = .FALSE. 
!                                                                       
!** Check input parameters                                              
!                                                                       
      YMAX = MAXVAL(ABS(YK)) !DVNORM (YK, Y0, N) 
      IF (EPSABS.LE.EPS2 * YMAX.AND.EPSREL.LE.EPS2) THEN 
         IERR = 1 
      ELSEIF (XEND.LT.XK) THEN 
         IERR = 2 
      ELSEIF (HK.LT.EPS2 * ABS (XK) ) THEN 
         IERR = 3 
      ELSEIF (N.LE.0.OR.N.GT.NDGL) THEN 
         IERR = 4 
      ENDIF 
      IF (IERR.NE.0) RETURN 
!                                                                       
!****  compute first integration   ****                                 
!                                                                       
      IF (XK + HK.GT.XEND) THEN 
         HK = XE-XK 
         DUMMY = HK 
         IEND = .TRUE. 
      ENDIF 
      DO 20 I = 1, N 
         HELP (I) = YK (I) 
   20 END DO 
      XKA = XK 
      XKE = XKA 
      HKA = 0.25 * HK 
      HK1 = HKA 
      DO 40 K = 1, 4 
         XKE = XKE+HKA 
         CALL IVP (XKA, HK1, HELP, N, DES, XKE, EPSABS, EPSREL, 1, NMAX &
         - NUSED, NANL, IERR)                                           
         NUSED = NUSED+NANL 
         IF (IERR.NE.0) RETURN 
         DO 30 I = 1, N 
            ZJP1 (K, I) = HELP (I) 
   30    END DO 
   40 END DO 
      CALL DES (XK, YK, F) 
      NUSED = NUSED+1 
!                                                                       
!** Determine first Gear-Nordsieck approximation                        
!                                                                       
      DO 50 I = 1, N 
         ZJ (0, I) = YK (I) 
         ZJ (1, I) = HK * F (I) 
         ZJ (2, I) = 1.0D0 / 24.0D0 * (35.0D0 * YK (I) - 104.0D0 * ZJP1 &
         (1, I) + 114.0D0 * ZJP1 (2, I) - 56.0D0 * ZJP1 (3, I) + 11.0D0 &
         * ZJP1 (4, I) )                                                
         ZJ (3, I) = 1.0D0 / 12.0D0 * ( - 5.0D0 * YK (I) + 18.0D0 *     &
         ZJP1 (1, I) - 24.0D0 * ZJP1 (2, I) + 14.0D0 * ZJP1 (3, I)      &
         - 3.0D0 * ZJP1 (4, I) )                                        
         ZJ (4, I) = 1.0D0 / 24.0D0 * (YK (I) - 4.0D0 * ZJP1 (1, I)     &
         + 6.0D0 * ZJP1 (2, I) - 4.0D0 * ZJP1 (3, I) + ZJP1 (4, I) )    
   50 END DO 
!                                                                       
!                                                                       
!****  S t e p  S i z e  A l g o r i t h m   ****                       
!                                                                       
!                                                                       
   75 CONTINUE 
!                                                                       
!** Compute implicit approximation using Newton method                  
!                                                                       
      DO 90 I = 1, N 
         YKP1 (I) = ZJ (0, I) + ZJ (1, I) + ZJ (2, I) + ZJ (3, I)       &
         + ZJ (4, I)                                                    
   90 END DO 
      CALL DES (XK + HK, YKP1, F) 
      DO 120 K = 1, N 
         DO 100 I = 1, N 
            HELP (I) = YKP1 (I) 
  100    END DO 
         HELP (K) = HELP (K) - HS 
         ! CALL DES (XK + HK, HELP, FS (1, K) ) 
         CALL DES (XK + HK, HELP, FS (:, K) ) 
         DO 110 I = 1, N 
            FS (I, K) = - HK * 0.48D0 * (F (I) - FS (I, K) ) / HS 
  110    END DO 
         FS (K, K) = FS (K, K) + 1.0D0 
  120 END DO 
      NUSED = NUSED+N + 1 
      DO 190 I = 1, N 
         CON (I) = YKP1 (I) - 0.48D0 * (ZJ (1, I) + 2.0D0 * ZJ (2, I)   &
         + 3.0D0 * ZJ (3, I) + 4.0D0 * ZJ (4, I) )                      
         DO 180 K = 1, N 
            FSG (K, I) = FS (K, I) 
  180    END DO 
  190 END DO 
      CALL GAUSSP (N, FSG, NDGL, IPIVOT, MARK, D) 
      IF (MARK.EQ.0) THEN 
         IERR = 6 
         RETURN 
      ENDIF 
      DO 220 ITER = 1, 3 
         DO 210 I = 1, N 
            HELP (I) = - YKP1 (I) 
            DO 200 K = 1, N 
               HELP (I) = HELP (I) + FS (I, K) * YKP1 (K) 
  200       END DO 
            HELP (I) = HK * 0.48D0 * F (I) + HELP (I) + CON (I) 
  210    END DO 
         CALL GAUSSS (N, FSG, NDGL, IPIVOT, HELP, YKP1) 
         CALL DES (XK + HK, YKP1, F) 
  220 END DO 
      NUSED = NUSED+3 
!                                                                       
!** Determine corresponding Gear-Nordsieck approximation                
!                                                                       
      DO 230 I = 1, N 
         HELP (I) = HK * F (I) - ZJ (1, I) - 2.0D0 * ZJ (2, I) - 3.0D0 *&
         ZJ (3, I) - 4.0D0 * ZJ (4, I)                                  
  230 END DO 
      DO 250 I = 1, N 
         ZJP1 (0, I) = YKP1 (I) 
         ZJP1 (1, I) = HK * F (I) 
         ZJP1 (2, I) = ZJ (2, I) + 3.0D0 * ZJ (3, I) + 6.0D0 * ZJ (4, I)&
         + 0.7D0 * HELP (I)                                             
         ZJP1 (3, I) = ZJ (3, I) + 4.0D0 * ZJ (4, I) + 0.2D0 * HELP (I) 
         ZJP1 (4, I) = ZJ (4, I) + 0.02D0 * HELP (I) 
  250 END DO 
!                                                                       
!** Determine whether the last step should be accepted                  
!                                                                       
      DO 260 I = 1, N 
         HELP (I) = ZJP1 (4, I) 
         CON (I) = ZJ (4, I) 
  260 END DO 
      DIFF = MAXVAL(ABS(HELP-CON)) !DVNORM (HELP, CON, N) 
      YMAX = MAXVAL(ABS(YKP1)) !DVNORM (YKP1, Y0, N) 
      EPS = (EPSABS + EPSREL * YMAX) / 6.0D0 
      Q = SQRT (SQRT (EPS / DIFF) ) / 1.2 
      IF (DIFF.LT.EPS) THEN 
!                                                                       
!** Accept last step; prepare for next integration step                 
!                                                                       
         XK = XK + HK 
         DO 270 I = 1, N 
            YK (I) = YKP1 (I) 
  270    END DO 
!                                                                       
!** Jump back if the interval endpoint XE has been reached or           
!** if the right hand side has been called too often.                   
!                                                                       
  275    IF (IEND) THEN 
            HK = DUMMY 
            RETURN 
         ELSEIF (NUSED.GT.NMAX) THEN 
            IERR = 5 
            RETURN 
         ENDIF 
!                                                                       
!** adapt step size for next step                                       
!                                                                       
         HALT = HK 
         HK = MIN (Q, 2.0D0) * HK 
         IF (XK + HK.GE.XEND) THEN 
            DUMMY = HK 
            HK = XE-XK 
            IEND = .TRUE. 
!                                                                       
!** jump back if sufficiently close to XE                               
!                                                                       
            IF (HK.LT.EPS1 * ABS (XE) ) GOTO 275 
         ENDIF 
!                                                                       
!** Set up the Gera-Nordsieck approximation for the next                
!** integration                                                         
!                                                                       
         QUOT1 = HK / HALT 
         QUOT2 = QUOT1 * QUOT1 
         QUOT3 = QUOT2 * QUOT1 
         QUOT4 = QUOT3 * QUOT1 
         DO 280 I = 1, N 
            ZJ (0, I) = ZJP1 (0, I) 
            ZJ (1, I) = QUOT1 * ZJP1 (1, I) 
            ZJ (2, I) = QUOT2 * ZJP1 (2, I) 
            ZJ (3, I) = QUOT3 * ZJP1 (3, I) 
            ZJ (4, I) = QUOT4 * ZJP1 (4, I) 
  280    END DO 
      ELSE 
!                                                                       
!** Repeat last step for a smaller step size                            
!** and modify the Gear-Nordsieck approximation accordingly             
!                                                                       
         HALT = HK 
         HK = MAX (0.5D0, Q) * HK 
         QUOT1 = HK / HALT 
         QUOT2 = QUOT1 * QUOT1 
         QUOT3 = QUOT2 * QUOT1 
         QUOT4 = QUOT3 * QUOT1 
         DO 290 I = 1, N 
            ZJ (1, I) = QUOT1 * ZJ (1, I) 
            ZJ (2, I) = QUOT2 * ZJ (2, I) 
            ZJ (3, I) = QUOT3 * ZJ (3, I) 
            ZJ (4, I) = QUOT4 * ZJ (4, I) 
  290    END DO 
         IEND = .FALSE. 
      ENDIF 
!                                                                       
      GOTO 75 
!                                                                       
    END SUBROUTINE GEAR4                          

    subroutine test_mod_ivp
    use mod_show_matrix
    
    integer, parameter :: neq = 2
    doubleprecision, parameter :: m = 1.0d0, k = 12.0d0, d = 0.3d0
    integer n, nmax, index, nused, ierr
    doubleprecision t0, te , h, y(neq)
    
    index = 1   ! ode method used
    nmax = 50    ! Fun evals allowed
    n  = 180
    
    t0 = 0.0d0
    te = 0.25d0
    h = (te-t0)/n
    
    y(1) = 1.0D0
    y(2) = 0.0D0
    
    print *, t0, y, h
    
    call IVP(t0, h, y, neq, DES, te, 1d-6, 1d-4, index, nmax, nused, ierr)
    
    if( ierr /= 0 ) then
        print *, 'ERROR #', ierr
        error stop
    end if
    print *, te, y, h
    
    contains
        subroutine DES(X,Y,YS)
        doubleprecision X, Y(:), YS(:)
        
            YS(1) = Y(2)
            YS(2) = (-k*Y(1) - d*Y(2))/m
        
        end subroutine
    
    end subroutine
    
    end module