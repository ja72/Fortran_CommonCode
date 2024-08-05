![BA*)
 
! Code converted by alexij using TO_F90 tool.
! Date: 2024-08-02  Time: 21:39:12
 
![LE*)
![LE*)
![LE*)
![FE{F 17.7.3}
![  {Gear's method for integrating stiff systems of DEs}
![  {Gear's method for integrating stiff systems of DEs}*)
![LE*)

SUBROUTINE gear4 (xk, hk, yk, n, des, xe, epsabs, epsrel, nmax, nused, ierr)
![IX{GEAR4}*)

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
![BE*)
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
!  Subroutines used:  IVP, GAUSSP, GAUSSS, DVNORM, MACHPD        *
!                                                                *
!*****************************************************************
!                                                                *
!  Author      : Klaus Niederdrenk                               *
!  Date        : 1.22.1996                                       *
!  Source code : FORTRAN 77                                      *
!                                                                *
![BA*)
!*****************************************************************
![BE*)

IMPLICIT DOUBLE PRECISION (a-h,o-z)
DOUBLE PRECISION :: yk(1:n)
PARAMETER ( ndgl = 20 )
DIMENSION zj(0:4,1:ndgl), zjp1(0:4,1:ndgl), f(1:ndgl)
DIMENSION fs(1:ndgl,1:ndgl), help(1:ndgl), y0(1:ndgl)
DIMENSION ykp1(1:ndgl), con(1:ndgl)
DIMENSION d(1:ndgl), ipivot(1:ndgl), fsg(1:ndgl,1:ndgl)
LOGICAL :: iend, lepsi
EXTERNAL des
SAVE eps1, eps2, lepsi, y0, hs
DATA lepsi / .true. / , y0 / ndgl * 0.0D0 /

!** Using the machine constant FMACHP, we determine EPS1 in order
!** to avoid excessively small final steps near XE, EPS2 to check
!** for zero and HS as the optimal step size for approximating the
!** Jacobi matrix. (This is done only once at the start.)

IF ( lepsi ) THEN
  fmachp = 1.0D0
  10    fmachp = 0.5 * fmachp
  IF ( machpd(1.0 + fmachp) == 1) GO TO 10
  fmachp = 2.0D0 * fmachp
  eps1   = fmachp ** 0.75D0
  eps2   = 100.0D0 * fmachp
  hs     = 10.0D0 * SQRT(fmachp)
  lepsi  = .false.
END IF

!** Initialize

sg    = DSIGN(1.0D0, xe)
xend  = (1.0D0 - sg*eps2) * xe
ierr  = 0
nused = 0
iend  = .false.

!** Check input parameters

ymax = dvnorm(yk, y0, n)
IF (epsabs <= eps2*ymax .AND. epsrel <= eps2) THEN
  ierr = 1
ELSE IF (xend < xk) THEN
  ierr = 2
ELSE IF (hk < eps2*DABS(xk)) THEN
  ierr = 3
ELSE IF ( n <= 0 .OR. n > ndgl ) THEN
  ierr = 4
END IF
IF (ierr /= 0)  r e t u r n

!****  compute first integration   ****

IF (xk+hk > xend) THEN
  hk = xe - xk
  dummy = hk
  iend = .true.
END IF
DO  i = 1, n
  help(i) = yk(i)
END DO
xka = xk
xke = xka
hka = 0.25*hk
hk1 = hka
DO  k = 1, 4
  xke = xke + hka
  CALL ivp (xka, hk1, help, n, des, xke, epsabs, epsrel,  &
      1, nmax-nused, nanl, ierr)
  nused = nused + nanl
  IF( ierr /= 0 ) r e t u r n
  DO  i = 1, n
    zjp1(k,i) = help(i)
  END DO
END DO
CALL des (xk, yk, n, f)
nused = nused + 1

!** Determine first Gear-Nordsieck approximation

DO  i = 1, n
  zj(0,i) = yk(i)
  zj(1,i) = hk*f(i)
  zj(2,i) = 1.0D0/24.0D0*(35.0D0*yk(i) - 104.0D0*zjp1(1,i)  &
      + 114.0D0*zjp1(2,i) - 56.0D0*zjp1(3,i) +  11.0D0*zjp1(4,i))
  zj(3,i) = 1.0D0/12.0D0*(-5.0D0*yk(i) + 18.0D0*zjp1(1,i)  &
      - 24.0D0*zjp1(2,i) + 14.0D0*zjp1(3,i) -  3.0D0*zjp1(4,i))
  zj(4,i) = 1.0D0/24.0D0*(yk(i) - 4.0D0*zjp1(1,i)  &
      + 6.0D0*zjp1(2,i) - 4.0D0*zjp1(3,i) + zjp1(4,i))
END DO


!****  S t e p  S i z e  A l g o r i t h m   ****


75  CONTINUE

!** Compute implicit approximation using Newton method

DO  i = 1, n
  ykp1(i) = zj(0,i)+zj(1,i)+zj(2,i)+zj(3,i)+zj(4,i)
END DO
CALL des (xk+hk, ykp1, n, f)
DO  k = 1, n
  DO  i = 1, n
    help(i) = ykp1(i)
  END DO
  help(k) = help(k) - hs
  CALL des (xk+hk, help, n, fs(1,k))
  DO  i = 1, n
    fs(i,k) = -hk * 0.48D0 * (f(i) - fs(i,k)) / hs
  END DO
  fs(k,k) = fs(k,k) + 1.0D0
END DO
nused = nused + n + 1
DO  i = 1, n
  con(i) = ykp1(i) - 0.48D0 * ( zj(1,i) + 2.0D0*zj(2,i)  &
      + 3.0D0*zj(3,i) + 4.0D0*zj(4,i) )
  DO  k = 1, n
    fsg(k,i) = fs(k,i)
  END DO
END DO
CALL gaussp (n, fsg, ndgl, ipivot, mark, d)
IF ( mark == 0 ) THEN
  ierr = 6
  r e t u r n
END IF
DO  iter = 1, 3
  DO  i = 1, n
    help(i) = -ykp1(i)
    DO  k = 1, n
      help(i) = help(i) + fs(i,k)*ykp1(k)
    END DO
    help(i) = hk*0.48D0*f(i) + help(i) + con(i)
  END DO
  CALL gausss (n, fsg, ndgl, ipivot, help, ykp1)
  CALL des (xk+hk, ykp1, n, f)
END DO
nused = nused + 3

!** Determine corresponding Gear-Nordsieck approximation

DO  i = 1, n
  help(i) = hk*f(i) - zj(1,i) - 2.0D0*zj(2,i) - 3.0D0*zj(3,i) - 4.0D0*zj(4,i)
END DO
DO  i = 1, n
  zjp1(0,i) = ykp1(i)
  zjp1(1,i) = hk*f(i)
  zjp1(2,i) = zj(2,i) + 3.0D0*zj(3,i) + 6.0D0*zj(4,i) + 0.7D0*help(i)
  zjp1(3,i) = zj(3,i) + 4.0D0*zj(4,i) + 0.2D0*help(i)
  zjp1(4,i) = zj(4,i) + 0.02D0*help(i)
END DO

!** Determine whether the last step should be accepted

DO  i = 1, n
  help(i) = zjp1(4,i)
  con(i)  = zj(4,i)
END DO
diff = dvnorm(help, con, n)
ymax = dvnorm(ykp1, y0, n)
eps = (epsabs + epsrel*ymax) / 6.0D0
q = DSQRT(DSQRT(eps/diff))/1.2
IF ( diff < eps ) THEN
  
!** Accept last step; prepare for next integration step
  
  xk = xk + hk
  DO  i = 1, n
    yk(i) = ykp1(i)
  END DO
  
!** Jump back if the interval endpoint XE has been reached or
!** if the right hand side has been called too often.
  
  275    IF ( iend ) THEN
    hk = dummy
    r e t u r n
  ELSE IF ( nused > nmax ) THEN
    ierr = 5
    r e t u r n
  END IF
  
!** adapt step size for next step
  
  halt = hk
  hk = DMIN1(q, 2.0D0) * hk
  IF ( xk + hk >= xend ) THEN
    dummy = hk
    hk = xe - xk
    iend = .true.
    
!** jump back if sufficiently close to XE
    
    IF ( hk < eps1*DABS(xe) ) GO TO 275
  END IF
  
!** Set up the Gera-Nordsieck approximation for the next
!** integration
  
  quot1 = hk / halt
  quot2 = quot1 * quot1
  quot3 = quot2 * quot1
  quot4 = quot3 * quot1
  DO  i = 1, n
    zj(0,i) = zjp1(0,i)
    zj(1,i) = quot1 * zjp1(1,i)
    zj(2,i) = quot2 * zjp1(2,i)
    zj(3,i) = quot3 * zjp1(3,i)
    zj(4,i) = quot4 * zjp1(4,i)
  END DO
ELSE
  
!** Repeat last step for a smaller step size
!** and modify the Gear-Nordsieck approximation accordingly
  
  halt = hk
  hk = DMAX1(0.5D0, q) * hk
  quot1 = hk / halt
  quot2 = quot1 * quot1
  quot3 = quot2 * quot1
  quot4 = quot3 * quot1
  DO  i = 1, n
    zj(1,i) = quot1 * zj(1,i)
    zj(2,i) = quot2 * zj(2,i)
    zj(3,i) = quot3 * zj(3,i)
    zj(4,i) = quot4 * zj(4,i)
  END DO
  iend = .false.
END IF

GO TO 75

END
