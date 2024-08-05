![BA*)
 
! Code converted by alexij using TO_F90 tool.
! Date: 2024-08-02  Time: 21:38:59
 
![LE*)

SUBROUTINE frkfsy(a,da,n,y,des,h,hmx,abserr,relerr,ierr)
![IX{FRKFSY}*)

!*****************************************************************
!                                                                *
!  A system of ordinary differential equations of 1st order is   *
!  integrated by applying the RUNGE-KUTTA-FEHLBERG method        *
!  [ of order O(H**5) ] with estimates for the local error and   *
!  step size control.                                            *
![BE*)
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
![BA*)
!*****************************************************************
![BE*)

IMPLICIT DOUBLE PRECISION (a-h,o-z)
DOUBLE PRECISION :: yt(10),t(10),r(10),k1(10),k2(10)
DOUBLE PRECISION :: k3(10),k4(10),k5(10),k6(10)
DOUBLE PRECISION :: y(n)

!     determine machine constant

fmachp = 1.0D0
2 fmachp = 0.5D0 * fmachp
IF (machpd(1.0D0+fmachp) == 1) GO TO 2
fmachp = fmachp * 2.0D0

!     check the input data

ierr=3
IF(relerr < 0.0D0 .OR. abserr < 0.0D0 .OR.  &
    relerr+abserr == 0.0D0 .OR. hmx <= 0.0D0) RETURN

ierr=4
b=a+da
IF(DABS(da) <= 13.0D0*fmachp*DMAX1(DABS(a),DABS(b))) RETURN

hmax=DMIN1(hmx,DABS(da))
IF(DABS(h) <= 13.0D0*fmachp*DABS(a)) h=hmax

!     Initialize counter for calls of SUBROUTINE DES

lfd=0
iad=0

!     H is bounded by HMAX and is chosen so that the
!     endpoint B is reached, if possible.

3  h=DSIGN(DMIN1(DABS(h),hmax),da)
IF(DABS(b-a) <= 1.25D0*DABS(h)) THEN
  hf=h
  
!        if IAD=1 and H=B-A acceptable, we stop after
!        the next integration step.
  
  iad=1
  h=b-a
END IF

!     an integration step is executed

CALL des(a,y,k1)
lfd=lfd+1
5  CONTINUE
x=0.25D0*h
DO  i=1,n
  yt(i)=y(i)+x*k1(i)
END DO
x=a+x
CALL des(x,yt,k2)
DO  i=1,n
  yt(i)=y(i)+h*(k1(i)*(3.0D0/32.0D0)+k2(i)*(9.0D0/32.0D0))
END DO
x=a+h*(3.0D0/8.0D0)
CALL des(x,yt,k3)
DO  i=1,n
  yt(i)=y(i)+h*(k1(i)*(1932.0D0/2197.0D0) -k2(i)*(7200.0D0/2197.0D0)  &
      +k3(i)*(7296.0D0/2197.0D0))
END DO
x=a+h*(12.0D0/13.0D0)
CALL des(x,yt,k4)
DO  i=1,n
  yt(i)=y(i)+h*(k1(i)*(439.0D0/216.0D0)-8.0D0*k2(i)  &
      +k3(i)*(3680.0D0/513.0D0) -k4(i)*(845.0D0/4104.0D0))
END DO
x=a+h
CALL des(x,yt,k5)
DO  i=1,n
  yt(i)=y(i)+h*(-k1(i)*(8.0D0/27.0D0)+2.0D0*k2(i)  &
      -k3(i)*(3544.0D0/2565.0D0) +k4(i)*(1859.0D0/4104.0D0)  &
      -k5(i)*(11.0D0/40.0D0))
END DO
x=a+0.5D0*h
CALL des(x,yt,k6)
DO  i=1,n
  t(i)=k1(i)*(25.0D0/216.0D0)+k3(i)*(1408.0D0/2565.0D0)  &
      +k4(i)*(2197.0D0/4104.0D0)-k5(i)*0.20D0
  yt(i)=y(i)+h*t(i)
END DO

!     YT(I) now represents the latest result of this pass.
!     Determine R(I), the estimate of the local
!     error, relative to the current step size.

DO  i=1,n
  r(i)=k1(i)/360.0D0-k3(i)*(128.0D0/4275.0D0)  &
      -k4(i)*(2197.0D0/75240.0D0)+k5(i)/50.0D0 +k6(i)*(2.0D0/55.0D0)
END DO

!     Check accuracy

quot=0.0D0
DO  i=1,n
  tr=DABS(r(i))/(relerr*DABS(yt(i))+abserr)
  quot=DMAX1(quot,tr)
END DO

!     If  QOUOT.LE.1.0   ==> integration step is accepted

IF(quot <= 1.0D0) THEN
  
!        result is accepted
  
  DO  i=1,n
    y(i)=yt(i)
  END DO
  a=a+h
  
!        if A=B ,  RETURN
  
  IF(iad == 1) THEN
    ierr=1
    h=hf
    RETURN
  END IF
  
!        prepare next step
  
  quot=DMAX1(quot,6.5536D-4)
END IF
quot=DMIN1(quot,4096.0D0)
h=0.8D0*h/DSQRT(DSQRT(quot))

!     We just achieved that H was increased by at most a factor of 5,
!     or alternatively, that it was decreased by a factor of 10
!     at most

IF(DABS(h) <= 13.0D0*fmachp*DABS(a)) THEN
  ierr=4
  RETURN
END IF
lfd=lfd+5
IF(lfd >= 2995) THEN
  ierr=2
  RETURN
END IF
IF(quot <= 1.0D0) THEN
  
!        the step was successful. Continue with another step.
  
  GO TO 3
ELSE
  
!        the step is repeated for a smaller H.
  
  iad=0
  GO TO 5
END IF
END
