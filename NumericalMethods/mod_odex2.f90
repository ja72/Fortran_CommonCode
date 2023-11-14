    module mod_odex2
    use, intrinsic :: iso_fortran_env
    
    integer,parameter :: r8 = real64, r4 = real32
    
    interface 
    SUBROUTINE f_deriv(n,x,y,yp,rpar,ipar)
    ! --- ODE SYSTEM DERIVATIVES
    import
    implicit double precision (a-h,o-z)
    integer, intent(in)                 :: n
    double precision, intent(in)        :: x
    double precision, intent(in)        :: y(n)
    double precision, intent(out)       :: yp(n)
    real, intent(in out)                :: rpar(*)
    integer, intent(in out)             :: ipar(*)

    end subroutine f_deriv
    end interface

    contains

    SUBROUTINE odex2(n,fcn,x,y,yp,xend,h,  &
        rtol,atol,itol,  &
        solout,iout,  &
        work,lwork,iwork,liwork,rpar,ipar,idid)

    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2023-11-13  Time: 08:50:32

    ! ----------------------------------------------------------
    !     NUMERICAL SOLUTION OF A SYSTEM OF SECOND 0RDER
    !     ORDINARY DIFFERENTIAL EQUATIONS  Y''=F(X,Y).
    !     THIS IS AN EXTRAPOLATION-ALGORITHM, BASED ON
    !     THE STOERMER RULE (WITH STEPSIZE CONTROL
    !     ORDER SELECTION AND DENSE OUTPUT).

    !     AUTHORS: E. HAIRER AND G. WANNER
    !              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
    !              CH-1211 GENEVE 24, SWITZERLAND
    !              E-MAIL:  Ernst.Hairer@math.unige.ch
    !                       Gerhard.Wanner@math.unige.ch

    !     THIS CODE IS DESCRIBED IN SECTION II.14 OF THE BOOK:
    !         E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
    !         DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
    !         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
    !         SPRINGER-VERLAG (1993)

    !     VERSION SEPTEMBER 30, 1995
    !         SMALL CORRECTIONS ON JUNE 11, 1999

    !     INPUT PARAMETERS
    !     ----------------
    !     N           DIMENSION OF THE SYSTEM

    !     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE
    !                 VALUE OF F(X,Y):
    !                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR)
    !                    DOUBLE PRECISION X,Y(N),F(N)
    !                    F(1)=...   ETC.

    !     X           INITIAL X-VALUE

    !     Y(N)        INITIAL VALUES FOR Y

    !     YP(N)       INITIAL VALUES FOR Y'

    !     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE)

    !     H           INITIAL STEP SIZE GUESS;
    !                 USUALLY 1.D-1 OR 1.D-3, IS GOOD. THIS CHOICE IS NOT
    !                 VERY IMPORTANT, THE CODE QUICKLY ADAPTS ITS STEPSIZE.
    !                 WHEN YOU ARE NOT SURE, THEN STUDY THE CHOSEN VALUES
    !                 FOR A FEW STEPS IN SUBROUTINE "SOLOUT".
    !                 (IF H=0.D0, THE CODE PUTS H=1.D-4).

    !     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY CAN
    !                 BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH 2*N.

    !     ITOL        SWITCH FOR RTOL AND ATOL:
    !                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
    !                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
    !                     Y(I)  BELOW  RTOL*ABS(Y(I))+ATOL
    !                     YP(I) BELOW  RTOL*ABS(YP(I))+ATOL
    !                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
    !                     THE CODE KEEPS THE LOCAL ERROR OF
    !                     Y(I)  BELOW  RTOL(I)*ABS(Y(I))+ATOL(I).
    !                     YP(I) BELOW  RTOL(I+N)*ABS(YP(I))+ATOL(I+N).

    !     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
    !                 NUMERICAL SOLUTION DURING INTEGRATION.
    !                 IF IOUT>=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
    !                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0.
    !                 IT MUST HAVE THE FORM
    !                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,YP,N,CON,NCON,ICOMP,ND,
    !                                       RPAR,IPAR,IRTRN)
    !                    DIMENSION Y(N),YP(N),CON(NCON),ICOMP(ND)
    !                    ....
    !                 SOLOUT FURNISHES THE SOLUTIONS "Y, YP" AT THE NR-TH
    !                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
    !                    THE FIRST GRID-POINT).
    !                 "XOLD" IS THE PRECEEDING GRID-POINT.
    !                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
    !                    IS SET <0, ODEX2 WILL RETURN TO THE CALLING PROGRAM.

    !          -----  CONTINUOUS OUTPUT (IF IOUT=2): -----
    !                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
    !                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
    !                 THE DOUBLE PRECISION FUNCTION
    !                    >>>   CONTX2(I,S,CON,NCON,ICOMP,ND)   <<<
    !                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
    !                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
    !                 S SHOULD LIE IN THE INTERVAL [XOLD,X].

    !     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
    !                    IOUT=0: SUBROUTINE IS NEVER CALLED
    !                    IOUT=1: SUBROUTINE IS USED FOR OUTPUT
    !                    IOUT=2: DENSE OUTPUT IS PERFORMED IN SOLOUT

    !     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
    !                 SERVES AS WORKING SPACE FOR ALL VECTORS.
    !                 "LWORK" MUST BE AT LEAST
    !                    N*(2*KM+6)+5*KM+20+(KM*(2*KM+5)+6)*NRDENS
    !                 WHERE NRDENS=IWORK(8) (SEE BELOW) AND
    !                        KM=9                IF IWORK(2)=0
    !                        KM=IWORK(2)         IF IWORK(2).GT.0
    !                 WORK(1),...,WORK(20) SERVE AS PARAMETERS
    !                 FOR THE CODE. FOR STANDARD USE, SET THESE
    !                 PARAMETERS TO ZERO BEFORE CALLING.

    !     LWORK       DECLARED LENGTH OF ARRAY "WORK".

    !     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
    !                 "LIWORK" MUST BE AT LEAST
    !                               KM+20+NRDENS
    !                 IWORK(1),...,IWORK(20) SERVE AS PARAMETERS
    !                 FOR THE CODE. FOR STANDARD USE, SET THESE
    !                 PARAMETERS TO ZERO BEFORE CALLING.

    !     LIWORK      DECLARED LENGTH OF ARRAY "IWORK".

    !     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH
    !                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
    !                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES.

    !-----------------------------------------------------------------------

    !     SOPHISTICATED SETTING OF PARAMETERS
    !     -----------------------------------
    !              SEVERAL PARAMETERS (WORK(1),...,IWORK(1),...) ALLOW
    !              TO ADAPT THE CODE TO THE PROBLEM AND TO THE NEEDS OF
    !              THE USER. FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES.

    !    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 2.3D-16.

    !    WORK(2)   MAXIMAL STEP SIZE, DEFAULT XEND-X.

    !    WORK(3)   STEP SIZE IS REDUCED BY FACTOR WORK(3), IF DURING THE
    !              COMPUTATION OF THE EXTRAPOLATION TABLEAU DIVERGENCE
    !              IS OBSERVED; DEFAULT 0.5.

    !    WORK(4), WORK(5)   PARAMETERS FOR STEP SIZE SELECTION
    !              THE NEW STEP SIZE FOR THE J-TH DIAGONAL ENTRY IS
    !              CHOSEN SUBJECT TO THE RESTRICTION
    !                 FACMIN/WORK(5) <= HNEW(J)/HOLD <= 1/FACMIN
    !              WHERE FACMIN=WORK(4)**(1/(2*J-1))
    !              DEFAULT VALUES: WORK(4)=0.02D0, WORK(5)=4.D0

    !    WORK(6), WORK(7)   PARAMETERS FOR THE ORDER SELECTION
    !              STEP SIZE IS DECREASED IF    W(K-1) <= W(K)*WORK(6)
    !              STEP SIZE IS INCREASED IF    W(K) <= W(K-1)*WORK(7)
    !              DEFAULT VALUES: WORK(6)=0.8D0, WORK(7)=0.9D0

    !    WORK(8), WORK(9)   SAFETY FACTORS FOR STEP CONTROL ALGORITHM
    !             HNEW=H*WORK(9)*(WORK(8)*TOL/ERR)**(1/(J-1))
    !             DEFAULT VALUES: WORK(8)=0.65D0,
    !                        WORK(9)=0.94D0  IF "HOPE FOR CONVERGENCE"

    !    IWORK(1)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
    !              THE DEFAULT VALUE (FOR IWORK(1)=0) IS 10000.

    !    IWORK(2)  THE MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION
    !              TABLE. THE DEFAULT VALUE (FOR IWORK(2)=0) IS 9.
    !              IF IWORK(2).NE.0 THEN IWORK(2) SHOULD BE .GE.3.

    !    IWORK(3)  SWITCH FOR THE STEP SIZE SEQUENCE (EVEN NUMBERS ONLY)
    !              IF IWORK(3).EQ.1 THEN 2,4,6,8,10,12,14,16,...
    !              IF IWORK(3).EQ.2 THEN 2,4,8,12,16,20,24,28,...
    !              IF IWORK(3).EQ.3 THEN 2,4,6,8,12,16,24,32,...
    !              IF IWORK(3).EQ.4 THEN 2,6,10,14,18,22,26,30,...
    !              THE DEFAULT VALUE IS IWORK(3)=1 IF IOUT.LE.1;
    !              THE DEFAULT VALUE IS IWORK(3)=4 IF IOUT.GE.2.

    !    IWORK(6)  IF  IWORK(6)=0  ERROR ESTIMATOR IN THE DENSE
    !              OUTPUT FORMULA IS ACTIVATED. IT CAN BE SUPPRESSED
    !              BY PUTTING IWORK(6)=1.
    !              DEFAULT IWORK(6)=0  (IF IOUT.GE.2).

    !    IWORK(7)  DETERMINES THE DEGREE OF INTERPOLATION FORMULA
    !              MU = 2 * KAPPA - IWORK(7) + 1
    !              IWORK(7) SHOULD LIE BETWEEN 1 AND 8
    !              DEFAULT IWORK(7)=6  (IF IWORK(7)=0).

    !    IWORK(8)  = NRDENS = NUMBER OF COMPONENTS, FOR WHICH DENSE OUTPUT
    !              IS REQUIRED

    !    IWORK(21),...,IWORK(NRDENS+20) INDICATE THE COMPONENTS, FOR WHICH
    !              DENSE OUTPUT IS REQUIRED. THEY NEED NOT BE SET IF
    !              NRDENS = N .

    !----------------------------------------------------------------------C
    !     OUTPUT PARAMETERS
    !     -----------------
    !     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
    !                 (AFTER SUCCESSFUL RETURN X=XEND).

    !     Y(N)        NUMERICAL SOLUTION AT X

    !     YP(N)       NUMERICAL DERIVATIVE OF SOLUTION AT X

    !     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP

    !     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
    !                   IDID=1  COMPUTATION SUCCESSFUL,
    !                   IDID=-1 COMPUTATION UNSUCCESSFUL.

    !   IWORK(17)  NFCN    NUMBER OF FUNCTION EVALUATIONS
    !   IWORK(18)  NSTEP   NUMBER OF COMPUTED STEPS
    !   IWORK(19)  NACCPT  NUMBER OF ACCEPTED STEPS
    !   IWORK(20)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
    !                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED)
    !-----------------------------------------------------------------------
    ! *** *** *** *** *** *** *** *** *** *** *** *** ***
    !          DECLARATIONS
    ! *** *** *** *** *** *** *** *** *** *** *** *** ***
    implicit double precision (a-h,o-z)    

    INTEGER, INTENT(IN)                  :: n
    procedure(f_deriv), intent(in), pointer :: fcn
    DOUBLE PRECISION, INTENT(IN OUT)     :: x
    DOUBLE PRECISION, INTENT(IN OUT)     :: y(n)
    DOUBLE PRECISION, INTENT(IN OUT)     :: yp(n)
    DOUBLE PRECISION, INTENT(IN)         :: xend
    DOUBLE PRECISION, INTENT(IN OUT)     :: h
    REAL, INTENT(IN)         :: rtol(*)
    REAL, INTENT(IN)         :: atol(*)
    INTEGER, INTENT(IN OUT)              :: itol
    DOUBLE PRECISION, INTENT(IN OUT)     :: solout
    INTEGER, INTENT(IN OUT)              :: iout
    INTEGER, INTENT(IN)                  :: lwork
    INTEGER, INTENT(IN)                  :: liwork
    DOUBLE PRECISION, INTENT(IN OUT)     :: work(lwork)
    INTEGER, INTENT(IN OUT)              :: iwork(liwork)
    REAL, INTENT(IN OUT)                 :: rpar(*)
    INTEGER, INTENT(IN OUT)              :: ipar(*)
    INTEGER, INTENT(OUT)                 :: idid

    integer :: nfcn, nstep, naccpt, nrejct, nmax, km, nsequ, iderr, mudif, nrdens
    integer :: i, lfsafe, iedy, ieyh1, ieyh2, iedz, iescal, iescp, iet, ietp, nrd
    integer :: iefs, ieys, iehh, iew, iea, iefac, ieco, istore, icom, ienj, ncom
    integer :: ieysp
    DOUBLE PRECISION :: uround, hmax, safe1, safe2, safe3, fac1, fac2, fac3, fac4
    LOGICAL :: arret    
    
    integer, parameter :: max_steps = 10000
    
    ! *** *** *** *** *** *** ***
    !        SETTING THE PARAMETERS
    ! *** *** *** *** *** *** ***
    nfcn=0
    nstep=0
    naccpt=0
    nrejct=0
    arret=.false.
    ! -------- NMAX , THE MAXIMAL NUMBER OF STEPS -----
    IF(iwork(1) == 0)THEN
        nmax=max_steps
    ELSE
        nmax=iwork(1)
        IF(nmax <= 0)THEN
            WRITE(6,*)' WRONG INPUT IWORK(1)=',iwork(1)
            arret=.true.
        END IF
    END IF
    ! -------- KM     MAXIMUM NUMBER OF COLUMNS IN THE EXTRAPOLATION
    IF(iwork(2) == 0)THEN
        km=9
    ELSE
        km=iwork(2)
        IF(km <= 2)THEN
            WRITE(6,*)' CURIOUS INPUT IWORK(2)=',iwork(2)
            arret=.true.
        END IF
    END IF
    ! -------- NSEQU     CHOICE OF STEP SIZE SEQUENCE
    nsequ=iwork(3)
    IF(iwork(3) == 0.AND.iout <= 1) nsequ=1
    IF(iwork(3) == 0.AND.iout >= 2) nsequ=4
    IF(nsequ <= 0.OR.nsequ >= 5)THEN
        WRITE(6,*)' CURIOUS INPUT IWORK(3)=',iwork(3)
        arret=.true.
    END IF
    IF (nsequ <= 3.AND.iout >= 2) THEN
        WRITE(6,*)' IWORK(3) NOT COMPATIBLE WITH IOUT'
        arret=.true.
    END IF
    ! -------- IDERR  PARAMETER FOR ERROR ESTIMATION IN DENSE OUTPUT
    IF(iwork(6) == 0)THEN
        IF(iout <= 1) iderr=1
        IF(iout >= 2) iderr=0
    ELSE
        iderr=iwork(6)
        IF(iout <= 1)THEN
            WRITE(6,*)' ERROR ESTIMATION IN DENSE OUTPUT',  &
                ' NOT POSSIBLE, WRONG IWORK(6)=',iwork(6)
            arret=.true.
        END IF
    END IF
    ! -------- MUDIF
    IF(iwork(7) == 0)THEN
        mudif=6
    ELSE
        mudif=iwork(7)
        IF(mudif <= 0.OR.mudif >= 9)THEN
            WRITE(6,*)' WRONG INPUT IWORK(7)=',iwork(7)
            arret=.true.
        END IF
    END IF
    ! -------- NRDENS   NUMBER OF DENSE OUTPUT COMPONENTS
    nrdens=iwork(8)
    IF(nrdens < 0.OR.nrdens > n)THEN
        WRITE(6,*)' CURIOUS INPUT IWORK(8)=',iwork(8)
        arret=.true.
    END IF
    IF (nrdens == n) THEN
        DO  i=1,nrdens
            iwork(20+i)=i
        END DO
    END IF
    ! -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0
    IF(work(1) == 0.d0)THEN
        uround=2.3D-16
    ELSE
        uround=work(1)
        IF(uround <= 1.d-35.OR.uround >= 1.d0)THEN
            WRITE(6,*)' WHICH MACHINE DO YOU HAVE? YOUR UROUND WAS:'  &
                ,work(1)
            arret=.true.
        END IF
    END IF
    ! -------- MAXIMAL STEP SIZE
    IF(work(2) == 0.d0)THEN
        hmax=xend-x
    ELSE
        hmax=ABS(work(2))
    END IF
    ! -------- STEP SIZE REDUCTION FACTOR
    IF(work(3) == 0.d0)THEN
        safe3=0.5D0
    ELSE
        safe3=work(3)
        IF(safe3 <= uround.OR.safe3 >= 1.d0)THEN
            WRITE(6,*)' CURIOUS INPUT WORK(3)=',work(3)
            arret=.true.
        END IF
    END IF
    ! -------  FAC1,FAC2     PARAMETERS FOR STEP SIZE SELECTION
    IF(work(4) == 0.d0)THEN
        fac1=0.02D0
    ELSE
        fac1=work(4)
    END IF
    IF(work(5) == 0.d0)THEN
        fac2=4.0D0
    ELSE
        fac2=work(5)
    END IF
    ! -------  FAC3, FAC4   PARAMETERS FOR THE ORDER SELECTION
    IF(work(6) == 0.d0)THEN
        fac3=0.8D0
    ELSE
        fac3=work(6)
    END IF
    IF(work(7) == 0.d0)THEN
        fac4=0.9D0
    ELSE
        fac4=work(7)
    END IF
    ! ------- SAFE1, SAFE2 SAFETY FACTORS FOR STEP SIZE PREDICTION
    IF(work(8) == 0.d0)THEN
        safe1=0.65D0
    ELSE
        safe1=work(8)
    END IF
    IF(work(9) == 0.d0)THEN
        safe2=0.94D0
    ELSE
        safe2=work(9)
    END IF
    ! ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK -----
    lfsafe=2*km*km+km
    iedy=21
    ieyh1=iedy+n
    ieyh2=ieyh1+n
    iedz=ieyh2+n
    iescal=iedz+n
    iescp=iescal+n
    iet=iescp+n
    ietp=iet+km*n
    iefs=ietp+km*n
    ieys=iefs+lfsafe*nrdens
    ieysp=ieys+km*nrdens
    iehh=ieysp+km*nrdens
    iew=iehh+km
    iea=iew+km
    iefac=iea+km
    ! ------ TOTAL STORAGE REQUIREMENT -----------
    ieco=iefac+2*km
    istore=ieco+(2*km+6)*nrdens-1
    IF(istore > lwork)THEN
        WRITE(6,*)' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',istore
        arret=.true.
    END IF
    ! ------- ENTRY POINTS FOR INTEGER WORKSPACE -----
    icom=21
    ienj=icom+nrdens
    ! --------- TOTAL REQUIREMENT ---------------
    istore=ienj+km-1
    IF(istore > liwork)THEN
        WRITE(6,*)' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=',istore
        arret=.true.
    END IF
    ! ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1
    IF (arret) THEN
        idid=-1
        RETURN
    END IF
    ! -------- CALL TO CORE INTEGRATOR ------------
    nrd=MAX(1,nrdens)
    ncom=MAX(1,(2*km+6)*nrdens)
    CALL odx2co(n,fcn,x,y,yp,xend,hmax,h,rtol,atol,itol,km,  &
        solout,iout,idid,nmax,uround,work(iedy),work(ieyh1),  &
        work(ieyh2),work(iedz),work(iescal),work(iescp),work(iefs),  &
        work(ieys),work(ieysp),work(iet),work(ietp),work(iehh),  &
        work(iew),work(iea),work(ieco),ncom,iwork(icom),iwork(ienj),  &
        nsequ,lfsafe,safe1,safe2,safe3,fac1,fac2,fac3,fac4,iderr,  &
        work(iefac),mudif,nrd,rpar,ipar,nfcn,nstep,naccpt,nrejct)
    iwork(17)=nfcn
    iwork(18)=nstep
    iwork(19)=naccpt
    iwork(20)=nrejct
    ! ----------- RETURN -----------
    RETURN
    END SUBROUTINE odex2



    !  ----- ... AND HERE IS THE CORE INTEGRATOR  ----------

    SUBROUTINE odx2co(n,fcn,x,y,yp,xend,hmax,h,rtol,atol,itol,km,  &
        solout,iout,idid,nmax,uround,dy,yh1,yh2,dz,scal,scalp,fsafe,  &
        ysafe,ysafep,t,tp,hh,w,a,dens,ncom,icomp,nj,nsequ,lfsafe,  &
        safe1,safe2,safe3,fac1,fac2,fac3,fac4,iderr,errfac,mudif,  &
        nrd,rpar,ipar,nfcn,nstep,naccpt,nrejct)
    ! ----------------------------------------------------------
    !     CORE INTEGRATOR FOR ODEX2
    !     PARAMETERS SAME AS IN ODEX2 WITH WORKSPACE ADDED
    ! ----------------------------------------------------------
    !         DECLARATIONS
    ! ----------------------------------------------------------

    IMPLICIT DOUBLE PRECISION (a-h,o-z)
    INTEGER, INTENT(IN)                      :: n
    procedure(f_deriv), pointer              :: fcn
    DOUBLE PRECISION, INTENT(IN OUT)         :: x
    DOUBLE PRECISION, INTENT(IN OUT)         :: y(n)
    DOUBLE PRECISION, INTENT(IN OUT)         :: yp(n)
    DOUBLE PRECISION, INTENT(IN)             :: xend
    DOUBLE PRECISION, INTENT(OUT)            :: hmax
    DOUBLE PRECISION, INTENT(OUT)            :: h
    REAL , INTENT(IN)                        :: rtol(*)
    REAL , INTENT(IN)                        :: atol(*)
    INTEGER, INTENT(IN OUT)                  :: itol
    INTEGER, INTENT(IN)                      :: km
    DOUBLE PRECISION, INTENT(IN OUT)         :: solout
    INTEGER, INTENT(IN OUT)                  :: iout
    INTEGER, INTENT(OUT)                     :: idid
    INTEGER, INTENT(IN OUT)                  :: nmax
    DOUBLE PRECISION, INTENT(IN OUT)         :: uround
    DOUBLE PRECISION, INTENT(IN OUT)         :: dy(n)
    DOUBLE PRECISION, INTENT(IN OUT)         :: yh1(n)
    DOUBLE PRECISION, INTENT(IN OUT)         :: yh2(n)
    DOUBLE PRECISION, INTENT(IN OUT)         :: dz(n)
    DOUBLE PRECISION, INTENT(IN OUT)         :: scal(n)
    DOUBLE PRECISION, INTENT(IN OUT)         :: scalp(n)
    INTEGER, INTENT(IN OUT)                  :: lfsafe
    INTEGER, INTENT(IN OUT)                  :: nrd
    DOUBLE PRECISION, INTENT(IN OUT)         :: fsafe(lfsafe,nrd)
    DOUBLE PRECISION, INTENT(OUT)            :: ysafe(km,nrd)
    DOUBLE PRECISION, INTENT(OUT)            :: ysafep(km,nrd)
    DOUBLE PRECISION, INTENT(IN OUT)         :: t(km,n)
    DOUBLE PRECISION, INTENT(IN OUT)         :: tp(km,n)
    DOUBLE PRECISION, INTENT(IN OUT)         :: hh(km)
    DOUBLE PRECISION, INTENT(OUT)            :: w(km)
    DOUBLE PRECISION, INTENT(OUT)            :: a(km)
    INTEGER, INTENT(IN OUT)                  :: ncom
    DOUBLE PRECISION, INTENT(OUT)            :: dens(ncom)
    INTEGER, INTENT(IN OUT)                  :: icomp(nrd)
    INTEGER, INTENT(IN OUT)                  :: nj(km)
    INTEGER, INTENT(IN OUT)                  :: nsequ
    DOUBLE PRECISION, INTENT(IN OUT)         :: safe1
    DOUBLE PRECISION, INTENT(IN OUT)         :: safe2
    DOUBLE PRECISION, INTENT(IN OUT)         :: safe3
    DOUBLE PRECISION, INTENT(IN OUT)         :: fac1
    DOUBLE PRECISION, INTENT(IN OUT)         :: fac2
    DOUBLE PRECISION, INTENT(IN OUT)         :: fac3
    DOUBLE PRECISION, INTENT(IN OUT)         :: fac4
    INTEGER, INTENT(IN OUT)                  :: iderr
    DOUBLE PRECISION, INTENT(OUT)            :: errfac(2*km)
    INTEGER, INTENT(IN)                      :: mudif
    REAL, INTENT(IN OUT)                     :: rpar(*)
    INTEGER, INTENT(IN OUT)                  :: ipar(*)
    INTEGER, INTENT(OUT)                     :: nfcn
    INTEGER, INTENT(OUT)                     :: nstep
    INTEGER, INTENT(OUT)                     :: naccpt
    INTEGER, INTENT(OUT)                     :: nrejct
    LOGICAL :: reject,last,atov
    integer :: i, j, k, l, kk, kc
    integer :: irtrn, mu, nnrd, ipt, kmit, kmi, kbeg, krn, kopt, lend, lbeg
    double precision :: posneg, errx, prod, xold, errold, hoptde, fac, err
    double precision :: h2, xoldd, hhh, dblenj, factor, facnj, errint

    COMMON /conod2/xoldd,hhh,kmit
    ! --- DEFINE THE STEP SIZE SEQUENCE
    IF (nsequ == 1) THEN
        DO  i=1,km
            nj(i)=2*i
        END DO
    END IF
    IF (nsequ == 2) THEN
        nj(1)=2
        DO  i=2,km
            nj(i)=4*i-4
        END DO
    END IF
    IF (nsequ == 3) THEN
        nj(1)=2
        nj(2)=4
        nj(3)=6
        DO  i=4,km
            nj(i)=2*nj(i-2)
        END DO
    END IF
    IF (nsequ == 4) THEN
        DO  i=1,km
            nj(i)=4*i-2
        END DO
    END IF
    ! --- DEFINE THE A(I) FOR ORDER SELECTION
    a(1)=1.d0+nj(1)/2
    DO  i=2,km
        a(i)=a(i-1)+nj(i)/2
    END DO
    ! --- INITIAL SCALING
    IF (itol == 0) THEN
        DO  i=1,n
            scal(i)=atol(1)+rtol(1)*ABS(y(i))
            scalp(i)=atol(1)+rtol(1)*ABS(yp(i))
        END DO
    ELSE
        DO  i=1,n
            scal(i)=atol(i)+rtol(i)*ABS(y(i))
            scalp(i)=atol(i+n)+rtol(i+n)*ABS(yp(i))
        END DO
    END IF
    ! --- INITIAL PREPARATIONS
    posneg=SIGN(1.d0,xend-x)
    k=MAX(2,MIN(km-1,INT(-LOG10(rtol(1))*0.6D0+1.5D0)))
    hmax=ABS(hmax)
    h=MAX(ABS(h),1.d-4)
    h=posneg*MIN(h,hmax,ABS(xend-x))
    IF (iout >= 1) THEN
        IF (iout >= 2) THEN
            DO  mu=1,km*2
                errx=SQRT(mu/(mu+6.d0))*0.5D0
                prod=(1.5D0/(mu+6.d0))**3
                DO  j=1,mu
                    prod=prod*errx/j
                END DO
                errfac(mu)=prod
            END DO
        END IF
        nnrd=nrd
        irtrn=0
        xold=x
        CALL sol_out(naccpt+1,xold,x,y,yp,n,dens,ncom,icomp,nrd,rpar,ipar,irtrn)
        
        IF (irtrn < 0) GO TO 120
    END IF
    ERR=0.d0
    errold=1.d10
    hoptde=posneg*hmax
    w(1)=0.d0
    reject=.false.
    last=.false.
10  atov=.false.
    ! --- IS XEND REACHED IN THE NEXT STEP?
    IF (0.1D0*ABS(xend-x) <= ABS(x)*uround)GO TO 110
    h=posneg*MIN(ABS(h),ABS(xend-x),hmax,ABS(hoptde))
    IF ((x+1.01D0*h-xend)*posneg > 0.d0) THEN
        h=xend-x
        last=.true.
    END IF
    IF (nstep == 0.OR.iout /= 2) CALL fcn(n,x,y,dz,rpar,ipar)
    nfcn=nfcn+1
    ! --- THE FIRST AND LAST STEP
    IF (nstep == 0.OR.last) THEN
        ipt=0
        nstep=nstep+1
        DO  j=1,k
            kc=j
            CALL stoerm(j,x,y,yp,h,hmax,n,fcn,dy,yh1,yh2,dz,t,tp,nj,hh,w,  &
                err,fac,a,safe1,uround,fac1,fac2,safe2,scal,scalp,atov,  &
                safe3,reject,km,rtol,atol,itol,errold,fsafe,  &
                lfsafe,iout,ipt,ysafe,ysafep,icomp,nrd,rpar,ipar,nfcn)
            IF (atov) GO TO 10
            IF (j > 1.AND.ERR <= 1.d0) GO TO 60
        END DO
        GO TO 55
    END IF
    ! --- BASIC INTEGRATION STEP
30  CONTINUE
    ipt=0
    nstep=nstep+1
    IF (nstep >= nmax) GO TO 120
    kc=k-1
    DO  j=1,kc
        CALL stoerm(j,x,y,yp,h,hmax,n,fcn,dy,yh1,yh2,dz,t,tp,nj,hh,w,  &
            ERR,fac,a,safe1,uround,fac1,fac2,safe2,scal,scalp,atov,safe3,  &
            reject,km,rtol,atol,itol,errold,fsafe,lfsafe,  &
            iout,ipt,ysafe,ysafep,icomp,nrd,rpar,ipar,nfcn)
        IF (atov) GO TO 10
    END DO
    ! --- CONVERGENCE MONITOR
    IF (k == 2.OR.reject) GO TO 50
    IF (ERR <= 1.d0) GO TO 60
    IF (ERR > ((nj(k+1)*nj(k))/4.d0)**2) GO TO 100
50  CONTINUE
    CALL stoerm(k,x,y,yp,h,hmax,n,fcn,dy,yh1,yh2,dz,t,tp,nj,hh,w,  &
        ERR,fac,a,safe1,uround,fac1,fac2,safe2,scal,scalp,atov,safe3,  &
        reject,km,rtol,atol,itol,errold,fsafe,lfsafe,  &
        iout,ipt,ysafe,ysafep,icomp,nrd,rpar,ipar,nfcn)
    IF (atov) GO TO 10
    kc=k
    IF (ERR <= 1.d0) GO TO 60
    ! --- HOPE FOR CONVERGENCE IN LINE K+1
55  CONTINUE
    IF (ERR > (nj(k+1)/2.d0)**2) GO TO 100
    kc=k+1
    CALL stoerm(kc,x,y,yp,h,hmax,n,fcn,dy,yh1,yh2,dz,t,tp,nj,hh,w,  &
        ERR,fac,a,safe1,uround,fac1,fac2,safe2,scal,scalp,atov,safe3,  &
        reject,km,rtol,atol,itol,errold,fsafe,lfsafe,  &
        iout,ipt,ysafe,ysafep,icomp,nrd,rpar,ipar,nfcn)
    IF (atov) GO TO 10
    IF (ERR > 1.d0) GO TO 100
    ! --- STEP IS ACCEPTED
60  xold=x
    x=x+h
    IF (iout >= 2) THEN
        ! ---  KMIT = MU OF THE PAPER
        kmit=MAX(1,2*kc-mudif+1)
        DO  i=1,n
            yh1(i)=t(1,i)
        END DO
        CALL fcn(n,x,yh1,yh2,rpar,ipar)
        h2=h*h
        DO  i=1,nrd
            dens(i)=y(icomp(i))
            dens(nrd+i)=h*yp(icomp(i))
            dens(2*nrd+i)=h2*dz(icomp(i))
            dens(3*nrd+i)=yh1(icomp(i))
            dens(4*nrd+i)=h*tp(1,icomp(i))
            dens(5*nrd+i)=h2*yh2(icomp(i))
        END DO
        xoldd=xold
        hhh=h
        ! --- COMPUTE SOLUTION AND FIRST DERIVATIVE AT MID-POINT ----
        DO  j=2,kc
            dblenj=nj(j)
            DO  l=j,2,-1
                factor=(dblenj/nj(l-1))**2-1
                DO  i=1,nrd
                    ysafe(l-1,i)=ysafe(l,i)+(ysafe(l,i)-ysafe(l-1,i))/factor
                    ysafep(l-1,i)=ysafep(l,i)+(ysafep(l,i)-ysafep(l-1,i))/factor
                END DO
            END DO
        END DO
        DO  i=1,nrd
            dens(6*nrd+i)=ysafe(1,i)
            dens(7*nrd+i)=h*ysafep(1,i)
        END DO

        DO  kmi=2,kmit
            ! --- COMPUTE KMI-TH DERIVATIVE AT MID-POINT ----
            kbeg=kmi/2
            IF (kmi == 2*kbeg) THEN
                IF (kmi == 2) THEN
                    DO  i=1,nrd
                        ysafe(1,i)=(dz(icomp(i))+fsafe(1,i))/2.d0
                    END DO
                    kbeg=2
                END IF
                DO  kk=kbeg,kc
                    facnj=0.5D0*(nj(kk)/2.d0)**(kmi-2)
                    ipt=(kk-1)**2+kk+kmi/2-2
                    DO  i=1,nrd
                        ysafe(kk,i)=(fsafe(ipt,i)+fsafe(ipt+1,i))*facnj
                    END DO
                END DO
            ELSE
                DO  kk=kbeg,kc
                    facnj=(nj(kk)/2.d0)**(kmi-2)
                    ipt=(kk-1)**2+kk+kbeg-1
                    DO  i=1,nrd
                        ysafe(kk,i)=fsafe(ipt,i)*facnj
                    END DO
                END DO
            END IF
            !------          EXTRAPOLATION
            DO  j=kbeg+1,kc
                dblenj=nj(j)
                DO  l=j,kbeg+1,-1
                    factor=(dblenj/nj(l-1))**2-1.d0
                    DO  i=1,nrd
                        ysafe(l-1,i)=ysafe(l,i)+(ysafe(l,i)-ysafe(l-1,i))/factor
                    END DO
                END DO
            END DO
            krn=(kmi+6)*nrd
            DO  i=1,nrd
                dens(krn+i)=ysafe(kbeg,i)*h2
            END DO
            IF (kmi == kmit) CYCLE
            ! ------          COMPUTE DIFFERENCES
            DO  kk=(kmi+1)/2,kc
                lbeg=(kk-1)**2+kmi-1
                lend=kk**2
                IF (kmi == 2) lbeg=lbeg+1
                DO  l=lend,lbeg,-1
                    DO  i=1,nrd
                        fsafe(l,i)=fsafe(l,i)-fsafe(l-1,i)
                    END DO
                END DO
                IF (kmi == 2) THEN
                    l=lbeg-1
                    DO  i=1,nrd
                        fsafe(l,i)=fsafe(l,i)-dz(icomp(i))
                    END DO
                END IF
            END DO
        END DO
        CALL intpo2(nrd,dens,kmit)
        ! --- ESTIMATION OF INTERPOLATION ERROR
        IF (iderr == 0.AND.kmit >= 1) THEN
            errint=0.d0
            DO  i=1,nrd
                errint=errint+(dens((kmit+6)*nrd+i)/scal(icomp(i)))**2
            END DO
            errint=SQRT(errint/nrd)*errfac(kmit)
            hoptde=h/MAX((errint)**(1.d0/(kmit+6)),0.01D0)
            IF (errint > 10.d0) THEN
                h=hoptde
                x=xold
                nrejct=nrejct+1
                reject=.true.
                GO TO 10
            END IF
        END IF
        DO  i=1,n
            dz(i)=yh2(i)
        END DO
    END IF
    DO  i=1,n
        y(i)=t(1,i)
        yp(i)=tp(1,i)
    END DO
    naccpt=naccpt+1
    IF (iout >= 1) THEN
        CALL sol_out (naccpt+1,xold,x,y,yp,n,dens,ncom,icomp,nrd, rpar,ipar,irtrn)
        IF (irtrn < 0) GO TO 120
    END IF
    ! --- COMPUTE OPTIMAL ORDER
    IF (kc == 2) THEN
        kopt=MIN(3,km-1)
        IF (reject) kopt=2
        GO TO 80
    END IF
    IF (kc <= k) THEN
        kopt=kc
        IF (w(kc-1) < w(kc)*fac3) kopt=kc-1
        IF (w(kc) < w(kc-1)*fac4) kopt=MIN(kc+1,km-1)
    ELSE
        kopt=kc-1
        IF (kc > 3.AND.w(kc-2) < w(kc-1)*fac3) kopt=kc-2
        IF (w(kc) < w(kopt)*fac4) kopt=MIN(kc,km-1)
    END IF
    ! --- AFTER A REJECTED STEP
80  IF (reject) THEN
        k=MIN(kopt,kc)
        h=posneg*MIN(ABS(h),ABS(hh(k)))
        reject=.false.
        GO TO 10
    END IF
    ! --- COMPUTE STEPSIZE FOR NEXT STEP
    IF (kopt <= kc) THEN
        h=hh(kopt)
    ELSE
        IF (kc < k.AND.w(kc) < w(kc-1)*fac4) THEN
            h=hh(kc)*a(kopt+1)/a(kc)
        ELSE
            h=hh(kc)*a(kopt)/a(kc)
        END IF
    END IF
    k=kopt
    h=posneg*ABS(h)
    IF (last) GO TO 110
    GO TO 10
    ! --- STEP IS REJECTED
100 CONTINUE
    k=MIN(k,kc,km-1)
    IF (k > 2.AND.w(k-1) < w(k)*fac3) k=k-1
    nrejct=nrejct+1
    h=posneg*hh(k)
    reject=.true.
    last=.false.
    GO TO 30
    ! --- SOLUTION EXIT
110 CONTINUE
    idid=1
    RETURN
    ! --- FAIL EXIT
120 WRITE (6,979) x,h
979 FORMAT(' EXIT OF ODEX2 AT X=',d14.7,'   H=',d14.7)
    idid=-1
    RETURN
    END SUBROUTINE odx2co

    SUBROUTINE stoerm(j,x,y,yp,h,hmax,n,fcn,dy,yh1,yh2,dz,t,tp,nj,  &
        hh,w,ERR,fac,a,safe1,uround,fac1,fac2,safe2,scal,scalp,atov,  &
        safe3,reject,km,rtol,atol,itol,errold,fsafe,  &
        lfsafe,iout,ipt,ysafe,ysafep,icomp,nrd,rpar,ipar,nfcn)
    ! --- THIS SUBROUTINE COMPUTES THE J-TH LINE OF THE
    ! --- EXTRAPOLATION TABLE AND PROVIDES AN ESTIMATION
    ! --- OF THE OPTIMAL STEPSIZE

    IMPLICIT DOUBLE PRECISION (a-h,o-z)
    INTEGER, INTENT(IN)                      :: j
    INTEGER, INTENT(IN)                      :: n
    DOUBLE PRECISION, INTENT(IN OUT)         :: x
    DOUBLE PRECISION, INTENT(IN OUT)         :: y(n)
    DOUBLE PRECISION, INTENT(IN OUT)         :: yp(n)
    DOUBLE PRECISION, INTENT(IN OUT)         :: h
    DOUBLE PRECISION, INTENT(IN OUT)         :: hmax
    INTEGER, INTENT(IN)                      :: km
    procedure(f_deriv), pointer              :: fcn
    DOUBLE PRECISION, INTENT(IN OUT)         :: dy(n)
    DOUBLE PRECISION, INTENT(OUT)            :: yh1(n)
    DOUBLE PRECISION, INTENT(OUT)            :: yh2(n)
    DOUBLE PRECISION, INTENT(IN)             :: dz(n)
    DOUBLE PRECISION, INTENT(OUT)            :: t(km,n)
    DOUBLE PRECISION, INTENT(OUT)            :: tp(km,n)
    INTEGER, INTENT(IN)                      :: nj(km)
    DOUBLE PRECISION, INTENT(OUT)            :: hh(km)
    DOUBLE PRECISION, INTENT(OUT)            :: w(km)
    DOUBLE PRECISION, INTENT(OUT)            :: ERR
    DOUBLE PRECISION, INTENT(OUT)            :: fac
    DOUBLE PRECISION, INTENT(IN)             :: a(km)
    DOUBLE PRECISION, INTENT(IN OUT)         :: safe1
    DOUBLE PRECISION, INTENT(IN)             :: uround
    DOUBLE PRECISION, INTENT(IN)             :: fac1
    DOUBLE PRECISION, INTENT(IN OUT)         :: fac2
    DOUBLE PRECISION, INTENT(IN OUT)         :: safe2
    DOUBLE PRECISION, INTENT(OUT)            :: scal(n)
    DOUBLE PRECISION, INTENT(OUT)            :: scalp(n)
    LOGICAL, INTENT(OUT)                     :: atov
    DOUBLE PRECISION, INTENT(IN)             :: safe3
    LOGICAL, INTENT(OUT)                     :: reject
    REAL , INTENT(IN)                        :: rtol(*)
    REAL , INTENT(IN)                        :: atol(*)
    INTEGER, INTENT(IN OUT)                  :: itol
    DOUBLE PRECISION, INTENT(OUT)            :: errold
    INTEGER, INTENT(IN OUT)                  :: lfsafe
    INTEGER, INTENT(IN)                      :: nrd
    DOUBLE PRECISION, INTENT(OUT)            :: fsafe(lfsafe,nrd)
    INTEGER, INTENT(IN OUT)                  :: iout
    INTEGER, INTENT(OUT)                     :: ipt
    DOUBLE PRECISION, INTENT(OUT)            :: ysafe(km,nrd)
    DOUBLE PRECISION, INTENT(OUT)            :: ysafep(km,nrd)
    INTEGER, INTENT(IN OUT)                  :: icomp(nrd)
    REAL, INTENT(IN OUT)                     :: rpar(*)
    INTEGER, INTENT(IN OUT)                  :: ipar(*)
    INTEGER, INTENT(OUT)                     :: nfcn
    
    integer :: i, m, mm, l
    double precision :: hj, hj2, dblenj, expo, facmin

    hj=h/nj(j)
    hj2=hj*2.d0
    ! --- EULER STARTING STEP
    DO  i=1,n
        yh1(i)=y(i)
        yh2(i)=yp(i)+hj*dz(i)
    END DO
    ! --- EXPLICIT MIDPOINT (STOERMER) RULE
    m=nj(j)/2
    DO  mm=1,m
        IF (iout >= 2.AND.mm == (m+1)/2) THEN
            DO  i=1,nrd
                ysafep(j,i)=yh2(icomp(i))
                ysafe(j,i) =yh1(icomp(i))+hj*yh2(icomp(i))
            END DO
        END IF
        DO  i=1,n
            yh1(i)=yh1(i)+hj2*yh2(i)
        END DO
        CALL fcn(n,x+hj2*mm,yh1,dy,rpar,ipar)
        IF (iout >= 2.AND.ABS(mm-(m+1)/2) < j) THEN
            ipt=ipt+1
            DO  i=1,nrd
                fsafe(ipt,i)=dy(icomp(i))
            END DO
        END IF
        IF (mm == m) EXIT
        DO  i=1,n
            yh2(i)=yh2(i)+hj2*dy(i)
        END DO
    END DO
    ! --- FINAL SMOOTHING STEP
37  CONTINUE
    DO  i=1,n
        t(j,i)=yh1(i)
        tp(j,i)=yh2(i)+hj*dy(i)
    END DO
    nfcn=nfcn+m
    ! --- POLYNOMIAL EXTRAPOLATION
    IF (j == 1) RETURN
    dblenj=nj(j)
    DO  l=j,2,-1
        fac=(dblenj/nj(l-1))**2-1.d0
        DO  i=1,n
            t(l-1,i)=t(l,i)+(t(l,i)-t(l-1,i))/fac
            tp(l-1,i)=tp(l,i)+(tp(l,i)-tp(l-1,i))/fac
        END DO
    END DO
    ERR=0.d0
    ! --- SCALING
    IF (itol == 0) THEN
        DO  i=1,n
            scal(i)=atol(1)+rtol(1)*MAX(ABS(y(i)),ABS(t(1,i)))
            scalp(i)=atol(1)+rtol(1)*MAX(ABS(yp(i)),ABS(tp(1,i)))
        END DO
    ELSE
        DO  i=1,n
            scal(i)=atol(i)+rtol(i)*MAX(ABS(y(i)),ABS(t(1,i)))
            scalp(i)=atol(i+n)+rtol(i+n)*MAX(ABS(yp(i)),ABS(tp(1,i)))
        END DO
    END IF
    DO  i=1,n
        ERR=ERR+((t(1,i)-t(2,i))/scal(i))**2 +((tp(1,i)-tp(2,i))/scalp(i))**2
    END DO
    ERR=SQRT(ERR/(2*n))
    IF (ERR*uround >= 1.d0) GO TO 79
    IF (j > 2.AND.ERR >= errold) GO TO 79
    errold=MAX(4*ERR,1.d0)
    ! --- COMPUTE OPTIMAL STEPSIZES
    expo=1.d0/(2*j-1)
    facmin=fac1**expo
    fac=MIN(fac2/facmin,MAX(facmin,(ERR/safe1)**expo/safe2))
    fac=1.d0/fac
    hh(j)=MIN(ABS(h)*fac,hmax)
    w(j)=a(j)/hh(j)
    RETURN
79  atov=.true.
    h=h*safe3
    reject=.true.
    RETURN
    END SUBROUTINE stoerm

    SUBROUTINE intpo2(n,y,imit)
    ! --- COMPUTES THE COEFFICIENTS OF THE INTERPOLATION FORMULA
    IMPLICIT DOUBLE PRECISION (a-h,o-z)

    integer, parameter :: max_count = 40
    INTEGER, INTENT(IN OUT)                  :: n
    INTEGER, INTENT(IN)                      :: imit
    DOUBLE PRECISION, INTENT(IN OUT)         :: y(n*(imit+7))
    DOUBLE PRECISION :: a(0:max_count)
    
    integer :: i, im
    DOUBLE PRECISION :: y0, y1, yp0, yp1, ypp0, ypp1, ydiff
    DOUBLE PRECISION :: ah, bh, ch, dh, eh, fh, gh, abh, gfh, gmf
    DOUBLE PRECISION :: ph0, ph1, ph2, ph3, ph4, ph5, fc1, fc2, fc3
    
    ! --- BEGIN WITH HERMITE INTERPOLATION (ORDER 5)
    DO  i=1,n
        y0=y(i)
        y1=y(3*n+i)
        yp0=y(n+i)
        yp1=y(4*n+i)
        ypp0=y(2*n+i)/2.d0
        ypp1=y(5*n+i)/2.d0
        ydiff=y1-y0
        ah=ydiff-yp0
        bh=yp1-ydiff
        ch=ah-ypp0
        dh=bh-ah
        eh=ypp1-bh
        fh=dh-ch
        gh=eh-dh
        y(n+i)=ydiff
        y(2*n+i)=-ah
        y(3*n+i)=-dh
        y(4*n+i)=gh
        y(5*n+i)=fh
        IF (imit < 0) CYCLE
        ! --- COMPUTE THE DERIVATIVES OF HERMITE AT MIDPOINT
        abh=ah+bh
        gfh=gh+fh
        gmf=gh-fh
        ph0=0.5D0*(y0+y1+0.25D0*(-abh+0.25D0*gfh))
        ph1=ydiff+0.25D0*(ah-bh+0.25D0*gmf)
        ph2=abh-0.5D0*gfh
        ph3=6.d0*(bh-ah)-3.d0*gmf
        ph4=12.d0*gfh
        ph5=120.d0*gmf
        ! --- COMPUTE THE FURTHER COEFFICIENTS
        IF (imit < 1) GO TO 20
        a(1)=64.d0*(y(7*n+i)-ph1)
        IF (imit < 3) GO TO 20
        a(3)=64.d0*(y(9*n+i)-ph3+a(1)*9.d0/8.d0)
        IF (imit < 5) GO TO 20
        a(5)=64.d0*(y(11*n+i)-ph5+a(3)*15.d0/4.d0-a(1)*90.d0)
        IF (imit < 7) GO TO 20
        DO  im=7,imit,2
            fc1=im*(im-1)*3.d0/16.d0
            fc2=fc1*(im-2)*(im-3)*4.d0
            fc3=im*(im-1)*(im-2)*(im-3)*(im-4)*(im-5)
            a(im)=64.d0*(y((im+6)*n+i)+fc1*a(im-2)-fc2*a(im-4)+fc3*a(im-6))
        END DO
20      CONTINUE
        a(0)=64.d0*(y(6*n+i)-ph0)
        IF (imit < 2) GO TO 60
        a(2)=64.d0*(y(n*8+i)-ph2+a(0)*3.d0/8.d0)
        IF (imit < 4) GO TO 60
        a(4)=64.d0*(y(n*10+i)-ph4+a(2)*9.d0/4.d0-a(0)*18.d0)
        IF (imit < 6) GO TO 60
        DO  im=6,imit,2
            fc1=im*(im-1)*3.d0/16.d0
            fc2=fc1*(im-2)*(im-3)*4.d0
            fc3=im*(im-1)*(im-2)*(im-3)*(im-4)*(im-5)
            a(im)=64.d0*(y(n*(im+6)+i)+a(im-2)*fc1-a(im-4)*fc2+a(im-6)*fc3)
        END DO
60      CONTINUE
        DO  im=0,imit
            y(n*(im+6)+i)=a(im)
        END DO
    END DO
    RETURN
    END SUBROUTINE intpo2

    FUNCTION contx2(ii,x,y,ncon,icomp,n)
    ! ----------------------------------------------------------
    !     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT IN CONECTION
    !     WITH THE OUTPUT-SUBROUTINE FOR ODEX2. IT PROVIDES AN
    !     APPROXIMATION TO THE II-TH COMPONENT OF THE SOLUTION AT X.
    ! ----------------------------------------------------------

    IMPLICIT DOUBLE PRECISION (a-h,o-z)
    DOUBLE PRECISION                         :: contx2
    INTEGER, INTENT(IN)                      :: ii
    DOUBLE PRECISION, INTENT(IN OUT)         :: x
    INTEGER, INTENT(IN OUT)                  :: ncon
    INTEGER, INTENT(IN)                      :: n
    DOUBLE PRECISION, INTENT(IN)             :: y(ncon)
    INTEGER, INTENT(IN OUT)                  :: icomp(n)
    
    integer :: i, j, imit, im, iout
    DOUBLE PRECISION :: theta, theta1, phthet, xold, h, thetah

    COMMON /conod2/xold,h,imit
    ! ----- COMPUTE PLACE OF II-TH COMPONENT
    i=0
    DO  j=1,n
        IF (icomp(j) == ii) i=j
    END DO
    IF (i == 0) THEN
        WRITE (6,*) ' NO DENSE OUTPUT AVAILABLE FOR COMP.',ii
        RETURN
    END IF
    ! ----- COMPUTE THE INTERPOLATED VALUE
    theta=(x-xold)/h
    theta1=1.d0-theta
    phthet=y(i)+theta*(y(n+i)+theta1*(y(2*n+i)+theta*(y(3*n+i)+  &
        theta1*(y(4*n+i)*theta+y(5*n+i)*theta1))))
    IF (imit < 0) THEN
        contx2=phthet
        RETURN
    END IF
    thetah=theta-0.5D0
    contx2=y(n*(imit+6)+i)
    DO  im=imit,1,-1
        contx2=y(n*(im+5)+i)+contx2*thetah/im
    END DO
    contx2=phthet+(theta*theta1)**3*contx2
    RETURN
    END FUNCTION contx2
    
    subroutine driver_odex2()
    ! * * * * * * * * * * * * * * * * * * * * * * * * *
 
    ! Code converted using TO_F90 by Alan Miller
    ! Date: 2023-11-13  Time: 08:56:21
 
    !                 DRIVER FOR ODEX2
    ! * * * * * * * * * * * * * * * * * * * * * * * * *
    !ompile odex2
    !feh dr_odex2 odex2
    IMPLICIT real(r8) (a-h,o-z)
    INTEGER, PARAMETER :: &
        ndgl=2, km=9, nrdens=2,  &
        lwork=ndgl*(2*km+6)+5*km+20+(km*(2*km+5)+6)*nrdens, &
        liwork=km+20+nrdens
    
    
    ! --- DIMENSION OF THE SYSTEM
    integer, parameter :: n=2
    
    integer :: iout, i, j, iwork, itol, idid, ipar(n)
    DOUBLE PRECISION :: x, y(n), yp(n), xend, tol, h, work ,solout
    DIMENSION work(lwork),iwork(liwork)
    REAL :: rtol(n), atol(n), rpar(n)
    
    ! --- OUTPUT ROUTINE AND DENSE OUTPUT IS USED DURING INTEGRATION
    iout=2
    ! --- INITIAL VALUES
    x=0.0D0
    y(1)=0.5D0
    y(2)=0.0D0
    yp(1)=0.0D0
    yp(2)=SQRT(3.d0)
    ! --- ENDPOINT OF INTEGRATION
    xend=20.0D0
    ! --- REQUIRED (RELATIVE) TOLERANCE
    tol=1.0D-10
    itol=0
    rtol=tol
    atol=tol
    ! --- DEFAULT VALUES FOR PARAMETERS
    DO  i=1,10
      iwork(i)=0
      work(i)=0.d0
    END DO
    h=0.01D0
    ! --- IF DENSE OUTPUT IS REQUIRED
    iwork(8)=n
    ! --- CALL OF THE SUBROUTINE DOPRI5
    CALL odex2(n,twob,x,y,yp,xend,h, rtol,atol,itol,  &
        solout,iout, work,lwork,iwork,liwork,rpar,ipar,idid)
    ! --- PRINT FINAL SOLUTION
    WRITE (6,99) x,y(1),y(2)
    99     FORMAT(1X,'X =',f5.2,'    Y =',2E18.10)
    ! --- PRINT STATISTICS
    WRITE (6,90) tol
    90     FORMAT('       tol=',d12.2)
    WRITE (6,91) (iwork(j),j=17,20)
    91 FORMAT(' fcn=',i5,' step=',i4,' accpt=',i4,' rejct=',i3)

    end subroutine
    
    SUBROUTINE sol_out (nr,xold,x,y,yp,n,con,ncon,icomp,nd,rpar,ipar,irtrn)
    ! --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
    ! --- BY USING "CONTD5", THE CONTINUOUS COLLOCATION SOLUTION

    IMPLICIT real(r8) (a-h,o-z)
    INTEGER, INTENT(IN)                  :: nr
    INTEGER, INTENT(IN)                  :: n
    DOUBLE PRECISION, INTENT(IN)         :: xold
    DOUBLE PRECISION, INTENT(IN)         :: x
    DOUBLE PRECISION, INTENT(IN)         :: y(n)
    DOUBLE PRECISION, INTENT(IN)         :: yp(n)
    INTEGER, INTENT(IN OUT)              :: ncon
    INTEGER, INTENT(IN)                  :: nd
    DOUBLE PRECISION, INTENT(IN)         :: con(ncon)
    INTEGER, INTENT(IN OUT)              :: icomp(nd)
    REAL, INTENT(IN)                     :: rpar(*)
    INTEGER, INTENT(IN)                  :: ipar(*)
    INTEGER, INTENT(IN)                  :: irtrn
    
    integer :: i
    double precision :: xout, sol1, sol2

    COMMON /intern/xout

    IF (nr == 1) THEN
      WRITE (6,99) x,y(1),y(2),nr-1
      xout=x+1.0D0
    ELSE
    !           WRITE (6,99) X,Y(1),Y(2),NR-1
      10        CONTINUE
      IF (x >= xout) THEN
        sol1=contx2(1,xout,con,ncon,icomp,nd)
        sol2=contx2(2,xout,con,ncon,icomp,nd)
        WRITE (6,99) xout,sol1,sol2,nr-1
        xout=xout+1.0D0
        GO TO 10
      END IF
    END IF
    99     FORMAT(1X,'X =',f5.2,'    Y =',2E18.10,'    NSTEP =',i4)
    RETURN
    END SUBROUTINE sol_out
    
    ! CALL fcn(n,x,y,yp,rpar,ipar)
    SUBROUTINE twob(n,x,y,yp,rpar,ipar)
    ! --- RIGHT-HAND SIDE OF VAN DER POL'S EQUATION

    IMPLICIT real(r8) (a-h,o-z)
    INTEGER, INTENT(IN)                 :: n
    DOUBLE PRECISION, INTENT(IN)        :: x
    DOUBLE PRECISION, INTENT(IN)        :: y(n)
    DOUBLE PRECISION, INTENT(OUT)       :: yp(n)
    REAL, INTENT(IN OUT)                :: rpar(*)
    INTEGER, INTENT(IN OUT)             :: ipar(*)
    DOUBLE PRECISION :: rad

    rad=y(1)**2+y(2)**2
    rad=rad*SQRT(rad)
    yp(1)=-y(1)/rad
    yp(2)=-y(2)/rad
    RETURN
    END SUBROUTINE twob
    
    end module