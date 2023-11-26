    !  Reference:
    !
    !    Erwin Fehlberg,
    !    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
    !    NASA Technical Report R-315, 1969.
    !
    !    Lawrence Shampine, Herman Watts, S Davenport,
    !    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
    !    SIAM Review,
    !    Volume 18, pages 376-411, 1976.


    module mod_nasa_ode
    use mod_common
    implicit none
    
    interface
        pure subroutine r8_deriv_scalar(t,y,yp)
        import
        real ( kind = 8 ), intent(in) :: t
        real ( kind = 8 ), intent(in) :: y
        real ( kind = 8 ), intent(out):: yp
        end subroutine
        
        pure subroutine r4_deriv_vector(t,n,y,yp)
        import
        real   ( kind = 4 ), intent(in) :: t
        integer( kind = 4 ), intent(in) :: n
        real   ( kind = 4 ), intent(in) :: y(n)
        real   ( kind = 4 ), intent(out):: yp(n)
        end subroutine
        
        pure subroutine r8_deriv_vector(t,n,y,yp)
        import
        real   ( kind = 8 ), intent(in) :: t
        integer( kind = 4 ), intent(in) :: n
        real   ( kind = 8 ), intent(in) :: y(n)
        real   ( kind = 8 ), intent(out):: yp(n)
        end subroutine
    end interface

    contains
    

    pure subroutine rk4 ( t0, u0, dt, f, u )

    !*****************************************************************************80
    !
    !! RK4 takes one Runge-Kutta step for a scalar ODE.
    !
    !  Discussion:
    !
    !    It is assumed that an initial value problem, of the form
    !
    !      du/dt = f ( t, u )
    !      u(t0) = u0
    !
    !    is being solved.
    !
    !    If the user can supply current values of t, u, a stepsize dt, and a
    !    function to evaluate the derivative, this function can compute the
    !    fourth-order Runge Kutta estimate to the solution at time t+dt.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    31 January 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T0, the current time.
    !
    !    Input, real ( kind = 8 ) U0, the solution estimate at the current time.
    !
    !    Input, real ( kind = 8 ) DT, the time step.
    !
    !    Input, external F, a subroutine of the form
    !      subroutine f ( t, u, uprime )
    !    which evaluates the derivative uprime given the time T and
    !    solution vector U.
    !
    !    Output, real ( kind = 8 ) U, the fourth-order Runge-Kutta solution
    !    estimate at time T0+DT.
    !
    implicit none
    
    real ( kind = 8 ) , intent(in) :: t0
    real ( kind = 8 ) , intent(in) :: u0
    real ( kind = 8 ) , intent(in) :: dt
    procedure(r8_deriv_scalar), pointer, intent(in) :: f
    real ( kind = 8 ) , intent(out) :: u
    
    real ( kind = 8 ) f0
    real ( kind = 8 ) f1
    real ( kind = 8 ) f2
    real ( kind = 8 ) f3
    real ( kind = 8 ) t1
    real ( kind = 8 ) t2
    real ( kind = 8 ) t3
    real ( kind = 8 ) u1
    real ( kind = 8 ) u2
    real ( kind = 8 ) u3
    !
    !  Get four sample values of the derivative.
    !
    call f ( t0, u0, f0 )

    t1 = t0 + dt / 2.0D+00
    u1 = u0 + dt * f0 / 2.0D+00
    call f ( t1, u1, f1 )

    t2 = t0 + dt / 2.0D+00
    u2 = u0 + dt * f1 / 2.0D+00
    call f ( t2, u2, f2 )

    t3 = t0 + dt
    u3 = u0 + dt * f2
    call f ( t3, u3, f3 )
    !
    !  Combine them to estimate the solution U at time T1.
    !
    u = u0 + dt * ( f0 + 2.0D+00 * f1 + 2.0D+00 * f2 + f3 ) / 6.0D+00

    return
    end
    pure subroutine rk4vec ( t0, m, u0, dt, f, u )

    !*****************************************************************************80
    !
    !! RK4VEC takes one Runge-Kutta step for a vector ODE.
    !
    !  Discussion:
    !
    !    Thanks  to Dante Bolatti for correcting the final function call to:
    !      call f ( t3, m, u3, f3 )
    !    18 August 2016.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 August 2016
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T0, the current time.
    !
    !    Input, integer ( kind = 4 ) M, the dimension of the system.
    !
    !    Input, real ( kind = 8 ) U0(M), the solution estimate at the current time.
    !
    !    Input, real ( kind = 8 ) DT, the time step.
    !
    !    Input, external F, a subroutine of the form
    !      subroutine f ( t, m, u, uprime )
    !    which evaluates the derivative UPRIME(1:M) given the time T and
    !    solution vector U(1:M).
    !
    !    Output, real ( kind = 8 ) U(M), the fourth-order Runge-Kutta solution
    !    estimate at time T0+DT.
    !
    implicit none
    !subroutine rk4vec ( t0, m, u0, dt, f, u )
    real ( kind = 8 ) , intent(in) :: t0
    integer ( kind = 4 ) , intent(in) :: m
    real ( kind = 8 ) , intent(in) :: u0(m)
    real ( kind = 8 ) , intent(in) :: dt
    procedure(r8_deriv_vector), pointer, intent(in) :: f
    real ( kind = 8 ) , intent(out) :: u(m)

    
    real ( kind = 8 ) f0(m)
    real ( kind = 8 ) f1(m)
    real ( kind = 8 ) f2(m)
    real ( kind = 8 ) f3(m)
    real ( kind = 8 ) t1
    real ( kind = 8 ) t2
    real ( kind = 8 ) t3
    real ( kind = 8 ) u1(m)
    real ( kind = 8 ) u2(m)
    real ( kind = 8 ) u3(m)
    !
    !  Get four sample values of the derivative.
    !
    call f ( t0, m, u0, f0 )

    t1 = t0 + dt / 2.0D+00
    u1(1:m) = u0(1:m) + dt * f0(1:m) / 2.0D+00
    call f ( t1, m, u1, f1 )

    t2 = t0 + dt / 2.0D+00
    u2(1:m) = u0(1:m) + dt * f1(1:m) / 2.0D+00
    call f ( t2, m, u2, f2 )

    t3 = t0 + dt
    u3(1:m) = u0(1:m) + dt * f2(1:m)
    call f ( t3, m, u3, f3 )
    !
    !  Combine them to estimate the solution U at time T1.
    !
    u(1:m) = u0(1:m) + ( dt / 6.0D+00 ) * ( &
        f0(1:m) &
        + 2.0D+00 * f1(1:m) &
        + 2.0D+00 * f2(1:m) &
        +           f3(1:m) )

    return
    end

    pure subroutine r4_fehl ( f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, s )

    !*****************************************************************************80
    !
    !! R4_FEHL takes one Fehlberg fourth-fifth order step (single precision).
    !
    !  Discussion:
    !
    !    This routine integrates a system of NEQN first order ordinary differential
    !    equations of the form
    !      dY(i)/dT = F(T,Y(1:NEQN))
    !    where the initial values Y and the initial derivatives
    !    YP are specified at the starting point T.
    !
    !    The routine advances the solution over the fixed step H and returns
    !    the fifth order (sixth order accurate locally) solution
    !    approximation at T+H in array S.
    !
    !    The formulas have been grouped to control loss of significance.
    !    The routine should be called with an H not smaller than 13 units of
    !    roundoff in T so that the various independent arguments can be
    !    distinguished.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    27 March 2004
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Erwin Fehlberg,
    !    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
    !    NASA Technical Report R-315, 1969.
    !
    !    Lawrence Shampine, Herman Watts, S Davenport,
    !    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
    !    SIAM Review,
    !    Volume 18, pages 376-411, 1976.
    !
    !  Parameters:
    !
    !    Input, external F, a user-supplied subroutine to evaluate the
    !    derivatives Y'(T), of the form:
    !
    !      subroutine f ( t, y, yp )
    !      real ( kind = 4 ) t
    !      real ( kind = 4 ) y(*)
    !      real ( kind = 4 ) yp(*)
    !
    !    Input, integer ( kind = 4 ) NEQN, the number of equations to be integrated.
    !
    !    Input, real ( kind = 4 ) Y(NEQN), the current value of the
    !    dependent variable.
    !
    !    Input, real ( kind = 4 ) T, the current value of the independent
    !    variable.
    !
    !    Input, real ( kind = 4 ) H, the step size to take.
    !
    !    Input, real ( kind = 4 ) YP(NEQN), the current value of the
    !    derivative of the dependent variable.
    !
    !    Output, real ( kind = 4 ) F1(NEQN), F2(NEQN), F3(NEQN), F4(NEQN),
    !    F5(NEQN), derivative values needed for the computation.
    !
    !    Output, real ( kind = 4 ) S(NEQN), the estimate of the solution at T+H.
    !
    implicit none
    ! subroutine r4_fehl ( f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, s )
    procedure(r4_deriv_vector), pointer, intent(in) :: f
    integer ( kind = 4 ) , intent(in) :: neqn
    real ( kind = 4 ) , intent(in) :: y(neqn)
    real ( kind = 4 ) , intent(in) :: t
    real ( kind = 4 ) , intent(in) :: h
    real ( kind = 4 ) , intent(in) :: yp(neqn)
    real ( kind = 4 ) , intent(out) :: f1(neqn)
    real ( kind = 4 ) , intent(out) :: f2(neqn)
    real ( kind = 4 ) , intent(out) :: f3(neqn)
    real ( kind = 4 ) , intent(out) :: f4(neqn)
    real ( kind = 4 ) , intent(out) :: f5(neqn)
    real ( kind = 4 ) , intent(out) :: s(neqn)

    real ( kind = 4 ) ch
    !external f

    ch = h / 4.0E+00

    f5(1:neqn) = y(1:neqn) + ch * yp(1:neqn)

    call f ( t + ch, neqn, f5, f1 )

    ch = 3.0E+00 * h / 32.0E+00

    f5(1:neqn) = y(1:neqn) + ch * ( yp(1:neqn) + 3.0E+00 * f1(1:neqn) )

    call f ( t + 3.0E+00 * h / 8.0E+00, neqn, f5, f2 )

    ch = h / 2197.0E+00

    f5(1:neqn) = y(1:neqn) + ch * &
        ( 1932.0E+00 * yp(1:neqn) &
        + ( 7296.0E+00 * f2(1:neqn) - 7200.0E+00 * f1(1:neqn) ) &
        )

    call f ( t + 12.0E+00 * h / 13.0E+00, neqn, f5, f3 )

    ch = h / 4104.0E+00

    f5(1:neqn) = y(1:neqn) + ch * &
        ( &
        ( 8341.0E+00 * yp(1:neqn) - 845.0E+00 * f3(1:neqn) ) &
        + ( 29440.0E+00 * f2(1:neqn) - 32832.0E+00 * f1(1:neqn) ) &
        )

    call f ( t + h, neqn, f5, f4 )

    ch = h / 20520.0E+00

    f1(1:neqn) = y(1:neqn) + ch * &
        ( &
        ( -6080.0E+00 * yp(1:neqn) &
        + ( 9295.0E+00 * f3(1:neqn) - 5643.0E+00 * f4(1:neqn) ) &
        ) &
        + ( 41040.0E+00 * f1(1:neqn) - 28352.0E+00 * f2(1:neqn) ) &
        )

    call f ( t + h / 2.0E+00, neqn, f1, f5 )
    !
    !  Ready to compute the approximate solution at T+H.
    !
    ch = h / 7618050.0E+00

    s(1:neqn) = y(1:neqn) + ch * &
        ( &
        ( 902880.0E+00 * yp(1:neqn) &
        + ( 3855735.0E+00 * f3(1:neqn) - 1371249.0E+00 * f4(1:neqn) ) ) &
        + ( 3953664.0E+00 * f2(1:neqn) + 277020.0E+00 * f5(1:neqn) ) &
        )

    return
    end
    pure subroutine r4_rkf45 ( f, neqn, y, yp, t, tout, relerr, abserr, flag, h)

    !*****************************************************************************80
    !
    !! R4_RKF45 carries out the Runge-Kutta-Fehlberg method (single precision).
    !
    !  Discussion:
    !
    !    This routine is primarily designed to solve non-stiff and mildly stiff
    !    differential equations when derivative evaluations are inexpensive.
    !    It should generally not be used when the user is demanding
    !    high accuracy.
    !
    !    This routine integrates a system of NEQN first-order ordinary differential
    !    equations of the form:
    !
    !      dY(i)/dT = F(T,Y(1),Y(2),...,Y(NEQN))
    !
    !    where the Y(1:NEQN) are given at T.
    !
    !    Typically the subroutine is used to integrate from T to TOUT but it
    !    can be used as a one-step integrator to advance the solution a
    !    single step in the direction of TOUT.  On return, the parameters in
    !    the call list are set for continuing the integration.  The user has
    !    only to call again (and perhaps define a new value for TOUT).
    !
    !    Before the first call, the user must
    !
    !    * supply the subroutine F(T,Y,YP) to evaluate the right hand side;
    !      and declare F in an EXTERNAL statement;
    !
    !    * initialize the parameters:
    !      NEQN, Y(1:NEQN), T, TOUT, RELERR, ABSERR, FLAG.
    !      In particular, T should initially be the starting point for integration,
    !      Y should be the value of the initial conditions, and FLAG should
    !      normally be +1.
    !
    !    Normally, the user only sets the value of FLAG before the first call, and
    !    thereafter, the program manages the value.  On the first call, FLAG should
    !    normally be +1 (or -1 for single step mode.)  On normal return, FLAG will
    !    have been reset by the program to the value of 2 (or -2 in single
    !    step mode), and the user can continue to call the routine with that
    !    value of FLAG.
    !
    !    (When the input magnitude of FLAG is 1, this indicates to the program
    !    that it is necessary to do some initialization work.  An input magnitude
    !    of 2 lets the program know that that initialization can be skipped,
    !    and that useful information was computed earlier.)
    !
    !    The routine returns with all the information needed to continue
    !    the integration.  If the integration reached TOUT, the user need only
    !    define a new TOUT and call again.  In the one-step integrator
    !    mode, returning with FLAG = -2, the user must keep in mind that
    !    each step taken is in the direction of the current TOUT.  Upon
    !    reaching TOUT, indicated by the output value of FLAG switching to 2,
    !    the user must define a new TOUT and reset FLAG to -2 to continue
    !    in the one-step integrator mode.
    !
    !    In some cases, an error or difficulty occurs during a call.  In that case,
    !    the output value of FLAG is used to indicate that there is a problem
    !    that the user must address.  These values include:
    !
    !    * 3, integration was not completed because the input value of RELERR, the
    !      relative error tolerance, was too small.  RELERR has been increased
    !      appropriately for continuing.  If the user accepts the output value of
    !      RELERR, then simply reset FLAG to 2 and continue.
    !
    !    * 4, integration was not completed because more than MAXNFE derivative
    !      evaluations were needed.  This is approximately (MAXNFE/6) steps.
    !      The user may continue by simply calling again.  The function counter
    !      will be reset to 0, and another MAXNFE function evaluations are allowed.
    !
    !    * 5, integration was not completed because the solution vanished,
    !      making a pure relative error test impossible.  The user must use
    !      a non-zero ABSERR to continue.  Using the one-step integration mode
    !      for one step is a good way to proceed.
    !
    !    * 6, integration was not completed because the requested accuracy
    !      could not be achieved, even using the smallest allowable stepsize.
    !      The user must increase the error tolerances ABSERR or RELERR before
    !      continuing.  It is also necessary to reset FLAG to 2 (or -2 when
    !      the one-step integration mode is being used).  The occurrence of
    !      FLAG = 6 indicates a trouble spot.  The solution is changing
    !      rapidly, or a singularity may be present.  It often is inadvisable
    !      to continue.
    !
    !    * 7, it is likely that this routine is inefficient for solving
    !      this problem.  Too much output is restricting the natural stepsize
    !      choice.  The user should use the one-step integration mode with
    !      the stepsize determined by the code.  If the user insists upon
    !      continuing the integration, reset FLAG to 2 before calling
    !      again.  Otherwise, execution will be terminated.
    !
    !    * 8, invalid input parameters, indicates one of the following:
    !      NEQN <= 0;
    !      T = TOUT and |FLAG| /= 1;
    !      RELERR < 0 or ABSERR < 0;
    !      FLAG == 0  or FLAG < -2 or 8 < FLAG.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 April 2011
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Erwin Fehlberg,
    !    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
    !    NASA Technical Report R-315, 1969.
    !
    !    Lawrence Shampine, Herman Watts, S Davenport,
    !    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
    !    SIAM Review,
    !    Volume 18, pages 376-411, 1976.
    !
    !  Parameters:
    !
    !    Input, external F, a user-supplied subroutine to evaluate the
    !    derivatives Y'(T), of the form:
    !
    !      subroutine f ( t, y, yp )
    !      real ( kind = 4 ) t
    !      real ( kind = 4 ) y(*)
    !      real ( kind = 4 ) yp(*)
    !
    !    Input, integer ( kind = 4 ) NEQN, the number of equations to be integrated.
    !
    !    Input/output, real ( kind = 4 ) Y(NEQN), the current solution vector at T.
    !
    !    Input/output, real ( kind = 4 ) YP(NEQN), the current value of the
    !    derivative of the dependent variable.  The user should not set or alter
    !    this information!
    !
    !    Input/output, real ( kind = 4 ) T, the current value of the independent
    !    variable.
    !
    !    Input, real ( kind = 4 ) TOUT, the output point at which solution is
    !    desired.  TOUT = T is allowed on the first call only, in which case
    !    the routine returns with FLAG = 2 if continuation is possible.
    !
    !    Input, real ( kind = 4 ) RELERR, ABSERR, the relative and absolute
    !    error tolerances for the local error test.  At each step the code
    !    requires:
    !      abs ( local error ) <= RELERR * abs ( Y ) + ABSERR
    !    for each component of the local error and the solution vector Y.
    !    RELERR cannot be "too small".  If the routine believes RELERR has been
    !    set too small, it will reset RELERR to an acceptable value and return
    !    immediately for user action.
    !
    !    Input/output, integer ( kind = 4 ) FLAG, indicator for status of
    !    integration.  On the first call, set FLAG to +1 for normal use, or to
    !    -1 for single step mode.  On return, a value of 2 or -2 indicates normal
    !    progress, while any other value indicates a problem that should be
    !    addressed.
    !
    implicit none
    
    procedure(r4_deriv_vector), pointer, intent(in) :: f
    integer ( kind = 4 ) , intent(in) :: neqn
    real ( kind = 4 ) , intent(inout) :: y(neqn)
    real ( kind = 4 ) , intent(inout) :: yp(neqn)
    real ( kind = 4 ) , intent(inout) :: t
    real ( kind = 4 ) , intent(in) :: tout
    real ( kind = 4 ) , intent(inout) :: relerr
    real ( kind = 4 ) , intent(in) :: abserr
    integer ( kind = 4 ), intent(out) :: flag
    real ( kind = 4 ), intent(out) :: h

    real ( kind = 4 ) ae
    real ( kind = 4 ) dt
    real ( kind = 4 ) ee
    real ( kind = 4 ) eeoet
    real ( kind = 4 ) eps
    real ( kind = 4 ) esttol
    real ( kind = 4 ) et
    
    real ( kind = 4 ) f1(neqn)
    real ( kind = 4 ) f2(neqn)
    real ( kind = 4 ) f3(neqn)
    real ( kind = 4 ) f4(neqn)
    real ( kind = 4 ) f5(neqn)    
    real ( kind = 4 ) s
    real ( kind = 4 ) scale
    real ( kind = 4 ) tol
    real ( kind = 4 ) toln
    real ( kind = 4 ) ypk
    
    integer ( kind = 4 ), parameter :: maxnfe = 3000
    real ( kind = 4 ), parameter :: remin = 1.0E-12
    
    logical hfaild
    real ( kind = 4 ) hmin
    integer ( kind = 4 ) k
    integer ( kind = 4 ) mflag
    logical output
    real ( kind = 4 ) relerr_min
    real ( kind = 4 ) :: relerr_save 
    real ( kind = 4 ) :: abserr_save 
    integer ( kind = 4 ) :: flag_save 
    integer ( kind = 4 ) :: init 
    integer ( kind = 4 ) :: kflag 
    integer ( kind = 4 ) :: kop 
    integer ( kind = 4 ) :: nfe 
    !
    ! Default Vaalues
    !
    h = -1.0E+00
    relerr_save = -1.0E+00
    abserr_save = -1.0E+00
    flag_save = -1000
    init = -1000
    kflag = -1000
    kop = -1
    nfe = -1
    
    !
    !  Check the input parameters.
    !
    eps = epsilon ( eps )

    if ( neqn < 1 ) then
        flag = 8
        return
    end if

    if ( relerr < 0.0E+00 ) then
        flag = 8
        return
    end if

    if ( abserr < 0.0E+00 ) then
        flag = 8
        return
    end if

    if ( flag == 0 .or. 8 < flag .or. flag < -2 ) then
        flag = 8
        return
    end if

    mflag = abs ( flag )
    !
    !  Is this a continuation call?
    !
    if ( mflag /= 1 ) then

        if ( t == tout .and. kflag /= 3 ) then
            flag = 8
            return
        end if

        if ( mflag == 2 ) then

            if ( kflag == 3 ) then

                flag = flag_save
                mflag = abs ( flag )

            else if ( init == 0 ) then

                flag = flag_save

            else if ( kflag == 4 ) then

                nfe = 0

            else if ( kflag == 5 .and. abserr == 0.0E+00 ) then

                error stop ' '//NEW_LINE('a') & 
                    //'R4_RKF45 - Fatal error!'//NEW_LINE('a') & 
                    //'  KFLAG = 5 and ABSERR = 0.0'//NEW_LINE('a') 
                
                !write ( *, '(a)' ) ' '
                !write ( *, '(a)' ) 'R4_RKF45 - Fatal error!'
                !write ( *, '(a)' ) '  KFLAG = 5 and ABSERR = 0.0'
                !stop

            else if ( &
                kflag == 6 .and. &
                relerr <= relerr_save .and. &
                abserr <= abserr_save ) then

                error stop ' '//NEW_LINE('a') & 
                    //'R4_RKF45 - Fatal error!'//NEW_LINE('a') & 
                    //'  KFLAG = 6 and'//NEW_LINE('a') & 
                    //'  RELERR <= RELERR_SAVE and'//NEW_LINE('a') & 
                    //'  ABSERR <= ABSERR_SAVE'//NEW_LINE('a') 
                
                !write ( *, '(a)' ) ' '
                !write ( *, '(a)' ) 'R4_RKF45 - Fatal error!'
                !write ( *, '(a)' ) '  KFLAG = 6 and'
                !write ( *, '(a)' ) '  RELERR <= RELERR_SAVE and'
                !write ( *, '(a)' ) '  ABSERR <= ABSERR_SAVE'
                !stop

            end if
            !
            !  FLAG = 3, 4, 5, 6, 7 or 8.
            !
        else

            if ( flag == 3 ) then

                flag = flag_save
                if ( kflag == 3 ) then
                    mflag = abs ( flag )
                end if

            else if ( flag == 4 ) then

                nfe = 0
                flag = flag_save
                if ( kflag == 3 ) then
                    mflag = abs ( flag )
                end if

            else if ( flag == 5 .and. 0.0E+00 < abserr ) then

                flag = flag_save
                if ( kflag == 3 ) then
                    mflag = abs ( flag )
                end if
                !
                !  Integration cannot be continued because the user did not respond to
                !  the instructions pertaining to FLAG = 5, 6, 7 or 8.
                !
            else
                
                error stop ' '//NEW_LINE('a') & 
                    //'R4_RKF45 - Fatal error!'//NEW_LINE('a') & 
                    //'  Integration cannot be continued.'//NEW_LINE('a') & 
                    //'  The user did not respond to the output'//NEW_LINE('a') & 
                    //'  value FLAG = 5, 6, 7, or 8.'//NEW_LINE('a') 
                
                !write ( *, '(a)' ) ' '
                !write ( *, '(a)' ) 'R4_RKF45 - Fatal error!'
                !write ( *, '(a)' ) '  Integration cannot be continued.'
                !write ( *, '(a)' ) '  The user did not respond to the output'
                !write ( *, '(a)' ) '  value FLAG = 5, 6, 7, or 8.'
                !stop
            end if

        end if

    end if
    !
    !  Save the input value of FLAG.
    !  Set the continuation flag KFLAG for subsequent input checking.
    !
    flag_save = flag
    kflag = 0
    !
    !  Save RELERR and ABSERR for checking input on subsequent calls.
    !
    relerr_save = relerr
    abserr_save = abserr
    !
    !  Restrict the relative error tolerance to be at least
    !
    !    2 * EPS + REMIN
    !
    !  to avoid limiting precision difficulties arising from impossible
    !  accuracy requests.
    !
    relerr_min = 2.0E+00 * epsilon ( relerr_min ) + remin
    !
    !  Is the relative error tolerance too small?
    !
    if ( relerr < relerr_min ) then
        relerr = relerr_min
        flag = 3
        kflag = 3
        return
    end if

    dt = tout - t
    !
    !  Initialization:
    !
    !  Set the initialization completion indicator, INIT;
    !  set the indicator for too many output points, KOP;
    !  evaluate the initial derivatives;
    !  set the counter for function evaluations, NFE;
    !  estimate the starting stepsize.
    !
    if ( mflag == 1 ) then

        init = 0
        kop = 0
        call f ( t, neqn, y, yp )
        nfe = 1

        if ( t == tout ) then
            flag = 2
            return
        end if

    end if

    if ( init == 0 ) then

        init = 1
        h = abs ( dt )
        toln = 0.0E+00

        do k = 1, neqn
            tol = relerr * abs ( y(k) ) + abserr
            if ( 0.0E+00 < tol ) then
                toln = tol
                ypk = abs ( yp(k) )
                if ( tol < ypk * h**5 ) then
                    h = ( tol / ypk )**0.2E+00
                end if
            end if
        end do

        if ( toln <= 0.0E+00 ) then
            h = 0.0E+00
        end if

        h = max ( h, 26.0E+00 * eps * max ( abs ( t ), abs ( dt ) ) )
        flag_save = sign ( 2, flag )

    end if
    !
    !  Set the stepsize for integration in the direction from T to TOUT.
    !
    h = sign ( h, dt )
    !
    !  Test to see if too may output points are being requested.
    !
    if ( 2.0E+00 * abs ( dt ) <= abs ( h ) ) then
        kop = kop + 1
    end if
    !
    !  Unnecessary frequency of output.
    !
    if ( kop == 100 ) then
        kop = 0
        flag = 7
        return
    end if
    !
    !  If we are too close to the output point, then simply extrapolate and return.
    !
    if ( abs ( dt ) <= 26.0E+00 * eps * abs ( t ) ) then
        t = tout
        y(1:neqn) = y(1:neqn) + dt * yp(1:neqn)
        call f ( t, neqn, y, yp )
        nfe = nfe + 1
        flag = 2
        return
    end if
    !
    !  Initialize the output point indicator.
    !
    output = .false.
    !
    !  To avoid premature underflow in the error tolerance function,
    !  scale the error tolerances.
    !
    scale = 2.0E+00 / relerr
    ae = scale * abserr
    !
    !  Step by step integration.
    !
    do

        hfaild = .false.
        !
        !  Set the smallest allowable stepsize.
        !
        hmin = 26.0E+00 * eps * abs ( t )
        !
        !  Adjust the stepsize if necessary to hit the output point.
        !
        !  Look ahead two steps to avoid drastic changes in the stepsize and
        !  thus lessen the impact of output points on the code.
        !
        dt = tout - t

        if ( 2.0E+00 * abs ( h ) <= abs ( dt ) ) then

        else
            !
            !  Will the next successful step complete the integration to the output point?
            !
            if ( abs ( dt ) <= abs ( h ) ) then
                output = .true.
                h = dt
            else
                h = 0.5E+00 * dt
            end if

        end if
        !
        !  Here begins the core integrator for taking a single step.
        !
        !  The tolerances have been scaled to avoid premature underflow in
        !  computing the error tolerance function ET.
        !  To avoid problems with zero crossings, relative error is measured
        !  using the average of the magnitudes of the solution at the
        !  beginning and end of a step.
        !  The error estimate formula has been grouped to control loss of
        !  significance.
        !
        !  To distinguish the various arguments, H is not permitted
        !  to become smaller than 26 units of roundoff in T.
        !  Practical limits on the change in the stepsize are enforced to
        !  smooth the stepsize selection process and to avoid excessive
        !  chattering on problems having discontinuities.
        !  To prevent unnecessary failures, the code uses 9/10 the stepsize
        !  it estimates will succeed.
        !
        !  After a step failure, the stepsize is not allowed to increase for
        !  the next attempted step.  This makes the code more efficient on
        !  problems having discontinuities and more effective in general
        !  since local extrapolation is being used and extra caution seems
        !  warranted.
        !
        !  Test the number of derivative function evaluations.
        !  If okay, try to advance the integration from T to T+H.
        !
        do
            !
            !  Have we done too much work?
            !
            if ( maxnfe < nfe ) then
                flag = 4
                kflag = 4
                return
            end if
            !
            !  Advance an approximate solution over one step of length H.
            !
            call r4_fehl ( f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, f1 )
            nfe = nfe + 5
            !
            !  Compute and test allowable tolerances versus local error estimates
            !  and remove scaling of tolerances.  The relative error is
            !  measured with respect to the average of the magnitudes of the
            !  solution at the beginning and end of the step.
            !
            eeoet = 0.0E+00

            do k = 1, neqn

                et = abs ( y(k) ) + abs ( f1(k) ) + ae

                if ( et <= 0.0E+00 ) then
                    flag = 5
                    return
                end if

                ee = abs &
                    ( ( -2090.0E+00 * yp(k) &
                    + ( 21970.0E+00 * f3(k) - 15048.0E+00 * f4(k) ) &
                    ) &
                    + ( 22528.0E+00 * f2(k) - 27360.0E+00 * f5(k) ) &
                    )

                eeoet = max ( eeoet, ee / et )

            end do

            esttol = abs ( h ) * eeoet * scale / 752400.0E+00

            if ( esttol <= 1.0E+00 ) then
                exit
            end if
            !
            !  Unsuccessful step.  Reduce the stepsize, try again.
            !  The decrease is limited to a factor of 1/10.
            !
            hfaild = .true.
            output = .false.

            if ( esttol < 59049.0E+00 ) then
                s = 0.9E+00 / esttol**0.2E+00
            else
                s = 0.1E+00
            end if

            h = s * h

            if ( abs ( h ) < hmin ) then
                flag = 6
                kflag = 6
                return
            end if

        end do
        !
        !  We exited the loop because we took a successful step.
        !  Store the solution for T+H, and evaluate the derivative there.
        !
        t = t + h
        y(1:neqn) = f1(1:neqn)
        call f ( t, neqn, y, yp )
        nfe = nfe + 1
        !
        !  Choose the next stepsize.  The increase is limited to a factor of 5.
        !  If the step failed, the next stepsize is not allowed to increase.
        !
        if ( 0.0001889568E+00 < esttol ) then
            s = 0.9E+00 / esttol**0.2E+00
        else
            s = 5.0E+00
        end if

        if ( hfaild ) then
            s = min ( s, 1.0E+00 )
        end if

        h = sign ( max ( s * abs ( h ), hmin ), h )
        !
        !  End of core integrator
        !
        !  Should we take another step?
        !
        if ( output ) then
            t = tout
            flag = 2
            return
        end if

        if ( flag <= 0 ) then
            exit
        end if

    end do
    !
    !  One step integration mode.
    !
    flag = -2

    return
    end
    pure subroutine r8_fehl ( f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, s )

    !*****************************************************************************80
    !
    !! R8_FEHL takes one Fehlberg fourth-fifth order step (double precision).
    !
    !  Discussion:
    !
    !    This routine integrates a system of NEQN first order ordinary differential
    !    equations of the form
    !      dY(i)/dT = F(T,Y(1:NEQN))
    !    where the initial values Y and the initial derivatives
    !    YP are specified at the starting point T.
    !
    !    The routine advances the solution over the fixed step H and returns
    !    the fifth order (sixth order accurate locally) solution
    !    approximation at T+H in array S.
    !
    !    The formulas have been grouped to control loss of significance.
    !    The routine should be called with an H not smaller than 13 units of
    !    roundoff in T so that the various independent arguments can be
    !    distinguished.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    27 March 2004
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Erwin Fehlberg,
    !    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
    !    NASA Technical Report R-315, 1969.
    !
    !    Lawrence Shampine, Herman Watts, S Davenport,
    !    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
    !    SIAM Review,
    !    Volume 18, pages 376-411, 1976.
    !
    !  Parameters:
    !
    !    Input, external F, a user-supplied subroutine to evaluate the
    !    derivatives Y'(T), of the form:
    !
    !      subroutine f ( t, y, yp )
    !      real ( kind = 8 ) t
    !      real ( kind = 8 ) y(*)
    !      real ( kind = 8 ) yp(*)
    !
    !    Input, integer ( kind = 4 ) NEQN, the number of equations to be integrated.
    !
    !    Input, real ( kind = 8 ) Y(NEQN), the current value of the
    !    dependent variable.
    !
    !    Input, real ( kind = 8 ) T, the current value of the independent
    !    variable.
    !
    !    Input, real ( kind = 8 ) H, the step size to take.
    !
    !    Input, real ( kind = 8 ) YP(NEQN), the current value of the
    !    derivative of the dependent variable.
    !
    !    Output, real ( kind = 8 ) F1(NEQN), F2(NEQN), F3(NEQN), F4(NEQN),
    !    F5(NEQN), derivative values needed for the computation.
    !
    !    Output, real ( kind = 8 ) S(NEQN), the estimate of the solution at T+H.
    !
    implicit none
    ! subroutine r8_fehl ( f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, s )
    procedure(r8_deriv_vector), pointer, intent(in) :: f
    integer ( kind = 4 ), intent(in) :: neqn
    real ( kind = 8 ) , intent(in) :: y(neqn)
    real ( kind = 8 ) , intent(in) :: t
    real ( kind = 8 ) , intent(in) ::  h
    real ( kind = 8 ) , intent(in) :: yp(neqn)

    real ( kind = 8 ) ch
    !external f
    real ( kind = 8 ) , intent(out) :: f1(neqn)
    real ( kind = 8 ) , intent(out) :: f2(neqn)
    real ( kind = 8 ) , intent(out) :: f3(neqn)
    real ( kind = 8 ) , intent(out) :: f4(neqn)
    real ( kind = 8 ) , intent(out) :: f5(neqn)
    real ( kind = 8 ) , intent(out) :: s(neqn)

    ch = h / 4.0D+00

    f5(1:neqn) = y(1:neqn) + ch * yp(1:neqn)

    call f ( t + ch, neqn, f5, f1 )

    ch = 3.0D+00 * h / 32.0D+00

    f5(1:neqn) = y(1:neqn) + ch * ( yp(1:neqn) + 3.0D+00 * f1(1:neqn) )

    call f ( t + 3.0D+00 * h / 8.0D+00, neqn, f5, f2 )

    ch = h / 2197.0D+00

    f5(1:neqn) = y(1:neqn) + ch * &
        ( 1932.0D+00 * yp(1:neqn) &
        + ( 7296.0D+00 * f2(1:neqn) - 7200.0D+00 * f1(1:neqn) ) &
        )

    call f ( t + 12.0D+00 * h / 13.0D+00, neqn, f5, f3 )

    ch = h / 4104.0D+00

    f5(1:neqn) = y(1:neqn) + ch * &
        ( &
        ( 8341.0D+00 * yp(1:neqn) - 845.0D+00 * f3(1:neqn) ) &
        + ( 29440.0D+00 * f2(1:neqn) - 32832.0D+00 * f1(1:neqn) ) &
        )

    call f ( t + h, neqn, f5, f4 )

    ch = h / 20520.0D+00

    f1(1:neqn) = y(1:neqn) + ch * &
        ( &
        ( -6080.0D+00 * yp(1:neqn) &
        + ( 9295.0D+00 * f3(1:neqn) - 5643.0D+00 * f4(1:neqn) ) &
        ) &
        + ( 41040.0D+00 * f1(1:neqn) - 28352.0D+00 * f2(1:neqn) ) &
        )

    call f ( t + h / 2.0D+00, neqn, f1, f5 )
    !
    !  Ready to compute the approximate solution at T+H.
    !
    ch = h / 7618050.0D+00

    s(1:neqn) = y(1:neqn) + ch * &
        ( &
        ( 902880.0D+00 * yp(1:neqn) &
        + ( 3855735.0D+00 * f3(1:neqn) - 1371249.0D+00 * f4(1:neqn) ) ) &
        + ( 3953664.0D+00 * f2(1:neqn) + 277020.0D+00 * f5(1:neqn) ) &
        )

    return
    end
    pure subroutine r8_rkf45 ( f, neqn, y, yp, t, tout, relerr, abserr, flag, h )

    !*****************************************************************************80
    !
    !! R8_RKF45 carries out the Runge-Kutta-Fehlberg method (double precision).
    !
    !  Discussion:
    !
    !    This routine is primarily designed to solve non-stiff and mildly stiff
    !    differential equations when derivative evaluations are inexpensive.
    !    It should generally not be used when the user is demanding
    !    high accuracy.
    !
    !    This routine integrates a system of NEQN first-order ordinary differential
    !    equations of the form:
    !
    !      dY(i)/dT = F(T,Y(1),Y(2),...,Y(NEQN))
    !
    !    where the Y(1:NEQN) are given at T.
    !
    !    Typically the subroutine is used to integrate from T to TOUT but it
    !    can be used as a one-step integrator to advance the solution a
    !    single step in the direction of TOUT.  On return, the parameters in
    !    the call list are set for continuing the integration.  The user has
    !    only to call again (and perhaps define a new value for TOUT).
    !
    !    Before the first call, the user must
    !
    !    * supply the subroutine F(T,Y,YP) to evaluate the right hand side;
    !      and declare F in an EXTERNAL statement;
    !
    !    * initialize the parameters:
    !      NEQN, Y(1:NEQN), T, TOUT, RELERR, ABSERR, FLAG.
    !      In particular, T should initially be the starting point for integration,
    !      Y should be the value of the initial conditions, and FLAG should
    !      normally be +1.
    !
    !    Normally, the user only sets the value of FLAG before the first call, and
    !    thereafter, the program manages the value.  On the first call, FLAG should
    !    normally be +1 (or -1 for single step mode.)  On normal return, FLAG will
    !    have been reset by the program to the value of 2 (or -2 in single
    !    step mode), and the user can continue to call the routine with that
    !    value of FLAG.
    !
    !    (When the input magnitude of FLAG is 1, this indicates to the program
    !    that it is necessary to do some initialization work.  An input magnitude
    !    of 2 lets the program know that that initialization can be skipped,
    !    and that useful information was computed earlier.)
    !
    !    The routine returns with all the information needed to continue
    !    the integration.  If the integration reached TOUT, the user need only
    !    define a new TOUT and call again.  In the one-step integrator
    !    mode, returning with FLAG = -2, the user must keep in mind that
    !    each step taken is in the direction of the current TOUT.  Upon
    !    reaching TOUT, indicated by the output value of FLAG switching to 2,
    !    the user must define a new TOUT and reset FLAG to -2 to continue
    !    in the one-step integrator mode.
    !
    !    In some cases, an error or difficulty occurs during a call.  In that case,
    !    the output value of FLAG is used to indicate that there is a problem
    !    that the user must address.  These values include:
    !
    !    * 3, integration was not completed because the input value of RELERR, the
    !      relative error tolerance, was too small.  RELERR has been increased
    !      appropriately for continuing.  If the user accepts the output value of
    !      RELERR, then simply reset FLAG to 2 and continue.
    !
    !    * 4, integration was not completed because more than MAXNFE derivative
    !      evaluations were needed.  This is approximately (MAXNFE/6) steps.
    !      The user may continue by simply calling again.  The function counter
    !      will be reset to 0, and another MAXNFE function evaluations are allowed.
    !
    !    * 5, integration was not completed because the solution vanished,
    !      making a pure relative error test impossible.  The user must use
    !      a non-zero ABSERR to continue.  Using the one-step integration mode
    !      for one step is a good way to proceed.
    !
    !    * 6, integration was not completed because the requested accuracy
    !      could not be achieved, even using the smallest allowable stepsize.
    !      The user must increase the error tolerances ABSERR or RELERR before
    !      continuing.  It is also necessary to reset FLAG to 2 (or -2 when
    !      the one-step integration mode is being used).  The occurrence of
    !      FLAG = 6 indicates a trouble spot.  The solution is changing
    !      rapidly, or a singularity may be present.  It often is inadvisable
    !      to continue.
    !
    !    * 7, it is likely that this routine is inefficient for solving
    !      this problem.  Too much output is restricting the natural stepsize
    !      choice.  The user should use the one-step integration mode with
    !      the stepsize determined by the code.  If the user insists upon
    !      continuing the integration, reset FLAG to 2 before calling
    !      again.  Otherwise, execution will be terminated.
    !
    !    * 8, invalid input parameters, indicates one of the following:
    !      NEQN <= 0;
    !      T = TOUT and |FLAG| /= 1;
    !      RELERR < 0 or ABSERR < 0;
    !      FLAG == 0  or FLAG < -2 or 8 < FLAG.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 April 2011
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Herman Watts, Lawrence Shampine.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Erwin Fehlberg,
    !    Low-order Classical Runge-Kutta Formulas with Stepsize Control,
    !    NASA Technical Report R-315, 1969.
    !
    !    Lawrence Shampine, Herman Watts, S Davenport,
    !    Solving Non-stiff Ordinary Differential Equations - The State of the Art,
    !    SIAM Review,
    !    Volume 18, pages 376-411, 1976.
    !
    !  Parameters:
    !
    !    Input, external F, a user-supplied subroutine to evaluate the
    !    derivatives Y'(T), of the form:
    !
    !      subroutine f ( t, y, yp )
    !      real ( kind = 8 ) t
    !      real ( kind = 8 ) y(*)
    !      real ( kind = 8 ) yp(*)
    !
    !    Input, integer ( kind = 4 ) NEQN, the number of equations to be integrated.
    !
    !    Input/output, real ( kind = 8 ) Y(NEQN), the current solution vector at T.
    !
    !    Input/output, real ( kind = 8 ) YP(NEQN), the current value of the
    !    derivative of the dependent variable.  The user should not set or alter
    !    this information!
    !
    !    Input/output, real ( kind = 8 ) T, the current value of the independent
    !    variable.
    !
    !    Input, real ( kind = 8 ) TOUT, the output point at which solution is
    !    desired.  TOUT = T is allowed on the first call only, in which case
    !    the routine returns with FLAG = 2 if continuation is possible.
    !
    !    Input, real ( kind = 8 ) RELERR, ABSERR, the relative and absolute
    !    error tolerances for the local error test.  At each step the code
    !    requires:
    !      abs ( local error ) <= RELERR * abs ( Y ) + ABSERR
    !    for each component of the local error and the solution vector Y.
    !    RELERR cannot be "too small".  If the routine believes RELERR has been
    !    set too small, it will reset RELERR to an acceptable value and return
    !    immediately for user action.
    !
    !    Input/output, integer ( kind = 4 ) FLAG, indicator for status of
    !    integration.  On the first call, set FLAG to +1 for normal use, or to -1
    !    for single step mode.  On return, a value of 2 or -2 indicates normal
    !    progress, while any other value indicates a problem that should
    !    be addressed.
    !
    implicit none
    
    procedure(r8_deriv_vector), pointer, intent(in) :: f
    integer ( kind = 4 ), intent(in) :: neqn
    real ( kind = 8 ), intent(inout) ::  t
    real ( kind = 8 ), intent(inout) ::  y(neqn)
    real ( kind = 8 ), intent(inout) ::  yp(neqn)
    real ( kind = 8 ), intent(in) ::  tout
    real ( kind = 8 ), intent(inout) ::  relerr
    real ( kind = 8 ), intent(in) ::  abserr
    integer ( kind = 4 ), intent(inout) :: flag
    real ( kind = 8 ), intent(out) :: h
    
    integer ( kind = 4 ), parameter :: maxnfe = 3000
    real ( kind = 8 ), parameter :: remin = 1.0E-12
    
    real ( kind = 8 ) :: relerr_save 
    real ( kind = 8 ) :: abserr_save 
    real ( kind = 8 ) ae
    real ( kind = 8 ) dt
    real ( kind = 8 ) ee
    real ( kind = 8 ) eeoet
    real ( kind = 8 ) eps
    real ( kind = 8 ) esttol
    real ( kind = 8 ) et
    real ( kind = 8 ) f1(neqn)
    real ( kind = 8 ) f2(neqn)
    real ( kind = 8 ) f3(neqn)
    real ( kind = 8 ) f4(neqn)
    real ( kind = 8 ) f5(neqn)
    logical hfaild
    real ( kind = 8 ) hmin
    integer ( kind = 4 ) k
    integer ( kind = 4 ) mflag
    integer ( kind = 4 ) :: flag_save 
    integer ( kind = 4 ) :: init 
    integer ( kind = 4 ) :: kflag 
    integer ( kind = 4 ) :: kop 
    integer ( kind = 4 ) :: nfe 
    logical output
    real ( kind = 8 ) relerr_min
    real ( kind = 8 ) s
    real ( kind = 8 ) scale
    real ( kind = 8 ) tol
    real ( kind = 8 ) toln
    real ( kind = 8 ) ypk
    !
    ! Default Values
    !
    h = -1.0D+00
    relerr_save = -1.0D+00
    abserr_save = -1.0D+00
    flag_save = -1000
    init = -1000
    kflag = -1000
    kop = -1
    nfe = -1
    
    !
    !  Check the input parameters.
    !
    eps = epsilon ( eps )

    if ( neqn < 1 ) then
        flag = 8
        return
    end if

    if ( relerr < 0.0D+00 ) then
        flag = 8
        return
    end if

    if ( abserr < 0.0D+00 ) then
        flag = 8
        return
    end if

    if ( flag == 0 .or. 8 < flag .or. flag < -2 ) then
        flag = 8
        return
    end if

    mflag = abs ( flag )
    !
    !  Is this a continuation call?
    !
    if ( mflag /= 1 ) then

        if ( t == tout .and. kflag /= 3 ) then
            flag = 8
            return
        end if

        if ( mflag == 2 ) then

            if ( kflag == 3 ) then

                flag = flag_save
                mflag = abs ( flag )

            else if ( init == 0 ) then

                flag = flag_save

            else if ( kflag == 4 ) then

                nfe = 0

            else if ( kflag == 5 .and. abserr == 0.0D+00 ) then

                error stop ' '//NEW_LINE('a') & 
                    //'R8_RKF45 - Fatal error!'//NEW_LINE('a') & 
                    //'  KFLAG = 5 and ABSERR = 0.0'//NEW_LINE('a')
                
                !write ( *, '(a)' ) ' '
                !write ( *, '(a)' ) 'R8_RKF45 - Fatal error!'
                !write ( *, '(a)' ) '  KFLAG = 5 and ABSERR = 0.0'
                !stop

            else if ( &
                kflag == 6 .and. &
                relerr <= relerr_save .and. &
                abserr <= abserr_save ) then

                error stop ' '//NEW_LINE('a') & 
                    //'R8_RKF45 - Fatal error!'//NEW_LINE('a') & 
                    //'  KFLAG = 6 and'//NEW_LINE('a') & 
                    //'  RELERR <= RELERR_SAVE and'//NEW_LINE('a') & 
                    //'  ABSERR <= ABSERR_SAVE'//NEW_LINE('a') 
                
                !write ( *, '(a)' ) ' '
                !write ( *, '(a)' ) 'R8_RKF45 - Fatal error!'
                !write ( *, '(a)' ) '  KFLAG = 6 and'
                !write ( *, '(a)' ) '  RELERR <= RELERR_SAVE and'
                !write ( *, '(a)' ) '  ABSERR <= ABSERR_SAVE'
                !stop

            end if
            !
            !  FLAG = 3, 4, 5, 6, 7 or 8.
            !
        else

            if ( flag == 3 ) then

                flag = flag_save
                if ( kflag == 3 ) then
                    mflag = abs ( flag )
                end if

            else if ( flag == 4 ) then

                nfe = 0
                flag = flag_save
                if ( kflag == 3 ) then
                    mflag = abs ( flag )
                end if

            else if ( flag == 5 .and. 0.0D+00 < abserr ) then

                flag = flag_save
                if ( kflag == 3 ) then
                    mflag = abs ( flag )
                end if
                !
                !  Integration cannot be continued because the user did not respond to
                !  the instructions pertaining to FLAG = 5, 6, 7 or 8.
                !
            else
                
                error stop ' '//NEW_LINE('a') & 
                    //'R8_RKF45 - Fatal error!'//NEW_LINE('a') & 
                    //'  Integration cannot be continued.'//NEW_LINE('a') & 
                    //'  The user did not respond to the output'//NEW_LINE('a') & 
                    //'  value FLAG = 5, 6, 7, or 8.'//NEW_LINE('a') 
                
                !write ( *, '(a)' ) ' '
                !write ( *, '(a)' ) 'R8_RKF45 - Fatal error!'
                !write ( *, '(a)' ) '  Integration cannot be continued.'
                !write ( *, '(a)' ) '  The user did not respond to the output'
                !write ( *, '(a)' ) '  value FLAG = 5, 6, 7, or 8.'
                !stop
            end if

        end if

    end if
    !
    !  Save the input value of FLAG.
    !  Set the continuation flag KFLAG for subsequent input checking.
    !
    flag_save = flag
    kflag = 0
    !
    !  Save RELERR and ABSERR for checking input on subsequent calls.
    !
    relerr_save = relerr
    abserr_save = abserr
    !
    !  Restrict the relative error tolerance to be at least
    !
    !    2 * EPS + REMIN
    !
    !  to avoid limiting precision difficulties arising from impossible
    !  accuracy requests.
    !
    relerr_min = 2.0D+00 * epsilon ( relerr_min ) + remin
    !
    !  Is the relative error tolerance too small?
    !
    if ( relerr < relerr_min ) then
        relerr = relerr_min
        flag = 3
        kflag = 3
        return
    end if

    dt = tout - t
    !
    !  Initialization:
    !
    !  Set the initialization completion indicator, INIT;
    !  set the indicator for too many output points, KOP;
    !  evaluate the initial derivatives;
    !  set the counter for function evaluations, NFE;
    !  estimate the starting stepsize.
    !
    if ( mflag == 1 ) then

        init = 0
        kop = 0
        call f ( t, neqn, y, yp )
        nfe = 1

        if ( t == tout ) then
            flag = 2
            return
        end if

    end if

    if ( init == 0 ) then

        init = 1
        h = abs ( dt )
        toln = 0.0D+00

        do k = 1, neqn
            tol = relerr * abs ( y(k) ) + abserr
            if ( 0.0D+00 < tol ) then
                toln = tol
                ypk = abs ( yp(k) )
                if ( tol < ypk * h**5 ) then
                    h = ( tol / ypk )**0.2D+00
                end if
            end if
        end do

        if ( toln <= 0.0D+00 ) then
            h = 0.0D+00
        end if

        h = max ( h, 26.0D+00 * eps * max ( abs ( t ), abs ( dt ) ) )
        flag_save = sign ( 2, flag )

    end if
    !
    !  Set the stepsize for integration in the direction from T to TOUT.
    !
    h = sign ( h, dt )
    !
    !  Test to see if too may output points are being requested.
    !
    if ( 2.0D+00 * abs ( dt ) <= abs ( h ) ) then
        kop = kop + 1
    end if
    !
    !  Unnecessary frequency of output.
    !  Seems like an error I'm willing to tolerate!
    !
    ! if ( kop == 100 ) then
    if ( kop == 10000 ) then
        kop = 0
        flag = 7
        return
    end if
    !
    !  If we are too close to the output point, then simply extrapolate and return.
    !
    if ( abs ( dt ) <= 26.0D+00 * eps * abs ( t ) ) then
        t = tout
        y(1:neqn) = y(1:neqn) + dt * yp(1:neqn)
        call f ( t, neqn, y, yp )
        nfe = nfe + 1
        flag = 2
        return
    end if
    !
    !  Initialize the output point indicator.
    !
    output = .false.
    !
    !  To avoid premature underflow in the error tolerance function,
    !  scale the error tolerances.
    !
    scale = 2.0D+00 / relerr
    ae = scale * abserr
    !
    !  Step by step integration.
    !
    do

        hfaild = .false.
        !
        !  Set the smallest allowable stepsize.
        !
        hmin = 26.0D+00 * eps * abs ( t )
        !
        !  Adjust the stepsize if necessary to hit the output point.
        !
        !  Look ahead two steps to avoid drastic changes in the stepsize and
        !  thus lessen the impact of output points on the code.
        !
        dt = tout - t

        if ( 2.0D+00 * abs ( h ) <= abs ( dt ) ) then

        else
            !
            !  Will the next successful step complete the integration to the output point?
            !
            if ( abs ( dt ) <= abs ( h ) ) then
                output = .true.
                h = dt
            else
                h = 0.5D+00 * dt
            end if

        end if
        !
        !  Here begins the core integrator for taking a single step.
        !
        !  The tolerances have been scaled to avoid premature underflow in
        !  computing the error tolerance function ET.
        !  To avoid problems with zero crossings, relative error is measured
        !  using the average of the magnitudes of the solution at the
        !  beginning and end of a step.
        !  The error estimate formula has been grouped to control loss of
        !  significance.
        !
        !  To distinguish the various arguments, H is not permitted
        !  to become smaller than 26 units of roundoff in T.
        !  Practical limits on the change in the stepsize are enforced to
        !  smooth the stepsize selection process and to avoid excessive
        !  chattering on problems having discontinuities.
        !  To prevent unnecessary failures, the code uses 9/10 the stepsize
        !  it estimates will succeed.
        !
        !  After a step failure, the stepsize is not allowed to increase for
        !  the next attempted step.  This makes the code more efficient on
        !  problems having discontinuities and more effective in general
        !  since local extrapolation is being used and extra caution seems
        !  warranted.
        !
        !  Test the number of derivative function evaluations.
        !  If okay, try to advance the integration from T to T+H.
        !
        do
            !
            !  Have we done too much work?
            !
            if ( maxnfe < nfe ) then
                flag = 4
                kflag = 4
                return
            end if
            !
            !  Advance an approximate solution over one step of length H.
            !
            call r8_fehl ( f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, f1 )
            nfe = nfe + 5
            !
            !  Compute and test allowable tolerances versus local error estimates
            !  and remove scaling of tolerances.  The relative error is
            !  measured with respect to the average of the magnitudes of the
            !  solution at the beginning and end of the step.
            !
            eeoet = 0.0D+00

            do k = 1, neqn

                et = abs ( y(k) ) + abs ( f1(k) ) + ae

                if ( et <= 0.0D+00 ) then
                    flag = 5
                    return
                end if

                ee = abs &
                    ( ( -2090.0D+00 * yp(k) &
                    + ( 21970.0D+00 * f3(k) - 15048.0D+00 * f4(k) ) &
                    ) &
                    + ( 22528.0D+00 * f2(k) - 27360.0D+00 * f5(k) ) &
                    )

                eeoet = max ( eeoet, ee / et )

            end do

            esttol = abs ( h ) * eeoet * scale / 752400.0D+00

            if ( esttol <= 1.0D+00 ) then
                exit
            end if
            !
            !  Unsuccessful step.  Reduce the stepsize, try again.
            !  The decrease is limited to a factor of 1/10.
            !
            hfaild = .true.
            output = .false.

            if ( esttol < 59049.0D+00 ) then
                s = 0.9D+00 / esttol**0.2D+00
            else
                s = 0.1D+00
            end if

            h = s * h

            if ( abs ( h ) < hmin ) then
                flag = 6
                kflag = 6
                return
            end if

        end do
        !
        !  We exited the loop because we took a successful step.
        !  Store the solution for T+H, and evaluate the derivative there.
        !
        t = t + h
        y(1:neqn) = f1(1:neqn)
        call f ( t, neqn, y, yp )
        nfe = nfe + 1
        !
        !  Choose the next stepsize.  The increase is limited to a factor of 5.
        !  If the step failed, the next stepsize is not allowed to increase.
        !
        if ( 0.0001889568D+00 < esttol ) then
            s = 0.9D+00 / esttol**0.2D+00
        else
            s = 5.0D+00
        end if

        if ( hfaild ) then
            s = min ( s, 1.0D+00 )
        end if

        h = sign ( max ( s * abs ( h ), hmin ), h )
        !
        !  End of core integrator
        !
        !  Should we take another step?
        !
        if ( output ) then
            t = tout
            flag = 2
            return
        end if

        if ( flag <= 0 ) then
            exit
        end if

    end do
    !
    !  One step integration mode.
    !
    flag = -2

    return
    end

    subroutine srk4_test()
    !call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RKF45_TEST'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) '  Test the RKF45 library.'

    call test01 ( )
    call test02 ( )
    call test03 ( )
    call test04 ( )
    call test05 ( )
    call test06 ( )
    !
    !  Terminate.
    !
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RKF45_TEST'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
    !call timestamp ( )

    contains

    subroutine test01 ( )

    !*****************************************************************************80
    !
    !! TEST01 solves a scalar ODE.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 June 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    implicit none

    integer ( kind = 4 ), parameter :: neqn = 1

    real ( kind = 4 ) abserr
    integer ( kind = 4 ) flag
    integer ( kind = 4 ) i_step
    integer ( kind = 4 ) n_step
    real ( kind = 4 ) relerr
    real ( kind = 4 ) t
    real ( kind = 4 ) t_out
    real ( kind = 4 ) t_start
    real ( kind = 4 ) t_stop
    real ( kind = 4 ) y(neqn)
    real ( kind = 4 ) yp(neqn)
    real ( kind = 4 ) h

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST01'
    write ( *, '(a)' ) '  Solve a scalar equation using R4_RKF45:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Y'' = 0.25 * Y * ( 1 - Y / 20 )'

    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    flag = 1

    t_start = 0.0E+00
    t_stop = 20.0E+00

    n_step = 5

    t_out = 0.0E+00
    t = t_out
    y(1) = 1.0E+00
    call r4_f1 ( t, neqn, y, yp )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  FLAG     T             Y            Y''           Y_Exact         Error'
    write ( *, '(a)' ) ' '
    write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r4_y1x ( t ), &
        y(1) - r4_y1x ( t )

    do i_step = 1, n_step

        t = ( real ( n_step - i_step + 1, kind = 4 ) * t_start  &
            + real (          i_step - 1, kind = 4 ) * t_stop ) &
            / real ( n_step,              kind = 4 )

        t_out = ( real ( n_step - i_step, kind = 4 ) * t_start  &
            + real (          i_step, kind = 4 ) * t_stop ) &
            / real ( n_step,          kind = 4 )

        call r4_rkf45 ( r4_f1, neqn, y, yp, t, t_out, relerr, abserr, flag, h)

        write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r4_y1x ( t ), &
            y(1) - r4_y1x ( t )

    end do

    return
    end
    subroutine test02 ( )

    !*****************************************************************************80
    !
    !! TEST02 solves a vector ODE.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 June 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    implicit none

    integer ( kind = 4 ), parameter :: neqn = 2

    real ( kind = 4 ) abserr
    integer ( kind = 4 ) flag
    integer ( kind = 4 ) i_step
    integer ( kind = 4 ) n_step
    real ( kind = 4 ) relerr
    real ( kind = 4 ) t
    real ( kind = 4 ) t_out
    real ( kind = 4 ) t_start
    real ( kind = 4 ) t_stop
    real ( kind = 4 ) y(neqn)
    real ( kind = 4 ) yp(neqn)
    real ( kind = 4 ) h

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST02'
    write ( *, '(a)' ) '  Solve a vector equation using R4_RKF45:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Y''(1) =  Y(2)'
    write ( *, '(a)' ) '  Y''(2) = -Y(1)'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This system is equivalent to the following'
    write ( *, '(a)' ) '  second order system:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Z" = - Z.'

    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    flag = 1

    t_start = 0.0E+00
    t_stop = 2.0E+00 * 3.14159265E+00

    n_step = 12

    t = 0.0E+00
    t_out = 0.0E+00

    y(1) = 1.0E+00
    y(2) = 0.0E+00
    call r4_f2 ( t, 2, y, yp )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  FLAG       T          Y(1)          Y(2)'
    write ( *, '(a)' ) ' '
    write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

    do i_step = 1, n_step

        t = ( real ( n_step - i_step + 1, kind = 4 ) * t_start &
            + real (          i_step - 1, kind = 4 ) * t_stop ) &
            / real ( n_step,              kind = 4 )

        t_out = ( real ( n_step - i_step, kind = 4 ) * t_start &
            + real (          i_step, kind = 4 ) * t_stop ) &
            / real ( n_step,          kind = 4 )

        call r4_rkf45 ( r4_f2, neqn, y, yp, t, t_out, relerr, abserr, flag, h )

        write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

    end do

    return
    end
    subroutine test03 ( )

    !*****************************************************************************80
    !
    !! TEST03 solves a scalar ODE and uses one-step integration.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 June 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    implicit none

    integer ( kind = 4 ), parameter :: neqn = 1

    real ( kind = 4 ) abserr
    integer ( kind = 4 ) flag
    integer ( kind = 4 ) i_step
    integer ( kind = 4 ) n_step
    real ( kind = 4 ) relerr
    real ( kind = 4 ) t
    real ( kind = 4 ) t_out
    real ( kind = 4 ) t_start
    real ( kind = 4 ) t_stop
    real ( kind = 4 ) y(neqn)
    real ( kind = 4 ) yp(neqn)
    real ( kind = 4 ) h

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST03'
    write ( *, '(a)' ) '  Solve a scalar equation using R4_RKF45:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Y'' = 0.25 * Y * ( 1 - Y / 20 )'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Use the special SINGLE_STEP mode'
    write ( *, '(a)' ) '  which returns after every step.'

    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    flag = -1

    t_start = 0.0E+00
    t_stop = 20.0E+00

    n_step = 5

    t = 0.0E+00
    t_out = 0.0E+00
    y(1) = 1.0E+00
    call r4_f1 ( t, neqn, y, yp )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  FLAG     T             Y           Y''         Y_Exact        Error'
    write ( *, '(a)' ) ' '
    write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r4_y1x ( t ), &
        y(1) - r4_y1x ( t )

    do i_step = 1, n_step

        t = ( real ( n_step - i_step + 1, kind = 4 ) * t_start  &
            + real (          i_step - 1, kind = 4 ) * t_stop ) &
            / real ( n_step,              kind = 4 )

        t_out = ( real ( n_step - i_step, kind = 4 ) * t_start  &
            + real (          i_step, kind = 4 ) * t_stop ) &
            / real ( n_step,          kind = 4 )
        !
        !  As long as FLAG is negative, we are heading towards T_OUT, but
        !  have not reached it!
        !
        do while ( flag < 0 )

            call r4_rkf45 ( r4_f1, neqn, y, yp, t, t_out, relerr, abserr, flag, h )

            write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r4_y1x ( t ), &
                y(1) - r4_y1x ( t )

        end do
        !
        !  FLAG is returned as +2 when we reach T_OUT.  Reset it to -2
        !  to continue to the next T_OUT in one step mode.
        !
        flag = -2

    end do

    return
    end
    subroutine test04 ( )

    !*****************************************************************************80
    !
    !! TEST04 solves a scalar ODE.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 June 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    implicit none

    integer ( kind = 4 ), parameter :: neqn = 1

    integer ( kind = 4 ) flag
    integer ( kind = 4 ) i_step
    integer ( kind = 4 ) n_step
    real ( kind = 8 ) abserr
    real ( kind = 8 ) relerr
    real ( kind = 8 ) t
    real ( kind = 8 ) t_out
    real ( kind = 8 ) t_start
    real ( kind = 8 ) t_stop
    real ( kind = 8 ) y(neqn)
    real ( kind = 8 ) yp(neqn)
    real ( kind = 8 ) h

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST04'
    write ( *, '(a)' ) '  Solve a scalar equation using R8_RKF45:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Y'' = 0.25 * Y * ( 1 - Y / 20 )'

    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    flag = 1

    t_start = 0.0D+00
    t_stop = 20.0D+00

    n_step = 5

    t_out = 0.0D+00
    t = t_out
    y(1) = 1.0D+00
    call r8_f1 ( t, neqn, y, yp )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  FLAG     T             Y            Y''           Y_Exact         Error'
    write ( *, '(a)' ) ' '
    write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r8_y1x ( t ), &
        y(1) - r8_y1x ( t )

    do i_step = 1, n_step

        t = ( real ( n_step - i_step + 1, kind = 8 ) * t_start  &
            + real (          i_step - 1, kind = 8 ) * t_stop ) &
            / real ( n_step,              kind = 8 )

        t_out = ( real ( n_step - i_step, kind = 8 ) * t_start  &
            + real (          i_step, kind = 8 ) * t_stop ) &
            / real ( n_step,          kind = 8 )

        call r8_rkf45 ( r8_f1, neqn, y, yp, t, t_out, relerr, abserr, flag, h )

        write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r8_y1x ( t ), &
            y(1) - r8_y1x ( t )

    end do

    return
    end
    subroutine test05 ( )

    !*****************************************************************************80
    !
    !! TEST05 solves a vector ODE.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 June 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    implicit none

    integer ( kind = 4 ), parameter :: neqn = 2

    real ( kind = 8 ) abserr
    integer ( kind = 4 ) flag
    integer ( kind = 4 ) i_step
    integer ( kind = 4 ) n_step
    real ( kind = 8 ) relerr
    real ( kind = 8 ) t
    real ( kind = 8 ) t_out
    real ( kind = 8 ) t_start
    real ( kind = 8 ) t_stop
    real ( kind = 8 ) y(neqn)
    real ( kind = 8 ) yp(neqn)
    real ( kind = 8 ) h

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST05'
    write ( *, '(a)' ) '  Solve a vector equation using R8_RKF45:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Y''(1) =  Y(2)'
    write ( *, '(a)' ) '  Y''(2) = -Y(1)'

    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    flag = 1

    t_start = 0.0D+00
    t_stop = 2.0D+00 * 3.14159265D+00

    n_step = 12

    t = 0.0D+00
    t_out = 0.0D+00

    y(1) = 1.0D+00
    y(2) = 0.0D+00
    call r8_f2 ( t, neqn, y, yp )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  FLAG       T          Y(1)          Y(2)'
    write ( *, '(a)' ) ' '
    write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

    do i_step = 1, n_step

        t = ( real ( n_step - i_step + 1, kind = 8 ) * t_start &
            + real (          i_step - 1, kind = 8 ) * t_stop ) &
            / real ( n_step,              kind = 8 )

        t_out = ( real ( n_step - i_step, kind = 8 ) * t_start &
            + real (          i_step, kind = 8 ) * t_stop ) &
            / real ( n_step,          kind = 8 )

        call r8_rkf45 ( r8_f2, neqn, y, yp, t, t_out, relerr, abserr, flag, h )

        write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

    end do

    return
    end
    subroutine test06 ( )

    !*****************************************************************************80
    !
    !! TEST06 solves a scalar ODE and uses one-step integration.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 June 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    implicit none

    integer ( kind = 4 ), parameter :: neqn = 1

    real ( kind = 8 ) abserr
    integer ( kind = 4 ) flag
    integer ( kind = 4 ) i_step
    integer ( kind = 4 ) n_step
    real ( kind = 8 ) relerr
    real ( kind = 8 ) t
    real ( kind = 8 ) t_out
    real ( kind = 8 ) t_start
    real ( kind = 8 ) t_stop
    real ( kind = 8 ) y(neqn)
    real ( kind = 8 ) yp(neqn)
    real ( kind = 8 ) h

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST06'
    write ( *, '(a)' ) '  Solve a scalar equation using R8_RKF45:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Y'' = 0.25 * Y * ( 1 - Y / 20 )'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Use the special SINGLE_STEP mode'
    write ( *, '(a)' ) '  which returns after every step.'

    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )

    flag = -1

    t_start = 0.0D+00
    t_stop = 20.0D+00

    n_step = 5

    t = 0.0D+00
    t_out = 0.0D+00
    y(1) = 1.0D+00
    call r8_f1 ( t, neqn, y, yp )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  FLAG     T             Y           Y''       Y_Exact        Error'
    write ( *, '(a)' ) ' '
    write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r8_y1x ( t ), &
        y(1) - r8_y1x ( t )

    do i_step = 1, n_step

        t = ( real ( n_step - i_step + 1, kind = 8 ) * t_start  &
            + real (          i_step - 1, kind = 8 ) * t_stop ) &
            / real ( n_step,              kind = 8 )

        t_out = ( real ( n_step - i_step, kind = 8 ) * t_start  &
            + real (          i_step, kind = 8 ) * t_stop ) &
            / real ( n_step,          kind = 8 )
        !
        !  As long as FLAG is negative, we are heading towards T_OUT, but
        !  have not reached it!
        !
        do while ( flag < 0 )

            call r8_rkf45 ( r8_f1, neqn, y, yp, t, t_out, relerr, abserr, flag, h )

            write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), yp(1), r8_y1x ( t ), &
                y(1) - r8_y1x ( t )

        end do
        !
        !  FLAG is returned as +2 when we reach T_OUT.  Reset it to -2
        !  to continue to the next T_OUT in one step mode.
        !
        flag = -2

    end do

    return
    end
    pure subroutine r4_f1 ( t, n, y, yp )

    !*****************************************************************************80
    !
    !! R4_F1 evaluates the derivative for the ODE.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 March 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 4 ) T, the value of the independent variable.
    !
    !    Input, real ( kind = 4 ) Y(NEQN), the value of the dependent variable.
    !
    !    Output, real ( kind = 4 ) YP(NEQN), the value of the derivative
    !    dY(1:NEQN)/dT.
    !
    implicit none

    real ( kind = 4 ), intent(in) :: t
    integer (kind = 4), intent(in) :: n
    real ( kind = 4 ), intent(in) :: y(n)
    real ( kind = 4 ), intent(out) :: yp(n)

    yp(1) = 0.25E+00 * y(1) * ( 1.0E+00 - y(1) / 20.0E+00 )

    return
    end
    
    pure function r4_y1x ( t )

    !*****************************************************************************80
    !
    !! R4_Y1X evaluates the exact solution of the ODE.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 March 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 4 ) T, the value of the independent variable.
    !
    !    Output, real ( kind = 4 ) Y1X_S, the exact solution.
    !
    implicit none

    real ( kind = 4 ), intent(in) :: t
    real ( kind = 4 ) :: r4_y1x

    r4_y1x = 20.0E+00 / ( 1.0E+00 + 19.0E+00 * exp ( - 0.25E+00 * t ) )

    return
    end
    pure subroutine r4_f2 ( t, n, y, yp )

    !*****************************************************************************80
    !
    !! R4_F2 evaluates the derivative for the ODE.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 March 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 4 ) T, the value of the independent variable.
    !
    !    Input, real ( kind = 4 ) Y(NEQN), the value of the dependent variable.
    !
    !    Output, real ( kind = 4 ) YP(NEQN), the value of the derivative
    !    dY(1:NEQN)/dT.
    !
    implicit none

    real ( kind = 4 ), intent(in) :: t
    integer ( kind = 4), intent(in) :: n
    real ( kind = 4 ), intent(in) :: y(n)
    real ( kind = 4 ), intent(out) :: yp(n)

    yp(1) =  y(2)
    yp(2) = -y(1)

    return
    end
    pure subroutine r8_f1 ( t, n, y, yp )

    !*****************************************************************************80
    !
    !! R8_F1 evaluates the derivative for the ODE.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 March 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T, the value of the independent variable.
    !
    !    Input, real ( kind = 8 ) Y(NEQN), the value of the dependent variable.
    !
    !    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative
    !    dY(1:NEQN)/dT.
    !
    implicit none

    real ( kind = 8 ), intent(in) :: t
    integer ( kind = 4), intent(in) :: n
    real ( kind = 8 ), intent(in) :: y(n)
    real ( kind = 8 ), intent(out) :: yp(n)

    yp(1) = 0.25D+00 * y(1) * ( 1.0D+00 - y(1) / 20.0D+00 )

    return
    end
    
    pure function r8_y1x ( t )

    !*****************************************************************************80
    !
    !! R8_Y1X evaluates the exact solution of the ODE.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 March 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T, the value of the independent variable.
    !
    !    Output, real ( kind = 8 ) Y1X_D, the exact solution.
    !
    implicit none

    real ( kind = 8 ), intent(in) :: t
    real ( kind = 8 ) :: r8_y1x

    r8_y1x = 20.0D+00 / ( 1.0D+00 + 19.0D+00 * exp ( - 0.25D+00 * t ) )

    return
    end
    pure subroutine r8_f2 ( t, n, y, yp )

    !*****************************************************************************80
    !
    !! R8_F2 evaluates the derivative for the ODE.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 March 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T, the value of the independent variable.
    !
    !    Input, real ( kind = 8 ) Y(NEQN), the value of the dependent variable.
    !
    !    Output, real ( kind = 8 ) YP(NEQN), the value of the derivative
    !    dY(1:NEQN)/dT.
    !
    implicit none

    real ( kind = 8 ), intent(in) :: t
    integer( kind = 4), intent(in) :: n
    real ( kind = 8 ), intent(in) :: y(n)
    real ( kind = 8 ), intent(out) :: yp(n)

    yp(1) =  y(2)
    yp(2) = -y(1)

    return
    end


    end subroutine


    end module