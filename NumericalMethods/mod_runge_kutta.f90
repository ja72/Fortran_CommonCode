!***************************************************************************
! PURPOSE:      : MODULE DEFININING ODE SCALAR INTEGRATORS.
! DEPENDENCIES  : MOD_COMMON
! HISTORY       : Oct 2014  Inittial Write                  [John Alexiou]
!               : Nov 2023  Add array suppport              [John Alexiou]

    MODULE mod_runge_kutta
    use mod_common
    IMPLICIT NONE
    
        integer, parameter :: wp = real64
        
        interface 
            pure function f_scalar(t, x) result(xp)
            import
            real(wp), intent(in) :: t,x
            real(wp) :: xp
            end function f_scalar
            
            pure function f_vec(t, x) result(xp)
            import
            real(wp), intent(in) :: t
            real(wp), intent(in) :: x(:)
            real(wp), allocatable :: xp(:)
            end function f_vec
        end interface 
        

    CONTAINS

    !-----------------------------------------------------------------------------
    !   RK45 - Runge-Kutta-Fehlberg method for solving an initial value problem
    !
    pure subroutine rk45_scalar(f,t,x,h,e)    
        procedure(f_scalar), pointer :: f
        real(wp), intent(inout):: t,x,h,e      
        real(wp) :: f1,f2,f3,f4,f5,f6,x5
        
        real(wp), parameter :: c20=0.25, c21=0.25, &        
        c30=0.375, c31=0.09375, c32=0.28125, &
        c40=0.92307692307692, &
        c41=0.87938097405553, c42=-3.2771961766045, c43=3.3208921256258, &
        c51=2.0324074074074, c52=-8.0, c53=7.1734892787524,&
        c54=-0.20589668615984,       & 
        c60=0.5, c61=-0.2962962962963, c62=2.0, &
        c63=-1.3816764132554, c64=0.45297270955166, c65=-0.275, &
        a1=0.11574074074074, a2=0, a3=0.54892787524366, &
        a4=0.5353313840156, a5=-0.2,        &
        b1=0.11851851851852, b2=0.0, b3=0.51898635477583, &
        b4=0.50613149034201, b5=-0.18,      &
        b6=0.036363636363636                    
   

        f1 = h*f(t,x) 
        f2 = h*f(t+ c20*h,x + c21*f1)  
        f3 = h*f(t+ c30*h,x + c31*f1 + c32*f2)   
        f4 = h*f(t+ c40*h,x + c41*f1 + c42*f2 + c43*f3)  
        f5 = h*f(t+h,x + c51*f1 + c52*f2 + c53*f3 + c54*f4)       
        f6 = h*f(t+ c60*h,x + c61*f1 + c62*f2 + c63*f3 + c64*f4 + c65*f5)
        x5 = x + b1*f1 + b3*f3 + b4*f4 + b5*f5 + b6*f6      
        x  = x + a1*f1 + a3*f3 + a4*f4 + a5*f5    
        t = t + h   
        e = abs(x - x5)     
    end subroutine rk45_scalar
    
    !----------------------------------------------------------------------------------
    ! RK45AD - Adaptive scheme based on Runge-Kutta-Fehlberg method
    ! url: http://www.ma.utexas.edu/CNA/cheney-kincaid/f90code/CHP10/rk45ad.f90
    !
    !    t= 1.0 
    !    x = 2.0 
    !    h = 7.8125e-3 
    !    tb = 1.5625 
    !    itmax = 100 
    !    emin = 1.0e-8 
    !    emax = 1.0e-4 
    !    hmin = 1.0e-6 
    !    hmax = 1.0
    pure subroutine rk45ad_scalar(f,t,x,h,tb,itmax,emin,emax,hmin,hmax,iflag) 
        procedure(f_scalar), pointer :: f
        real(wp), intent(in):: tb,emin,emax,hmin,hmax     
        integer, intent(out):: iflag
        integer, intent(in):: itmax
        real(wp), intent(inout) :: h, t, x
        real(wp) ::delta, d,xsave,tsave,e
        integer :: k

        delta=0.5e-5
        iflag = 1 
        k = 0
        do 
            k = k + 1
            if (k > itmax) exit
            if(abs(h) < hmin) h = sign(1.0,h)*hmin
            if(abs(h) > hmax) h = sign(1.0,h)*hmax
            d = abs(tb - t)      

            if(d <= abs(h))  then
                iflag = 0 
                if(d <= delta*max(abs(tb),abs(t))) exit    
                h = sign(1.0,h)*d     
            end if      
            xsave = x   
            tsave = t   
            call rk45_scalar(f,t,x,h,e)
            if (iflag == 0) exit
            if(e < emin)  then 
                h = 2.0*h 
            end if
            if (e > emax) then
                h = 0.5*h 
                x = xsave       
                t = tsave
                k = k - 1
            end if      
        end do 
    end subroutine rk45ad_scalar
    
    !-----------------------------------------------------------------------------
    !   RK45_ARRAY - Runge-Kutta-Fehlberg method for solving an initial value problem
    !
    pure subroutine rk45_vec(f,t,x,h,e)    
        procedure(f_vec), pointer :: f
        real(wp), intent(inout):: t,h,e      
        real(wp), intent(inout) :: x(:)
        real(wp), dimension(:), allocatable :: f1,f2,f3,f4,f5,f6,x5
        
        real(wp), parameter :: c20=0.25, c21=0.25, &        
        c30=0.375, c31=0.09375, c32=0.28125, &
        c40=0.92307692307692, &
        c41=0.87938097405553, c42=-3.2771961766045, c43=3.3208921256258, &
        c51=2.0324074074074, c52=-8.0, c53=7.1734892787524,&
        c54=-0.20589668615984,       & 
        c60=0.5, c61=-0.2962962962963, c62=2.0, &
        c63=-1.3816764132554, c64=0.45297270955166, c65=-0.275, &
        a1=0.11574074074074, a2=0, a3=0.54892787524366, &
        a4=0.5353313840156, a5=-0.2,        &
        b1=0.11851851851852, b2=0.0, b3=0.51898635477583, &
        b4=0.50613149034201, b5=-0.18,      &
        b6=0.036363636363636                    
   

        f1 = h*f(t,x) 
        f2 = h*f(t+ c20*h,x + c21*f1)  
        f3 = h*f(t+ c30*h,x + c31*f1 + c32*f2)   
        f4 = h*f(t+ c40*h,x + c41*f1 + c42*f2 + c43*f3)  
        f5 = h*f(t+h,x + c51*f1 + c52*f2 + c53*f3 + c54*f4)       
        f6 = h*f(t+ c60*h,x + c61*f1 + c62*f2 + c63*f3 + c64*f4 + c65*f5)
        x5 = x + b1*f1 + b3*f3 + b4*f4 + b5*f5 + b6*f6      
        x  = x + a1*f1 + a3*f3 + a4*f4 + a5*f5    
        t = t + h   
        e = maxval( abs(x - x5) )
    end subroutine rk45_vec

    !----------------------------------------------------------------------------------
    ! RK45AD_ARRAY - Adaptive scheme based on Runge-Kutta-Fehlberg method
    ! url: http://www.ma.utexas.edu/CNA/cheney-kincaid/f90code/CHP10/rk45ad.f90
    !
    !    t= 1.0 
    !    x = 2.0 
    !    h = 7.8125e-3 
    !    tb = 1.5625 
    !    itmax = 100 
    !    emin = 1.0e-8 
    !    emax = 1.0e-4 
    !    hmin = 1.0e-6 
    !    hmax = 1.0
    pure subroutine rk45ad_vec(f,t,x,h,tb,itmax,emin,emax,hmin,hmax,iflag) 
        procedure(f_vec), pointer :: f
        real(wp), intent(in):: tb,emin,emax,hmin,hmax     
        integer, intent(out):: iflag
        integer, intent(in):: itmax
        real(wp), intent(inout) :: h, t
        real(wp), intent(inout) :: x(:)
        real(wp) ::delta, d,tsave,e
        real(wp), dimension(:), allocatable ::xsave
        integer :: k

        delta=0.5e-5
        iflag = 1 
        k = 0
        do 
            k = k + 1
            if (k > itmax) exit
            if(abs(h) < hmin) h = sign(1.0,h)*hmin
            if(abs(h) > hmax) h = sign(1.0,h)*hmax
            d = abs(tb - t)      

            if(d <= abs(h))  then
                iflag = 0 
                if(d <= delta*max(abs(tb),abs(t))) exit    
                h = sign(1.0,h)*d     
            end if      
            xsave = x   
            tsave = t   
            call rk45_vec(f,t,x,h,e)
            if (iflag == 0) exit
            if(e < emin)  then 
                h = 2.0*h 
            end if
            if (e > emax) then
                h = 0.5*h 
                x = xsave       
                t = tsave
                k = k - 1
            end if      
        end do 
    end subroutine rk45ad_vec

    
    END MODULE