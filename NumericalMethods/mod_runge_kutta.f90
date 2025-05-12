!***************************************************************************
! PURPOSE:      : MODULE DEFININING ODE SCALAR INTEGRATORS.
! DEPENDENCIES  : MOD_COMMON
! HISTORY       : Oct 2014  Inittial Write                  [John Alexiou]
!               : Nov 2023  Add array suppport              [John Alexiou]

    MODULE mod_runge_kutta
    use mod_common
    IMPLICIT NONE
    
        interface 
            pure function f_scalar(t, x) result(xp)
            import
            real(real64), intent(in) :: t,x
            real(real64) :: xp
            end function f_scalar
            
            pure function f_vec(t, x) result(xp)
            import
            real(real64), intent(in) :: t
            real(real64), intent(in) :: x(:)
            real(real64), allocatable :: xp(:)
            end function f_vec
        end interface 
        

    CONTAINS

    !-----------------------------------------------------------------------------
    !   RK45 - Runge-Kutta-Fehlberg method for solving an initial value problem
    !
    pure subroutine rk45_scalar(f,t,x,h,e)    
        procedure(f_scalar), pointer :: f
        real(real64), intent(inout):: t,x,h,e      
        real(real64) :: f1,f2,f3,f4,f5,f6, xb, xa
        real(real64) :: t1,t2,t3,t4,t5,t6
        real(real64) :: x1,x2,x3,x4,x5,x6
        
        real(real64), parameter :: c20=0.25, c21=0.25, &        
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
   
        ! f1 = h*f(t ,x)
        ! f2 = h*f(t + c20*h ,x + c21*f1)
        ! f3 = h*f(t + c30*h ,x + c31*f1 + c32*f2 )
        ! f4 = h*f(t + c40*h ,x + c41*f1 + c42*f2 + c43*f3 )
        ! f5 = h*f(t + h, x + c51*f1 + c52*f2 + c53*f3 + c54*f4 )
        ! f6 = h*f(t + c60*h, x + c61*f1 + c62*f2 + c63*f3 + c64*f4 + c65*f5 )
        ! xb = x + b1*f1 + b3*f3 + b4*f4 + b5*f5 + b6*f6
        ! xa  = x + a1*f1 + a3*f3 + a4*f4 + a5*f5    
        
        ! Prepare time steps
        t1 = t
        t2 = ieee_fma(c20, h, t)
        t3 = ieee_fma(c30, h, t)
        t4 = ieee_fma(c40, h, t)
        t5 = t + h
        t6 = ieee_fma(c60, h, t)

        x1 = x
        
        ! Stage 1
        f1 = h*f(t1 ,x1) 
        x2 = ieee_fma(c21, f1, x1)                
        
        ! Stage 2
        f2 = h*f(t2 ,x2)  
        x3 = ieee_fma(c31, f1, x1)
        x3 = ieee_fma(c32, f2, x3)                    
        
        ! Stage 3
        f3 = h*f(t3 ,x3)   
        x4 = ieee_fma(c42, f2, x4)        
        x4 = ieee_fma(c41, f1, x1)
        x4 = ieee_fma(c43, f3, x4)                        
        
        ! Stage 4
        f4 = h*f(t4, x4)
        x5 = ieee_fma(c51, f1, x1)
        x5 = ieee_fma(c52, f2, x5)
        x5 = ieee_fma(c53, f3, x5)
        x5 = ieee_fma(c54, f4, x5)
        
        ! Stage 5
        f5 = h*f(t5, x5)
        x6 = ieee_fma(c61, f1, x1)
        x6 = ieee_fma(c62, f2, x6)
        x6 = ieee_fma(c63, f3, x6)
        x6 = ieee_fma(c64, f4, x6)
        x6 = ieee_fma(c65, f5, x6)
        
        ! Stage 6
        f6 = h*f(t6, x6)
                
        xb = ieee_fma( b1, f1, x1)
        xb = ieee_fma( b3, f3, xb)
        xb = ieee_fma( b4, f4, xb)
        xb = ieee_fma( b5, f5, xb)
        xb = ieee_fma( b6, f6, xb)
        
        xa = ieee_fma( a1, f1, x1)
        xa = ieee_fma( a3, f3, xa)
        xa = ieee_fma( a4, f4, xa)
        xa = ieee_fma( a5, f5, xa)
        
        t = t + h 
        x = xa
        e = abs(xa - xb)     
    end subroutine rk45_scalar
    
    
    !-----------------------------------------------------------------------------
    !   RK45_VEC - Runge-Kutta-Fehlberg method for solving an initial value problem
    !
    pure subroutine rk45_vec(f,t,x,h,e)    
        procedure(f_vec), pointer :: f
        real(real64), intent(inout):: t,h,e
        real(real64), intent(inout) :: x(:)
        real(real64), dimension(size(x)) :: f1,f2,f3,f4,f5,f6,xb,xa
        real(real64) :: t1, t2, t3, t4, t5, t6
        real(real64), dimension(size(x)) :: x1,x2,x3,x4,x5,x6
        
        real(real64), parameter :: c20=0.25, c21=0.25, &        
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
   
        ! Prepare time steps
        t1 = t
        t2 = t + c20*h
        t3 = t + c30*h
        t4 = t + c40*h
        t5 = t +     h
        t6 = t + c60*h
        
        !f1 = h*f(t1 , x) 
        !f2 = h*f(t2 , x + c21*f1)  
        !f3 = h*f(t3 , x + c31*f1 + c32*f2)   
        !f4 = h*f(t4 , x + c41*f1 + c42*f2 + c43*f3)  
        !f5 = h*f(t5 , x + c51*f1 + c52*f2 + c53*f3 + c54*f4)       
        !f6 = h*f(t6 , x + c61*f1 + c62*f2 + c63*f3 + c64*f4 + c65*f5)
        
        x1 = x
        
        ! Stgage 1
        f1 = h*f(t1 , x1) 
        x2 = x + c21*f1
        
        ! Stgage 2
        f2 = h*f(t2 , x2)  
        x3 = x + c31*f1 + c32*f2
        
        ! Stgage 3
        f3 = h*f(t3 , x3)   
        x4 = x + c41*f1 + c42*f2 + c43*f3
        
        ! Stgage 4
        f4 = h*f(t4 , x4)  
        x5 = x + c51*f1 + c52*f2 + c53*f3 + c54*f4
        
        ! Stgage 5
        f5 = h*f(t5 , x5)       
        x6 = x + c61*f1 + c62*f2 + c63*f3 + c64*f4 + c65*f5
        
        ! Stgage 6
        f6 = h*f(t6 , x6)
        
        xb = x1 + b1*f1 + b3*f3 + b4*f4 + b5*f5 + b6*f6      
        xa = x1 + a1*f1 + a3*f3 + a4*f4 + a5*f5    
        
        t = t + h   
        x = xa
        e = maxval( abs(xa - xb) )
    end subroutine rk45_vec

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
        real(real64), intent(in):: tb,emin,emax,hmin,hmax     
        integer, intent(out):: iflag
        integer, intent(in):: itmax
        real(real64), intent(inout) :: h, t, x
        real(real64) ::delta, d,xsave,tsave,e
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
    
    !----------------------------------------------------------------------------------
    ! RK45AD_VEC - Adaptive scheme based on Runge-Kutta-Fehlberg method
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
        real(real64), intent(in):: tb,emin,emax,hmin,hmax     
        integer, intent(out):: iflag
        integer, intent(in):: itmax
        real(real64), intent(inout) :: h, t
        real(real64), intent(inout) :: x(:)
        real(real64) ::delta, d,tsave,e
        real(real64), dimension(size(x)) :: xsave
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