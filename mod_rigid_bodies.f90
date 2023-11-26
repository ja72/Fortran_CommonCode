    module mod_rigid_bodies
    use mod_native_vectors
    
    integer, parameter :: n_dims=3, n_state=13
    
    type :: rigid_body
        real(wp) :: mass
        real(wp) :: mmoi(3)
        real(wp) :: cg(3)
    contains
        procedure :: get_mmoi => rb_get_mmoi_state
        procedure :: get_motion => rb_get_motion
        procedure :: get_kin_energy => rb_get_kin_energy
    end type
    
    
    contains
    
    pure function rb_ellipsoid(mass, dx,dy,dz, cg) result(rb)
    type(rigid_body) :: rb
    real(wp), intent(in) :: mass, dx,dy,dz
    real(wp), optional, intent(in) :: cg(3)
    real(wp), parameter :: f= 1.0_wp/10
    
        rb%mass = mass
        if(present(cg)) then
            rb%cg = cg
        else
            rb%cg = o_
        end if
        
        rb%mmoi = f*mass*[ (dy**2+dz**2), (dz**2+dx**2), (dx**2+dy**2)]
    
    end function
    
    pure function rb_box(mass, dx,dy,dz, cg) result(rb)
    type(rigid_body) :: rb
    real(wp), intent(in) :: mass, dx,dy,dz
    real(wp), optional, intent(in) :: cg(3)
    real(wp), parameter :: f= 1.0_wp/12
    
        rb%mass = mass
        if(present(cg)) then
            rb%cg = cg
        else
            rb%cg = o_
        end if
        
        rb%mmoi = f*mass*[ (dy**2+dz**2), (dz**2+dx**2), (dx**2+dy**2)]
    
    end function
    
    pure function rb_cylinder(mass, diameter, height, cg) result(rb)
    type(rigid_body) :: rb
    real(wp), intent(in) :: mass, diameter, height
    real(wp), optional, intent(in) :: cg(3)
    real(wp), parameter :: fh= 1.0_wp/10, fd=1.0_wp/16
    
        rb%mass = mass
        if(present(cg)) then
            rb%cg = cg
        else
            rb%cg = o_
        end if
        
        rb%mmoi = mass*[ fh*height**2 + fd*diameter**2, fh*height**2 + fd*diameter**2, 2*fd*diameter**2]
    
    end function
    pure function rb_disk(mass, diameter, cg) result(rb)
    type(rigid_body) :: rb
    real(wp), intent(in) :: mass, diameter
    real(wp), optional, intent(in) :: cg(3)
        if(present(cg)) then
            rb = rb_cylinder(mass, diameter, 0.0_wp, cg)
        else
            rb = rb_cylinder(mass, diameter, 0.0_wp)
        end if
    end function
    pure function rb_rod(mass, length, cg) result(rb)
    type(rigid_body) :: rb
    real(wp), intent(in) :: mass, length
    real(wp), optional, intent(in) :: cg(3)
        if(present(cg)) then
            rb = rb_cylinder(mass, 0.0_wp, length, cg)
        else
            rb = rb_cylinder(mass, 0.0_wp, length)
        end if
    end function
    
    pure function rb_get_mmoi_rot(rb, R, inv) result(I)
    ! Get mass moment of inertia tensor for rigid body
    ! given a the rotation matrix
    class(rigid_body),intent(in) :: rb
    real(wp), intent(in) :: R(3,3)
    logical, intent(in), optional :: inv
    real(wp) :: I(3,3)
    integer :: n
        if( present(inv) .and. inv) then
            forall(n=1:3)
                I(n,:) = (1/rb%mmoi(n)) * R(:,n)
            end forall
        else
            forall(n=1:3)
                I(n,:) = rb%mmoi(n) * R(:,n)
            end forall
        end if
        I = matmul(R, I)
    end function
    
    pure function rb_get_mmoi_state(rb, state, inv) result(I)
    class(rigid_body),intent(in) :: rb
    real(wp), intent(in) :: state(n_state)
    logical, intent(in), optional :: inv
    real(wp) :: I(3,3), R(3,3), q(4)
        q = state(4:7)
        R = rot(q)
        if(present(inv)) then
            I = rb_get_mmoi_rot(rb, R, inv)
        else
            I = rb_get_mmoi_rot(rb, R)
        end if
    end function
    
    pure function rb_get_state(rb, r_A, q, v_A, omega) result(state)
    real(wp) :: state(n_state)
    class(rigid_body),intent(in) :: rb
    real(wp), intent(in) :: r_A(3), q(4), v_A(3), omega(3)
    real(wp) :: R(3,3), I_C(3,3), p(3), L_A(3), c(3)
    !tex: State vector composition
    !$$\boldsymbol{y}:\begin{Bmatrix}\boldsymbol{r}_{A}=(x_{A},y_{A},z_{A})\\
    !\boldsymbol{q}=(\boldsymbol{\hat{z}}\sin\left(\tfrac{\theta}{2}\right)|\cos\left(\tfrac{\theta}{2}\right))\\
    !\boldsymbol{p}=m\left(\boldsymbol{v}_{A}+\boldsymbol{\omega}\times\boldsymbol{c}\right)\\
    !\boldsymbol{L}_{A}={\rm I}_{C}\boldsymbol{\omega}+\boldsymbol{c}\times\boldsymbol{p}
    !\end{Bmatrix}$$
        R = rot(q)
        I_C = rb_get_mmoi_rot(rb, R)
        c = matmul(R, rb%cg)
        p = rb%mass*(v_A + cross(omega,c))
        L_A = matmul(I_C, omega) + cross(c,p)
        state = [r_A, q, p, L_A]
    end function
    
    pure subroutine rb_get_motion(rb, state, c, v_A, omega) 
    ! Get motion variables from rigid body state y=(r_A,q,p,L_A)
    class(rigid_body),intent(in) :: rb
    real(wp), intent(in) :: state(n_state)
    real(wp), intent(out) :: c(3), v_A(3), omega(3)
    real(wp) :: R(3,3), I_inv(3,3), q(4), p(3), L_A(3)
    !tex: Rigid body motion derived from momentum
    !$$\begin{Bmatrix}\boldsymbol{v}_{A}=\tfrac{1}{m}\boldsymbol{p}-\boldsymbol{\omega}\times\boldsymbol{c}\\
    !\boldsymbol{\omega}={\rm I}_{C}^{-1}\left(\boldsymbol{L}_{A}-\boldsymbol{c}\times\boldsymbol{p}\right)
    !\end{Bmatrix}$$
    
        ! quaternion for orientation
        q = state(4:7)                              
        ! 3×3 rotation matrix from quaternion
        R = rot(q)
        ! relative position of cg
        c = matmul(R, rb%cg)     
        ! 3×3 inverse interia matrix at the CG 
        ! from principal mass moment of inertia and rotation matrix
        I_inv = rb_get_mmoi_rot(rb, R, .true.)
        ! linear momentum vector
        p = state(8:10)
        ! angular momentum vector at the handle
        L_A = state(11:13)
        ! angular velocity from angular momentum
        omega = matmul(I_inv, L_A - cross(c, p))
        ! linear velocity of the handle from momentum
        v_A = p / rb%mass - cross(omega, c)        
    end subroutine
    
    pure function rb_get_state_rate(rb, state, force, torque_A) result(rate)
    class(rigid_body),intent(in) :: rb
    real(wp), intent(in) :: state(n_state)
    real(wp), intent(in) :: force(3), torque_A(3)
    real(wp) :: rate(n_state), p(3), L_A(3), c(3)
    real(wp) :: r_A(3), q(4), v_A(3), omega(3)
    real(wp) :: rp_A(3), qp(4), pp(3), Lp_A(3)
    !tex: State vector rate of change
    !$$\tfrac{{\rm d}}{{\rm d}t}\boldsymbol{y}:\begin{Bmatrix}\dot{\boldsymbol{r}}_{A}=\tfrac{1}{m}\boldsymbol{p}-\boldsymbol{\omega}\times\boldsymbol{c}\\
    !\dot{\boldsymbol{q}}=\tfrac{1}{2}(\boldsymbol{q}_{s}\boldsymbol{\omega}+\boldsymbol{\omega}\times\boldsymbol{q}_{v}|-\boldsymbol{\omega}\cdot\boldsymbol{q}_{v})\\
    !\dot{\boldsymbol{p}}=\boldsymbol{F}\\
    !\dot{\boldsymbol{L}}_{A}=\boldsymbol{\tau}_{A}+\boldsymbol{p}\times\boldsymbol{v}_{A}
    !\end{Bmatrix}$$
        r_A = state(1:3)
        q = state(4:7)
        p = state(8:10)
        L_A = state(11:13)
        call rb%get_motion(state, c, v_A, omega)
        rp_A = p/rb%mass - cross(omega,c)
        qp = quat_product( [0.0_wp, omega/2], q)
        pp = force
        Lp_A = torque_A + cross(p, v_A)
        rate = [rp_A, qp, pp, Lp_A]
    end function
    
    pure function rb_get_kin_energy(rb, state) result(KE)
    class(rigid_body),intent(in) :: rb
    real(wp), intent(in) :: state(n_state)
    real(wp) :: KE, c(3), v_A(3), omega(3), p(3), L_A(3)
    !tex: Kinetic Energy 
    !$$KE=\tfrac{1}{2}\left(\boldsymbol{v}_{A}\cdot\boldsymbol{p}+\boldsymbol{\omega}\cdot\boldsymbol{L}_{A}\right)$$
        call rb%get_motion(state, c, v_A, omega)
        p = state(8:10)
        L_A = state(11:13)
        KE = (dot_product(v_A,p) + dot_product(omega,L_A))/2
    end function
    
    subroutine test_rb_sim()
    use mod_show_matrix
    use mod_nasa_ode
    type(rigid_body) :: rb
    real(wp) :: r_A(3), q(4), v_A(3), omega(3), c(3), I_C(3,3)
    real(wp) :: KE, omg_max, h_min, p(3), L_A(3), F(3), tau(3)
    real(wp), dimension(n_state) :: state, rate
    integer :: i, n, flag
    real(wp) :: t,t_end, t_next, h, y(n_state), yp(n_state)
    real ( kind = 8 ) abserr
    real ( kind = 8 ) relerr
    
    real(wp), parameter :: gravity_(3) = 0.0_wp*j_
    real(wp), parameter :: dia = 1, ht=0.125, m = 1
    real(wp), parameter :: cg(3) = o_
        
        rb = rb_cylinder(m, dia, ht, cg)
        print *, "RB MASS=", rb%mass
        print *, "RB MMOI="
        call show(rb%mmoi)
        print *, "RB CG="
        call show(rb%cg)
        
        r_A = o_
        print *, "r_A="
        call show(r_A)
        q = quat(i_, pi/2)
        print *, "q="
        call show(q)
        omega = 10*k_
        print *, "omega="
        call show(omega)
        v_A = cross(r_A + cg, omega)
        print *, "v_A="
        call show(v_A)  ! -5*j_
        
        state = rb_get_state(rb, r_A, q, v_A, omega)
        !print *, "state="
        !call show(state)
        
        I_C = rb%get_mmoi(state)
        print *, "I_C="
        call show(I_C)
        
        call rb%get_motion(state, c, v_A, omega)
        print *, "c="
        call show(c)
        p = state(8:10)
        print *, "p="
        call show(p)
        L_A = state(11:13)
        print *, "L_A="
        call show(L_A)
        
        KE = rb_get_kin_energy(rb, state)
        
        print *, "KE (state)=", KE
        q = state(4:7)
        print *, "q (mag) = ", quat_magnitude(q), NEW_LINE('a')
        
        omg_max = maxval(abs(omega))
        h_min = 5/(deg*omg_max)
        print *, "Ideal Time Step = ", real(h_min), NEW_LINE('a')
        tau = o_
        F = o_        
                
        t = 0.0_wp
        t_end = 10.0_wp
        y = state
        yp = rb_get_state_rate(rb, y, F, tau)
        abserr = sqrt ( epsilon ( t_end ) )
        relerr = sqrt ( epsilon ( t_end ) )
        flag = 1
        write (*, '(a4,1x,a8,1x,*(a12,1x))') "i", "t", "KE", "omg_x", "omg_y", "omg_z", "L_A_x", "L_A_y", "L_A__z", "h"
        call rb_get_motion(rb, y, c, v_A, omega)
        KE = rb_get_kin_energy(rb, y)
        L_A = y(11:13)
        write (*, '(i4,1x,f8.3,1x,*(f12.6,1x))') 0, t, KE, omega, L_A
        i = 1
        do while(t<=t_end .and. ( abs(flag)==1 .or. abs(flag)==2))
            q = y(4:7)
            y(4:7) = quat_normalize(q)
            t_next = t+45*h_min
            call r8_rkf45(rb_rate, n_state, y, yp, t, t_next, relerr, abserr, flag, h)            
            call rb_get_motion(rb, y, c, v_A, omega)
            KE = rb_get_kin_energy(rb, y)
            L_A = y(11:13)
            write (*, '(i4,1x,f8.3,1x,*(f12.6,1x))') i, t, KE, omega, L_A, h
            i = i + 1
        end do
        
    contains
    
        pure subroutine rb_rate(t, n, y, yp)
        real (wp), intent(in) :: t
        integer, intent(in) :: n
        real (wp), intent(in) :: y(n)
        real (wp), intent(out) :: yp(n)
        integer :: n_bodies, i, k
        real(wp), dimension(n_state) :: state, rate
            
            n_bodies = n/n_state
            
            do i=1, n_bodies
                k = 1 + n_state * (i-1)
                state(1:n_state) = y(k:k+n_state-1)
                rate = rb_get_state_rate(rb, state, F, tau)
                yp(k:k+n_state-1) = rate(1:n_state)
            end do
        
        end subroutine
    
    end subroutine
    
    end module