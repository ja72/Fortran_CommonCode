    module mod_native_vectors
    use mod_common
    implicit none
    
    integer, parameter :: wp = real64
    
    enum, bind(c)
        enumerator :: origin = 0
        enumerator :: x_axis = 1
        enumerator :: y_axis = 2
        enumerator :: z_axis = 3
    end enum
    
    real(wp), dimension(3), parameter :: o_ = [0._wp,0._wp,0._wp]
    real(wp), dimension(3), parameter :: i_ = [1._wp,0._wp,0._wp]
    real(wp), dimension(3), parameter :: j_ = [0._wp,1._wp,0._wp]
    real(wp), dimension(3), parameter :: k_ = [0._wp,0._wp,1._wp]
    real(wp), dimension(4), parameter :: q_eye = [1._wp, 0._wp, 0._wp, 0._wp]
    
    real(wp), dimension(3,3), parameter :: zero_ = reshape( [0._wp,0._wp,0._wp, 0._wp,0._wp,0._wp, 0._wp,0._wp,0._wp], [3,3])
    real(wp), dimension(3,3), parameter :: eye_ = reshape( [1._wp,0._wp,0._wp, 0._wp,1._wp,0._wp, 0._wp,0._wp,1._wp], [3,3])

    interface vector
        module procedure vector_from_axis
    end interface
    
    interface sumsq
        module procedure vector_sumsq
    end interface
    interface norm
        module procedure vector_norm
    end interface
    interface unit
        module procedure vector_unit
    end interface
    
    interface operator (.x.)
        module procedure vector_cross_product, vector_cross_product_matrix
    end interface

    interface mmoi
        module procedure vector_mmoi_matrix
    end interface
    
    interface dot
        module procedure vector_dot_product
    end interface
        
    interface cross
        module procedure vector_cross_product_matrix
        module procedure vector_cross_product
    end interface
    
    interface outer
        module procedure vector_outer_product
    end interface        
        
    interface rotx
        module procedure vector_rotate_x, matrix_rotate_x
    end interface 
    interface roty
        module procedure vector_rotate_y, matrix_rotate_y
    end interface 
    interface rotz
        module procedure vector_rotate_z, matrix_rotate_z
    end interface 
    
    interface rot
        module procedure vector_rotate_axis, matrix_rotate_axis
        module procedure vector_rotate_vector, matrix_rotate_vector
        module procedure quat_rot_matrix, quat_rot_vector
    end interface
    
    interface quat
        module procedure quat_from_vector, quat_axis_angle, quat_vector_angle
    end interface
    
    contains    
    
    function vector_from_axis(axis) result(v)
    integer, intent(in) :: axis
    real(wp) :: v(3)
    
        select case (axis)
        case (x_axis)
            v = i_
        case (y_axis)
            v = j_
        case (z_axis)
            v = k_
        case default
            v = o_
        end select    
    end function
    
    function vector_sumsq(a) result(m)
    real(wp), intent(in) :: a(3)
    real(wp) :: m
        m = a(1)**2+a(2)**2+a(3)**2
    end function
    
    function vector_norm(a) result(m)
    real(wp), intent(in) :: a(3)
    real(wp) :: m
        m = sqrt( a(1)**2+a(2)**2+a(3)**2 )
    end function
    
    function vector_unit(a) result(b)
    real(wp), intent(in) :: a(3)
    real(wp) :: b(3), m
        m = norm(a)
        if( m>0 .and. m /= 1.0 ) then
            b = a/m
        else
            b = a        
        end if
    end function
    
    function vector_project_vector(vector,axis) result(res)
    real(wp),  intent(in) :: vector(3), axis(3)
    real(wp) :: res(3), t
         t = dot_product(axis,vector) / dot_product(axis,axis)
         res = axis * t         
    end function
    
    function vector_perpendicular_vector(vector,axis) result(res)
    real(wp),  intent(in) :: vector(3), axis(3)
    real(wp) :: res(3)
         res = vector - vector_project_vector(vector, axis)
    end function
    
    function vector_project_axis(vector,axis) result(res)
    real(wp),  intent(in) :: vector(3)
    integer, intent(in) :: axis
    real(wp) :: res(3)
        select case(axis)
        case (x_axis)
            res = [ vector(1), 0._wp, 0._wp]
        case (y_axis)
            res = [ 0._wp, vector(2), 0._wp]
        case (z_axis)
            res = [ 0._wp, 0._wp, vector(3)]
        case default
            res = o_
        end select
    end function
    
    function vector_perpendicular_axis(vector,axis) result(res)
    real(wp),  intent(in) :: vector(3)
    integer, intent(in) :: axis
    real(wp) :: res(3)
        res = vector - vector_project_axis(vector, axis)
    end function
    
    pure function vector_dot_product(a,b) result(c)
    real(wp), intent(in) :: a(:), b(:)
    real(wp) :: c
        c = dot_product(a,b)
    end function
    
    pure function vector_outer_product(a,b) result(c)
    real(wp), intent(in) :: a(:), b(:)
    real(wp) :: c(size(a),size(b))
        integer :: i,j,n,m
        n = size(a)
        m = size(b)
        forall(i=1:n, j=1:m)
            c(i,j) = a(i)*b(j)
        end forall
    end function
    
    pure function vector_cross_product(a,b) result(c)
    real(wp),  intent(in) :: a(3), b(3)
    real(wp) ::  c(3)
    
        c = [ a(2)*b(3) - a(3)*b(2), &
              a(3)*b(1) - a(1)*b(3), &
              a(1)*b(2) - a(2)*b(1) ]
    end function
    
    pure function vector_cross_product_matrix(a) result(c)
    real(wp),  intent(in) :: a(3)
    real(wp) ::  c(3,3)
    
        c = reshape( &
            [0._wp, a(3), -a(2), &
            -a(3), 0._wp, a(1), &
            a(2), -a(1), 0._wp], [3,3])
        
    end function
        
    function vector_mmoi_matrix(a) result(c)
    real(wp),  intent(in) :: a(3)
    real(wp) ::  c(3,3)
    
        c = reshape( &
            [   a(2)**2+a(3)**2, -a(1)*a(2), -a(1)*a(3), &
                -a(1)*a(2), a(1)**2+a(3)**2, -a(2)*a(3), &
                -a(1)*a(3), -a(2)*a(3), a(1)**2+a(2)**2  ], [3,3])
        
    end function
    
    function vector_rotate_x(vector, angle) result(res)
    real(wp),  intent(in) :: vector(3), angle
    real(wp) ::  res(3), c, s
        c = cos(angle)
        s = sin(angle)
        res = [ &
            vector(1), &
            vector(2)*c - vector(3)*s, &
            vector(2)*s + vector(3)*c]        
    end function
    
    function matrix_rotate_x(angle) result(res)
    real(wp),  intent(in) :: angle
    real(wp) ::  res(3,3), c, s
        c = cos(angle)
        s = sin(angle)
        res = reshape( &
            [1._wp, 0._wp, 0._wp, &
            0._wp, c,s, &
            0._wp, -s,c], [3,3] )        
    end function
    
    function vector_rotate_y(vector, angle) result(res)
    real(wp),  intent(in) :: vector(3), angle
    real(wp) ::  res(3), c, s
        c = cos(angle)
        s = sin(angle)
        res = [ &
            vector(1)*c + vector(3)*s, &
            vector(2), &
            -vector(1)*s + vector(3)*c]        
    end function
    
    function matrix_rotate_y(angle) result(res)
    real(wp),  intent(in) :: angle
    real(wp) ::  res(3,3), c, s
        c = cos(angle)
        s = sin(angle)
        res = reshape( &
            [ c, 0._wp, -s, &
            0._wp, 1._wp, 0._wp, &
            s, 0._wp, c], [3,3] )        
    end function
    
    function vector_rotate_z(vector, angle) result(res)
    real(wp),  intent(in) :: vector(3), angle
    real(wp) ::  res(3), c, s
        c = cos(angle)
        s = sin(angle)
        res = [ &
            vector(1)*c - vector(2)*s, &
            vector(1)*s + vector(2)*c, &
            vector(3)]        
    end function
    
    function matrix_rotate_z(angle) result(res)
    real(wp),  intent(in) :: angle
    real(wp) ::  res(3,3), c, s
        c = cos(angle)
        s = sin(angle)
        res = reshape( &
            [c,s,0._wp, &
            -s,c,0._wp, &
            0._wp,0._wp,1._wp], [3,3] )        
    end function
    
    function vector_rotate_axis(vector, axis, angle) result(res)
    integer, intent(in) :: axis
    real(wp),  intent(in) :: vector(3), angle
    real(wp) ::  res(3)
        select case(axis)
        case (x_axis)
            res = vector_rotate_x(vector, angle)
        case (y_axis)
            res = vector_rotate_y(vector, angle)
        case (z_axis)
            res = vector_rotate_z(vector, angle)
        case default
            res = o_
        end select
    end function
    
    function matrix_rotate_axis(axis, angle) result(res)
    integer, intent(in) :: axis
    real(wp),  intent(in) :: angle
    real(wp) ::  res(3,3)
        select case(axis)
        case (x_axis)
            res = matrix_rotate_x(angle)
        case (y_axis)
            res = matrix_rotate_y(angle)
        case (z_axis)
            res = matrix_rotate_z(angle)
        case default
            res = zero_
        end select
    end function
    
    function vector_rotate_vector(vector, axis, angle) result(res)
    real(wp),  intent(in) :: vector(3), axis(3), angle
    real(wp) ::  a(3), res(3), axv(3), axaxv(3), c, s
        c = cos(angle)
        s = sin(angle)
        a = unit(axis)
        axv = cross(a, vector)
        axaxv = cross(a, axv)
        res = vector + s*axv + (1d0-c)*axaxv
    end function
    
    function matrix_rotate_vector(axis, angle) result(res)
    real(wp),  intent(in) :: axis(3), angle
    real(wp) ::  a(3), res(3,3), ax(3,3), axax(3,3),c,s
        c = cos(angle)
        s = sin(angle)
        a = unit(axis)
        ax = cross(a)
        axax = matmul(ax,ax)
        res = eye_ + s*ax + (1d0-c)*axax
    end function
    
    pure function quat_axis_angle(axis, angle) result(q)
    integer, intent(in) :: axis
    real(wp),  intent(in) :: angle
    real(wp) ::  q(4)
        select case(axis)
        case(x_axis)
            q = quat_vector_angle(i_, angle)
        case(y_axis)
            q = quat_vector_angle(j_, angle)
        case(z_axis)
            q = quat_vector_angle(k_, angle)
        case default
            q = q_eye    
        end select
    end function
    
    pure function quat_vector_angle(axis, angle) result(q)
    real(wp),  intent(in) :: axis(3), angle
    real(wp) ::  q(4), c, s, m
        c = cos(angle/2)
        s = sin(angle/2)    
        m = norm2(axis)
        q = [c, (s/m)*axis]
    end function

    pure function quat_from_vector(vector) result(q)
    real(wp),  intent(in) :: vector(3)
    real(wp) ::  q(4)
        q = [0._wp, vector]
    end function
    
    pure function quat_rot_matrix(q, inv) result(R)
    real(wp),  intent(in) :: q(4)
    logical, intent(in), optional :: inv
    real(wp) ::  q_s, q_v(3), R(3,3), qx(3,3), qxqx(3,3)
    
        q_s = q(1)
        q_v = q(2:4)
        qx = cross(q_v)
        qxqx = matmul(qx,qx)
        
        if( present(inv) .and. inv) then
            q_s = -q_s
        end if
        
        R = eye_ + 2*q_s*qx + 2*qxqx
    end function
    
    pure function quat_rot_vector(q, v, inv) result(r)
    real(wp),  intent(in) :: q(4), v(3)
    logical, intent(in), optional :: inv
    real(wp) ::  q_s, q_v(3), r(3), qxv(3), qxqxv(3)
        q_s = q(1)
        q_v = q(2:4)
        qxv = cross(q_v,v)
        qxqxv = cross(q_v, qxv)
        
        if( present(inv) .and. inv) then
            q_s = -q_s
        end if
        
        r = v + 2*q_s*qxv + 2*qxqxv
    end function
    
    pure function quat_conjugate(q) result(p)
    real(wp), intent(in) :: q(4)
    real(wp) :: p(4)
        p = [ q(1), -q(2), -q(3), -q(4) ]
    end function
    
    pure function quat_product(q_1,q_2) result(p)
    real(wp), intent(in) :: q_1(4), q_2(4)
    real(wp) :: p(4), s_1, s_2, v_1(3), v_2(3)
        s_1 = q_1(1)
        s_2 = q_2(1)
        v_1 = q_1(2:4)
        v_2 = q_2(2:4)
        
        p = [s_1*s_2 - dot_product(v_1,v_2), &
            s_1*v_2 + s_2*v_1 + cross(v_1,v_2)]
        
    end function
    
    pure function quat_cross_product(q_1,q_2) result(p)
    real(wp), intent(in) :: q_1(4), q_2(4)
    real(wp) :: p(4), v_1(3), v_2(3)
        v_1 = q_1(2:4)
        v_2 = q_2(2:4)        
        p = [0.0_wp, &
            cross(v_1,v_2)]        
    end function

    pure function quat_dot_product(q_1,q_2) result(p)
    real(wp), intent(in) :: q_1(4), q_2(4)
    real(wp) :: p, s_1, s_2, v_1(3), v_2(3)
        s_1 = q_1(1)
        s_2 = q_2(1)
        v_1 = q_1(2:4)
        v_2 = q_2(2:4)
        
        p = s_1*s_2 + dot_product(v_1,v_2)
        
    end function
    
    pure function quat_magnitude(q) result(m)
    real(wp),  intent(in) :: q(4)
    real(wp) :: m    
         m = sqrt(dot_product(q,q))
    end function

    pure function quat_normalize(q) result(p)
    real(wp),  intent(in) :: q(4)
    real(wp) :: p(4), m2
        m2 = dot_product(q, q)
        if( m2 >= 0.0_wp) then
            p = q/sqrt(m2)
        else
            p = q
        end if
    end function
    
    pure function quat_inv(q) result(p)
    real(wp),  intent(in) :: q(4)
    real(wp) :: p(4), m2
        m2 = dot_product(q, q)
        if( m2 /= 1 .and. m2/=0) then
            p = [ q(1)/m2, -q(2)/m2, -q(3)/m2, -q(4)/m2 ]
        else if(m2 /= 0) then
            p = [ q(1), -q(2), -q(3), -q(4) ]
        else
            error stop "Cannot invert a zero quaternion."
        end if
    end function
    
    function vector_angle(a,b) result(t)
    real(wp), intent(in) :: a(3), b(3)
    reaL(wp) :: t, ma, mb, ab
    
        ! |a.b| = |a| |b| cos(t)
        ! |a×b| = |a| |b| sin(t)
        ma = norm(a)
        mb = norm(b)
    
        if( ma == 0._wp .or. mb == 0._wp ) then
            t = 0._wp
            return
        end if
    
        ab = dot_product(a, b)
        t = acos(ab/(ma*mb))
        
    end function
    
    end module