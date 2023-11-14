    module mod_spatial_vectors
    use mod_array_inv
    implicit none
    
    !integer, parameter :: wp = real64
    
    type vector2
        real(wp) :: data(2)
    end type
    type vector3
        real(wp) :: data(3)
    contains
        procedure :: vec3_cross_op
        generic :: cross=> vec3_cross_op
    end type
    type matrix2
        real(wp) :: data(2,2)
    end type
    type matrix3
        real(wp) :: data(3,3)
    end type
            
    type vector
        integer :: size
        real(wp), allocatable :: data(:)
    contains    
        procedure :: item => vec_item, vec_items
        procedure :: write_vector
        generic :: write(formatted) => write_vector
    end type vector

    type matrix
        integer :: rows, columns
        real(wp), allocatable :: data(:,:)
    contains
        procedure :: item => mat_item, mat_items
        procedure :: write_matrix
        generic :: write(formatted) => write_matrix
    end type matrix
    
    type lu_matrix
        integer :: size
        real(wp), allocatable :: data(:,:)
        integer, allocatable :: indx(:)
        integer :: sgn
        integer :: ierr
    contains
        procedure :: solve => lu_solve_vec, lu_solve_mat
    end type lu_matrix
    
    interface vector
        module procedure vec_from_array, vec_zeros, vec_elemental
    end interface
    interface matrix
        module procedure mat_from_array_1, mat_zeros, mat_elemental
        module procedure mat_from_array_2
    end interface
    
    interface vector2
        module procedure :: vec2_zero
    end interface
    interface matrix2
        module procedure :: mat2_zero
    end interface
    interface vector3
        module procedure :: vec3_zero
    end interface
    interface matrix3
        module procedure :: mat3_zero
    end interface
    
    interface operator (+)
        module procedure add_vec_vec
        module procedure add_vec2_vec2
        module procedure add_vec3_vec3
        module procedure add_mat_mat
        module procedure add_mat2_mat2
        module procedure add_mat3_mat3
    end interface
    
    interface operator (-)
        module procedure sub_vec_vec
        module procedure sub_vec2_vec2
        module procedure sub_vec3_vec3
        module procedure sub_mat_mat
        module procedure sub_mat2_mat2
        module procedure sub_mat3_mat3
    end interface
    
    interface operator (*)
        module procedure mul_vec_scalar, mul_scalar_vec
        module procedure mul_vec_array, mul_array_vec
        module procedure mul_vec_mat, mul_mat_vec
        module procedure mul_mat_mat, mul_scalar_mat, mul_mat_scalar
    end interface
    
    interface operator (/)
        module procedure div_vec_scalar 
        module procedure div_vec2_scalar
        module procedure div_vec3_scalar
        module procedure div_mat_scalar 
        module procedure div_mat2_scalar
        module procedure div_mat3_scalar
    end interface

    interface assignment (=)
        module procedure array_to_vec , vec_to_array
        module procedure array_to_vec2, vec2_to_array
        module procedure array_to_vec3, vec3_to_array
        module procedure array_to_mat , mat_to_array
        module procedure array_to_mat2, mat2_to_array
        module procedure array_to_mat3, mat3_to_array
    end interface
    
    interface operator (.i.)
        module procedure vec_dot_vec
        module procedure vec2_dot_vec2
        module procedure vec3_dot_vec3
    end interface
    interface operator (.x.)
        module procedure s_cross_vec2
        module procedure vec2_cross_s
        module procedure vec2_cross_vec2
        module procedure vec3_cross_vec3
    end interface
    interface operator (.o.)
        module procedure vec_outer_vec
        module procedure vec2_outer_vec2
        module procedure vec3_outer_vec3
    end interface
    
    interface cross
        module procedure s_cross_vec2
        module procedure vec2_cross_s
        module procedure vec2_cross_vec2
        module procedure vec3_cross_vec3
        module procedure vec3_cross_op
    end interface
    
    interface zeros
        module procedure :: mat_zeros
    end interface
    interface ident
        module procedure :: mat_identity
    end interface
    interface rand
        module procedure :: vec_random_unit, vec_random
        module procedure :: mat_random_unit, mat_random
    end interface
    
    interface det
        module procedure :: det_mat, det_mat2, det_mat3
    end interface
    interface inv
        module procedure :: inv_mat, inv_mat2, inv_mat3
    end interface
    interface solve
        module procedure :: solve_mat2_vec2, solve_mat3_vec3
        module procedure :: solve_mat2_mat2, solve_mat3_mat3
        module procedure :: solve_lu_vec, solve_lu_mat
    end interface
    interface lu
        module procedure :: lu_decompose
    end interface

    contains
    
    !! Assignments
            
    pure function vec_from_array(a) result(v)
    real(wp), intent(in) :: a(:)
    type(vector) :: v
        v = vector(size(a), a)
    end function
    
    pure function mat_from_array_1(n, m, a) result(g)
    integer, intent(in) :: n, m
    real(wp), intent(in) :: a(:)
    type(matrix) :: g
        g = matrix(n,m, reshape(a, [n,m]))
    end function
    
    pure function mat_from_array_2(a) result(g)
    real(wp), intent(in) :: a(:,:)
    type(matrix) :: g
        g%rows = size(a,1)
        g%columns = size(a,2)
        g%data = a
    end function
    
    pure subroutine array_to_vec(v,a)
    type(vector), intent(out) :: v
    real(wp), intent(in) :: a(:)
    integer :: n
        v%size = size(a)
        v%data = a
    end subroutine
    
    pure subroutine vec_to_array(a,v)
    real(wp), intent(out), allocatable :: a(:)
    type(vector), intent(in) :: v
        allocate(a(v%size))
        a(:) = v%data(:)
    end subroutine
    
    pure subroutine array_to_vec2(v,a)
    type(vector2), intent(out) :: v
    real(wp), intent(in) :: a(2)
        v%data = a
    end subroutine
    
    pure subroutine vec2_to_array(a,v)
    real(wp), intent(out) :: a(2)
    type(vector2), intent(in) :: v
        a = v%data
    end subroutine
    
    pure subroutine array_to_vec3(v,a)
    type(vector3), intent(out) :: v
    real(wp), intent(in) :: a(3)
        v%data = a
    end subroutine

    pure subroutine vec3_to_array(a,v)
    real(wp), intent(out) :: a(3)
    type(vector3), intent(in) :: v
        a = v%data
    end subroutine
    
    pure subroutine array_to_mat(mx,a)
    type(matrix), intent(out) :: mx
    real(wp), intent(in) :: a(:,:)
    integer :: i, n , m
        n = size(a,1)
        m = size(a,2)
        mx%rows = n
        mx%columns = m
        mx%data = reshape(a, [n,m])
    end subroutine
    
    pure subroutine mat_to_array(a, mx)
    real(wp), intent(out), allocatable :: a(:,:)
    type(matrix), intent(in) :: mx
    integer :: i, n, m
        n = mx%rows
        m = mx%columns
        a = mx%data
    end subroutine
    
    pure subroutine array_to_mat2(mx,a)
    type(matrix2), intent(out) :: mx
    real(wp), intent(in) :: a(2,2)
        mx%data = a
    end subroutine
    
    pure subroutine mat2_to_array(a, mx)
    real(wp), intent(out) :: a(2,2)
    type(matrix2), intent(in) :: mx
        a = mx%data
    end subroutine

    pure subroutine array_to_mat3(mx,a)
    type(matrix3), intent(out) :: mx
    real(wp), intent(in) :: a(3,3)
        mx%data = a
    end subroutine

    pure subroutine mat3_to_array(a, mx)
    real(wp), intent(out) :: a(3,3)
    type(matrix3), intent(in) :: mx
        a = mx%data
    end subroutine
    
    pure function vec2_zero() result(x)
    type(vector2) :: x
        x%data = 0.0_wp
    end function
    pure function vec3_zero() result(x)
    type(vector3) :: x
        x%data = 0.0_wp
    end function
    pure function mat2_zero() result(x)
    type(matrix2) :: x
        x%data = 0.0_wp
    end function
    pure function mat3_zero() result(x)
    type(matrix3) :: x
        x%data = 0.0_wp
    end function
    
    pure function vec_zeros(n,e) result(x)
    type(vector) :: x
    integer, intent(in) :: n
    real(wp), optional, intent(in) :: e
        x%size = n
        allocate(x%data(n))        
        if( present(e) ) then
            x%data(1:n) = e
        else
            x%data(1:n) = 0._wp
        endif
    end function
    pure function vec_elemental(n,i) result(x)
    type(vector) :: x
    integer, intent(in) :: n, i
        x = vec_zeros(n)
        if(i>= 1 .and. i<=n) then
            x%data(i) = 1._wp
        end if
    end function
    
    function vec_random_unit(n) result(x)
    type(vector) :: x
    integer, intent(in) :: n
        x%size = n
        allocate(x%data(n))
        call random_number(x%data)
    end function
    
    function vec_random(n,x_min,x_max) result(x)
    type(vector) :: x
    integer, intent(in) :: n
    real(wp), intent(in) :: x_min, x_max
        x%size = n
        allocate(x%data(n))
        call random_number(x%data)
        x%data = x_min + (x_max - x_min) * x%data        
    end function
    
    pure function mat_zeros(n,m,e) result(x)
    type(matrix) :: x
    integer, intent(in) :: n, m
    real(wp), optional, intent(in) :: e
        x%rows = n
        x%columns = m
        allocate(x%data(n,m))
        if( present(e) ) then
            x%data(1:n,1:m) = e
        else
            x%data(1:n,1:m) = 0._wp
        endif
    end function
    pure function mat_elemental(n,m,i,j) result(x)
    type(matrix) :: x
    integer, intent(in) :: n, m, i, j
        x = mat_zeros(n,m)
        if(i>= 1 .and. i<=n .and. j>= 1 .and. j<=m) then
            x%data(i,j) = 1._wp
        end if        
    end function
    function mat_random_unit(n,m) result(x)
    type(matrix) :: x
    integer, intent(in) :: n, m
        x%rows = n
        x%columns = m
        allocate(x%data(n,m))
        call random_number(x%data)
    end function
    
    function mat_random(n,m,x_min,x_max) result(x)
    type(matrix) :: x
    integer, intent(in) :: n, m
    real(wp), intent(in) :: x_min, x_max
        x%rows = n
        x%columns = m
        allocate(x%data(n,m))
        call random_number(x%data)
        x%data = x_min + (x_max - x_min) * x%data        
    end function
    
    pure function mat_identity(n,m) result(x)
    type(matrix) :: x
    integer, intent(in) :: n
    integer, intent(in), optional :: m
    integer :: i
        if( present(m) ) then
            x = mat_zeros(n,m)
        else
            x = mat_zeros(n,n)
        end if
        forall(i=1:min(n,m))
            x%data(i,i) = 1._wp
        end forall
    end function
    
    pure function vec_item(a,i) result(x)
    real(wp) :: x
    class(vector), intent(in) :: a
    integer, intent(in) :: i
        x = a%data(i)
    end function
    pure function vec_items(a,i_start,i_end) result(x)
    real(wp), allocatable :: x(:)
    class(vector), intent(in) :: a
    integer, intent(in) :: i_start,i_end
        x = a%data(i_start:i_end)
    end function
    
    pure function mat_item(a,i,j) result(x)
    real(wp) :: x
    class(matrix), intent(in) :: a
    integer, intent(in) :: i,j
        x = a%data(i,j)
    end function
    pure function mat_items(a,i_start,j_start,i_end,j_end) result(x)
    real(wp), allocatable :: x(:,:)
    class(matrix), intent(in) :: a
    integer, intent(in) :: i_start, j_start, i_end, j_end
        x = a%data(i_start:i_end, j_start:j_end)
    end function

    !! Algebra
    
    elemental function add_vec_vec(a, b) result(r)
        type(vector), intent(in) :: a,b
        type(vector) :: r
        r%size = a%size
        r%data = a%data + b%data
    end function
    elemental function sub_vec_vec(a, b) result(r)
        type(vector), intent(in) :: a,b
        type(vector) :: r
        r%size = a%size
        r%data = a%data - b%data
    end function
    elemental function add_vec2_vec2(a, b) result(r)
        type(vector2), intent(in) :: a,b
        type(vector2) :: r
        r%data = a%data + b%data
    end function
    elemental function sub_vec2_vec2(a, b) result(r)
        type(vector2), intent(in) :: a,b
        type(vector2) :: r
        r%data = a%data - b%data
    end function
    elemental function add_vec3_vec3(a, b) result(r)
        type(vector3), intent(in) :: a,b
        type(vector3) :: r
        r%data = a%data + b%data
    end function
    elemental function sub_vec3_vec3(a, b) result(r)
        type(vector3), intent(in) :: a,b
        type(vector3) :: r
        r%data = a%data - b%data
    end function
    
    elemental function add_mat_mat(a, b) result(r)
        type(matrix), intent(in) :: a,b
        type(matrix) :: r
        r%rows = a%rows
        r%columns = a%columns
        r%data = a%data + b%data
    end function
    elemental function sub_mat_mat(a, b) result(r)
        type(matrix), intent(in) :: a,b
        type(matrix) :: r
        r%rows = a%rows
        r%columns = a%columns
        r%data = a%data - b%data
    end function
    
    elemental function add_mat2_mat2(a, b) result(r)
        type(matrix2), intent(in) :: a,b
        type(matrix2) :: r
        r%data = a%data + b%data
    end function
    elemental function sub_mat2_mat2(a, b) result(r)
        type(matrix2), intent(in) :: a,b
        type(matrix2) :: r
        r%data = a%data - b%data
    end function
    elemental function add_mat3_mat3(a, b) result(r)
        type(matrix3), intent(in) :: a,b
        type(matrix3) :: r
        r%data = a%data + b%data
    end function
    elemental function sub_mat3_mat3(a, b) result(r)
        type(matrix3), intent(in) :: a,b
        type(matrix3) :: r
        r%data = a%data - b%data
    end function
    
    elemental function mul_vec_scalar(vt,f) result(ut)
    type(vector), intent(in) :: vt
    real(wp), intent(in) :: f
    type(vector) :: ut
        ut % size = vt % size        
        ut % data = vt%data * f
    end function    
        
    elemental function mul_scalar_vec(f, vt) result(ut)
    real(wp), intent(in) :: f
    type(vector), intent(in) :: vt
    type(vector) :: ut
        ut % size = vt % size        
        ut % data = f * vt%data
    end function    
    
    pure function mul_vec_array(vt,a) result(ut)
    type(vector), intent(in) :: vt
    real(wp), intent(in) :: a(:,:)
    type(vector) :: ut
        ut % size = size(a,2)
        ut % data = matmul(transpose(a), vt%data)
    end function
    
    
    pure function mul_array_vec(a, vt) result(ut)
    real(wp), intent(in) :: a(:,:)
    type(vector), intent(in) :: vt
    type(vector) :: ut
        ut%size = size(a,1)
        ut%data = matmul(a, vt%data)
    end function
    
    
    pure function mul_vec_mat(vt, mx) result(ut)
    type(vector), intent(in) :: vt
    type(matrix), intent(in) :: mx
    type(vector) :: ut
        ut%size = mx%columns
        ut%data = matmul(transpose(mx%data), vt%data)
    end function
    
    pure function mul_mat_vec(mx, vt) result(ut)
    type(matrix), intent(in) :: mx
    type(vector), intent(in) :: vt
    type(vector) :: ut
        ut%size = mx%rows
        ut%data = matmul(mx%data, vt%data)
    end function
    
    pure function mul_mat_mat(mx, mz) result(mt)
    type(matrix), intent(in) :: mx, mz
    type(matrix) :: mt
        mt%rows = mx%rows
        mt%columns = mz%columns
        mt%data = matmul(mx%data, mz%data)
    end function    
    
    elemental function mul_scalar_mat(f, a) result(b)
    real(wp), intent(in) :: f
    type(matrix), intent(in) :: a
    type(matrix) :: b
        b = matrix( f*a%data )
    end function
    elemental function mul_mat_scalar(a, f) result(b)
    type(matrix), intent(in) :: a
    real(wp), intent(in) :: f
    type(matrix) :: b
        b = matrix( f*a%data )
    end function
    
    elemental function div_vec_scalar(v,d) result(u)
    type(vector), intent(in) :: v
    real(wp), intent(in) :: d
    type(vector) :: u
        u%size = v%size
        u%data = (1/d)*v%data
    end function
    
    elemental function div_vec2_scalar(v,d) result(u)
    type(vector2), intent(in) :: v
    real(wp), intent(in) :: d
    type(vector2) :: u
        u%data = (1/d)*v%data
    end function
    elemental function div_vec3_scalar(v,d) result(u)
    type(vector3), intent(in) :: v
    real(wp), intent(in) :: d
    type(vector3) :: u
        u%data = (1/d)*v%data
    end function
    
    elemental function div_mat_scalar(m,d) result(u)
    type(matrix), intent(in) :: m
    real(wp), intent(in) :: d
    type(matrix) :: u
        u%rows = m%rows
        u%columns = m%columns
        u%data = (1/d)*m%data
    end function
    
    elemental function div_mat2_scalar(m,d) result(u)
    type(matrix2), intent(in) :: m
    real(wp), intent(in) :: d
    type(matrix2) :: u
        u%data = (1/d)*m%data
    end function
    
    elemental function div_mat3_scalar(m,d) result(u)
    type(matrix3), intent(in) :: m
    real(wp), intent(in) :: d
    type(matrix3) :: u
        u%data = (1/d)*m%data
    end function
    
    pure function det_mat(a) result(d)
    type(matrix), intent(in) :: a
    real(wp) :: d
    type(lu_matrix) :: lu
    integer :: i, n
        lu = lu_decompose(a)
        d = lu%sgn
        n = lu%size
        do i=1, n
            d = d * lu%data(i,i)
        end do
    end function
    pure function det_mat2(a) result(d)
    type(matrix2), intent(in) :: a
    real(wp) :: d
        d = mat2_det(a%data)
    end function
    pure function det_mat3(a) result(d)
    type(matrix3), intent(in) :: a
    real(wp) :: d
        d = mat3_det(a%data)
    end function
    function inv_mat(a) result(b)
    type(matrix), intent(in) :: a
    type(matrix) :: b
        b%data = inv(a%data)
    end function
    pure function inv_mat2(a) result(b)
    type(matrix2), intent(in) :: a
    type(matrix2) :: b
        b%data = mat2_inv(a%data)
    end function
    pure function inv_mat3(a) result(b)
    type(matrix3), intent(in) :: a
    type(matrix3) :: b
        b%data = mat3_inv(a%data)
    end function
    
    function solve_mat_vec(a,b) result(x)
    type(matrix), intent(in) :: a
    type(vector), intent(in) :: b
    type(vector) :: x        
        x = vector( mat_solve_vec(a%data, b%data))
    end function
    pure function solve_mat2_vec2(a,b) result(x)
    type(matrix2), intent(in) :: a
    type(vector2), intent(in) :: b
    type(vector2) :: x
        x = vector2( mat2_solve_vec(a%data, b%data) )
    end function
    pure function solve_mat3_vec3(a,b) result(x)
    type(matrix3), intent(in) :: a
    type(vector3), intent(in) :: b
    type(vector3) :: x
        x = vector3( mat3_solve_vec(a%data, b%data) )
    end function
    
    function solve_mat_mat(a,b) result(x)
    type(matrix), intent(in) :: a, b
    type(matrix) :: x        
        x = matrix( mat_solve_mat(a%data, b%data))
    end function
    pure function solve_mat2_mat2(a,b) result(x)
    type(matrix2), intent(in) :: a,b
    type(matrix2) :: x
        x = matrix2( matmul( mat2_inv(a%data), b%data) )
    end function
    pure function solve_mat3_mat3(a,b) result(x)
    type(matrix3), intent(in) :: a, b
    type(matrix3) :: x
        x = matrix3( matmul( mat3_inv(a%data), b%data) )
    end function
    
    pure function lu_decompose(a) result(lu)
    !  ***************************************************************
    !  * Given an N x N matrix A, this routine replaces it by the LU *
    !  * decomposition of a rowwise permutation of itself. A and N   *
    !  * are input. INDX is an output vector which records the row   *
    !  * permutation effected by the partial pivoting; D is output   *
    !  * as -1 or 1, depending on whether the number of row inter-   *
    !  * changes was even or odd, respectively. This routine is used *
    !  * in combination with LUBKSB to solve linear equations or to  *
    !  * invert a matrix. Return code is 1, if matrix is singular.   *
    !  ***************************************************************
    type(matrix), intent(in) :: a
    type(lu_matrix) :: lu
    real(wp)  :: amax, dum, summ
    real(wp), allocatable  :: vv(:)
    integer :: i, j, k, imax, n , m, sgn, ierr
        
        n = a%rows
        m = a%columns
        if(n/=m) then
            error stop "Expecting square matrix for LU decomposition."
        end if
        lu%size = a%rows
        lu%data = a%data
        allocate(lu%indx(n))
        allocate(vv(n))
        lu%sgn=1; lu%ierr=0

        do i=1,n
            amax=0_wp
            do j=1,n
                if (abs(lu%data(i,j)) > amax) amax=abs(lu%data(i,j))
            end do ! j loop
            if(amax < tiny) then
                lu%ierr = 1
                return
            end if
            vv(i) = 1 / amax
        end do ! i loop

        do j=1,n
            do i=1,j-1
                summ = lu%data(i,j)
                do k=1,i-1
                    summ = summ - lu%data(i,k)*lu%data(k,j)
                end do ! k loop
                lu%data(i,j) = summ
            end do ! i loop
            amax = 0._wp
            do i=j,n
                summ = lu%data(i,j)
                do k=1,j-1
                    summ = summ - lu%data(i,k)*lu%data(k,j)
                end do ! k loop
                lu%data(i,j) = summ
                dum = vv(i)*abs(summ)
                if(dum >= amax) then
                    imax = i
                    amax = dum
                end if
            end do ! i loop

            if(j /= imax) then
                do k=1,n
                    dum = lu%data(imax,k)
                    lu%data(imax,k) = lu%data(j,k)
                    lu%data(j,k) = dum
                end do ! k loop
                lu%sgn = -lu%sgn
                vv(imax) = vv(j)
            end if

            lu%indx(j) = imax
            if(abs(lu%data(j,j)) < tiny) lu%data(j,j) = tiny

            if(j /= n) then
                dum = 1 / lu%data(j,j)
                do i=j+1,n
                    lu%data(i,j) = lu%data(i,j)*dum
                end do ! i loop
            end if
        end do ! j loop
    
    end function
    
    pure function lu_solve_vec(lu,b) result(x)
    !  ******************************************************************
    !  * Solves the set of N linear equations A . X = B.  Here A is     *
    !  * input, not as the matrix A but rather as its LU decomposition, *
    !  * determined by the routine LUDCMP. INDX is input as the permuta-*
    !  * tion vector returned by LUDCMP. B is input as the right-hand   *
    !  * side vector B, and returns with the solution vector X. A, N and*
    !  * INDX are not modified by this routine and can be used for suc- *
    !  * cessive calls with different right-hand sides. This routine is *
    !  * also efficient for plain matrix inversion.                     *
    !  ******************************************************************
    class(lu_matrix), intent(in) :: lu
    type(vector), intent(in) :: b
    type(vector) :: x
        
    integer :: i,j,ii,ll,n
    real(wp) :: summ
    
    n = lu%size
    
    if( n/=b%size) then
        error stop "Incompatible vector size."
    end if

    ii = 0
    x = b
    do i=1,n
        ll = lu%indx(i)
        summ = x%data(ll)
        x%data(ll) = x%data(i)
        if(ii /= 0) then
            do j=ii,i-1
                summ = summ - lu%data(i,j)*x%data(j)
            end do ! j loop
        else if(summ /= 0.d0) then
            ii = i
        end if
        x%data(i) = summ
    end do ! i loop

    do i=n,1,-1
        summ = x%data(i)
        if(i < n) then
            do j=i+1,n
                summ = summ - lu%data(i,j)*x%data(j)
            end do ! j loop
        end if
        x%data(i) = summ / lu%data(i,i)
    end do ! i loop
    
    end function
    
    pure function lu_solve_mat(lu,b) result(x)
    class(lu_matrix), intent(in) :: lu
    type(matrix), intent(in) :: b
    type(matrix) :: x
        
    integer :: i,j,ii,ll,n,m,k
    real(wp) :: summ
    
    n = lu%size
    m = b%columns
    
    if( n/=b%rows) then
        error stop "Incompatible matrix size."
    end if

    ii = 0
    x = b
    
    do k=1,m
        do i=1,n
            ll = lu%indx(i)
            summ = x%data(ll,k)
            x%data(ll,k) = x%data(i,k)
            if(ii /= 0) then
                do j=ii,i-1
                    summ = summ - lu%data(i,j)*x%data(j,k)
                end do ! j loop
            else if(summ /= 0.d0) then
                ii = i
            end if
            x%data(i,k) = summ
        end do ! i loop

        do i=n,1,-1
            summ = x%data(i,k)
            if(i < n) then
                do j=i+1,n
                    summ = summ - lu%data(i,j)*x%data(j,k)
                end do ! j loop
            end if
            x%data(i,k) = summ / lu%data(i,i)
        end do ! i loop
    end do ! k loop
    end function
    
    pure function solve_lu_vec(a,b) result(x)
    type(matrix), intent(in) :: a
    type(vector), intent(in) :: b
    type(vector) :: x
    type(lu_matrix) :: lu
        lu = lu_decompose(a)
        x = lu_solve_vec(lu, b)
    end function
    
    pure function solve_lu_mat(a,b) result(x)
    type(matrix), intent(in) :: a
    type(matrix), intent(in) :: b
    type(matrix) :: x
    type(lu_matrix) :: lu
        lu = lu_decompose(a)
        x = lu_solve_mat(lu, b)
    end function
    
    
    pure function vec_dot_vec(a,b) result(s)
    type(vector), intent(in) :: a, b
    real(wp) :: s
        s = dot_product(a%data, b%data)
    end function
    pure function vec2_dot_vec2(a,b) result(s)
    type(vector2), intent(in) :: a, b
    real(wp) :: s
        s = dot_product(a%data, b%data)
    end function
    pure function vec3_dot_vec3(a,b) result(s)
    type(vector3), intent(in) :: a, b
    real(wp) :: s
        s = dot_product(a%data, b%data)
    end function
    
    pure function vec_outer_vec(a,b) result(r)
    type(vector), intent(in) :: a, b
    type(matrix) :: r
    integer :: i, j, n, m
        n = a%size
        m = b%size
        r%rows = n
        r%columns = m
        allocate(r%data(n,m))
        forall (i=1:n)
          forall(j=1:m) r%data(i,j) = a%data(i)*b%data(j)
        end forall        
    end function
    
    pure function vec2_outer_vec2(a,b) result(r)
    type(vector2), intent(in) :: a, b
    type(matrix2) :: r
    integer :: i, j
        forall (i=1:2)
          forall(j=1:2) r%data(i,j) = a%data(i)*b%data(j)
        end forall        
    end function
    
    pure function vec3_outer_vec3(a,b) result(r)
    type(vector3), intent(in) :: a, b
    type(matrix3) :: r
    integer :: i, j
        forall (i=1:3)
          forall(j=1:3) r%data(i,j) = a%data(i)*b%data(j)
        end forall        
    end function
                
    pure function vec2_cross_vec2(a,b) result(c)
    type(vector2), intent(in) ::a
    type(vector2), intent(in) ::b
    real(wp) :: c
        c = a%data(1)*b%data(2)-a%data(2)*b%data(1)
    end function

    pure function vec2_cross_s(a,b) result(c)
    type(vector2), intent(in) ::a
    real(wp), intent(in) ::b
    type(vector2) :: c                
        c%data = [ &
            a%data(2)*b, &
            -a%data(1)*b ]        
    end function
    
    pure function s_cross_vec2(a,b) result(c)
    real(wp), intent(in) ::a
    type(vector2), intent(in) ::b
    type(vector2) :: c                
        c%data = [ &
            -a*b%data(2), &
            a*b%data(2) ]        
    end function
    
    pure function vec3_cross_vec3(a,b) result(c)
    type(vector3), intent(in) ::a, b
    type(vector3) :: c                
        c%data = [ &
            a%data(2)*b%data(3)-a%data(3)*b%data(2), &
            a%data(3)*b%data(1)-a%data(1)*b%data(3), &
            a%data(1)*b%data(2)-a%data(2)*b%data(1) ]        
    end function

    pure function vec3_cross_op(a) result(c)
    class(vector3), intent(in) ::a
    type(matrix3) :: c
    real(wp) :: x,y,z
    real(wp), parameter :: O = 0
    
        !tex:Constructs the 3×3 skew symmetric cross product
        ! operator from a vector $\boldsymbol{v}=\pmatrix{x&y&z}$
        ! $$ \mathrm{cr}(\boldsymbol{v})=\pmatrix{0 & -z & y \\ z & 0 & -x \\ -y & x & 0}$$    
    
        x = a%data(1)
        y = a%data(2)
        z = a%data(3)
        c%data = reshape([O,z,-y, -z,O,x, y,-x,O], [3,3])    
    end function
    
    !! Display

    subroutine write_vector(v, unit, iotype, v_list, iostat, iomsg)
    class(vector), intent(in) :: v
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    integer :: i, n
    character(LEN=32) :: fmt, item
    if( iotype(1:2)=='DT' ) then
        fmt = iotype(3:)
    else
        fmt = "g0.4"
    end if
    if (allocated(v%data)) then
        n = size(v%data)
        write (unit, '(a)',iostat=iostat, iomsg=iomsg) '['
        do i = 1, n            
            write (item, "("//trim(fmt)//")", iostat=iostat, iomsg=iomsg)  &
                v%data(i)
            write(unit,*) trim(item)
            if (iostat /= 0) return
            if( i<n ) then
                write (unit, '(a)') ','
            end if
        end do
        write (unit, '(a)', iostat=iostat, iomsg=iomsg) ']'
        !write (unit, "(/)", iostat=iostat, iomsg=iomsg)
    else
        write (unit, "(l1,/)", iostat=iostat, iomsg=iomsg) .false.
    end if
    end subroutine 
    
    subroutine write_matrix(mx, unit, iotype, v_list, iostat, iomsg)
    class(matrix), intent(in) :: mx
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg
    integer :: i, j, n, m
    character(LEN=32) :: fmt, item
    if( iotype(1:2)=='DT' ) then
        fmt = iotype(3:)
    else
        fmt = "g0.4"
    end if
    if (allocated(mx%data)) then
        n = mx%rows
        m = mx%columns
        write (unit, '(a)',iostat=iostat, iomsg=iomsg) '['
        do i = 1, n
            write (unit, '(a)',iostat=iostat, iomsg=iomsg) '['
            do j = 1, m
                write (item, "("//trim(fmt)//")", iostat=iostat, iomsg=iomsg)  &
                    mx%data(i, j)
                write(unit,*) trim(item)
                if (iostat /= 0) return
                if( j<m ) then
                    write (unit, '(a)') ','
                end if
            end do
            write (unit, '(a)', iostat=iostat, iomsg=iomsg) ']'
            if (iostat /= 0) return
            if( i<n ) then
                write (unit, '(a)') ','
            end if
        end do
        write (unit, '(a)', iostat=iostat, iomsg=iomsg) ']'
        !write (unit, "(/)", iostat=iostat, iomsg=iomsg)
    else
        write (unit, "(l1,/)", iostat=iostat, iomsg=iomsg) .false.
    end if
    end subroutine 
    

    end module
    