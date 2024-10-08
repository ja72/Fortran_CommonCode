module mod_array_inv
use mod_common
use mod_show_matrix
implicit none

    integer, parameter :: wp = real64
    
    type, public :: t_matrix(rows, cols, knd)
       integer, len :: rows, cols
       integer, kind :: knd = real64
       real(kind=knd), dimension(rows, cols) :: data
    contains
        procedure :: solve=> t_mat_solve
        procedure :: inv=> t_mat_invert
    end type    
    
    
    interface det
        module procedure :: mat_det
    end interface
    interface inv
        module procedure :: mat_inv
    end interface
    interface solve
        module procedure :: mat_solve_vec, mat_solve_mat
    end interface

    contains
    
    pure function vec_inner(a,b) result(ab)
    real(real64), intent(in) :: a(:), b(:)
    real(real64) :: ab
        ab = dot_product(a,b)
    end function
    
    pure function vec_outer(a,b) result(ab)
    real(real64), intent(in) :: a(:), b(:)
    real(real64), allocatable :: ab(:,:)
    integer :: n,m,i,j
        n = size(a)
        m = size(b)
        allocate(ab(n,m))
        forall(i=1:n, j=1:m)
            ab(i,j) = a(i)*b(j)
        end forall
    end function
        
    function mat_det(A) result(d)
    real(real64), intent(in) :: A(:,:)
    real(real64) :: d
    integer :: n, m
    
        n = size(A, 1)
        m = size(A, 2)
        
        if(n /= m) then
            error stop "Expecting square matrix."
        end if
        
        select case(n)
        case(1)  
            d = A(1,1)
        case(2)
            d = mat2_det(A)
        case(3)
            d = mat3_det(A)
        case(4)
            d = mat4_det(A)
        case default
            d = lu_mat_det(A)
        end select
    end function
    
    function mat_inv(A) result(B)
    real(real64), intent(in) :: A(:,:)    
    real(real64) :: B(size(A,1),size(A,2))
    integer :: n, m
    
        n = size(A, 1)
        m = size(A, 2)
        
        if(n /= m) then
            error stop "Expecting square matrix."
        end if
        
        select case(n)
        case (1)
            B = 1/A
        case (2)
            B = mat2_inv(A)
        case (3)
            B = mat3_inv(A)
        case (4)
            B = mat4_inv(A)
        case default
            B = lu_mat_invert(A)
        end select
    end function
    
    function mat_solve_vec(A, b) result(x)
    real(real64), intent(in) :: A(:,:), b(:)
    real(real64), allocatable :: x(:)
    integer :: n,m
    
        n = size(A, 1)
        m = size(A, 2)
        
        if(n /= m) then
            error stop "Expecting square matrix."
        end if

        select case(n)
        case (1)
            x = b/A(1,1)
        case (2)   
            x = mat2_solve_vec(A,b)
        case (3)
            x = mat3_solve_vec(A,b)
        case (4)
            x = mat4_solve_vec(A,b)
        case default
            x = lu_mat_solve_vec(A,b)
        end select
    end function
    
    function mat_solve_mat(A, B) result(X)
    real(real64), intent(in) :: A(:,:), B(:,:)
    real(real64), allocatable :: X(:,:)
    real(real64), allocatable :: A_inv(:,:)
        A_inv = mat_inv(A)
        X = matmul(A_inv, B)
    end function
    
    pure function mat2_det(A) result(d)
    real(real64) :: d    
    real(real64), intent(in) :: A(2,2)
    real(real64) :: t2, t5
        t2 = A(1,1)*A(2,2)
        t5 = A(1,2)*A(2,1)
        d = t2-t5
    end function
    
    pure function mat2_inv(A) result(B)
    real(real64) :: B(2,2), d_inv
    real(real64), intent(in) :: A(2,2)        
        d_inv = 1/mat2_det(A)
        B(1,1) = A(2,2)*d_inv
        B(1,2) = -A(1,2)*d_inv
        B(2,1) = -A(2,1)*d_inv
        B(2,2) = A(1,1)*d_inv        
    end function
    
    pure function mat2_solve_vec(A,b) result(x)
    real(real64) :: x(2), d_inv
    real(real64), intent(in) :: A(2,2), b(2)
        d_inv = 1/mat2_det(A)
        x(1) = d_inv*(A(2,2)*b(1) - A(1,2)*b(2))
        x(2) = d_inv*(A(1,1)*b(2) - A(2,1)*b(1))    
    end function
    
    pure function mat3_det(A) result(d)
    real(real64) :: d
    real(real64), intent(in) :: A(3,3)
    real(real64) :: t2, t3, t4, t7, t8, t9
        t2 = A(1,1)*A(2,2)*A(3,3)
        t3 = A(1,2)*A(2,3)*A(3,1)
        t4 = A(1,3)*A(2,1)*A(3,2)
        t7 = A(1,1)*A(2,3)*A(3,2)
        t8 = A(1,2)*A(2,1)*A(3,3)
        t9 = A(1,3)*A(2,2)*A(3,1)
        d = t2+t3+t4-t7-t8-t9
    end function
    
    pure function mat3_inv(A) result(B)
    real(real64) :: B(3,3), d_inv
    real(real64), intent(in) :: A(3,3)
    
        d_inv = 1/mat3_det(A)
        B(1,1) = d_inv*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
        B(1,2) = -d_inv*(A(1,2)*A(3,3)-A(1,3)*A(3,2))
        B(1,3) = d_inv*(A(1,2)*A(2,3)-A(1,3)*A(2,2))
        B(2,1) = -d_inv*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
        B(2,2) = d_inv*(A(1,1)*A(3,3)-A(1,3)*A(3,1))
        B(2,3) = -d_inv*(A(1,1)*A(2,3)-A(1,3)*A(2,1))
        B(3,1) = d_inv*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
        B(3,2) = -d_inv*(A(1,1)*A(3,2)-A(1,2)*A(3,1))
        B(3,3) = d_inv*(A(1,1)*A(2,2)-A(1,2)*A(2,1))    
        
    end function
    
    pure function mat3_solve_vec(A,b) result(x)
    real(real64) :: x(3), d_inv
    real(real64), intent(in) :: A(3,3), b(3)
        d_inv = 1/mat3_det(A)
        
        x(1) = d_inv*(A(1,2)*(A(2,3)*b(3)-A(3,3)*b(2))+A(1,3)*(A(3,2)*b(2)-A(2,2)*b(3))+b(1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)))
        x(2) = d_inv*(A(1,1)*(A(3,3)*b(2)-A(2,3)*b(3))+A(1,3)*(A(2,1)*b(3)-A(3,1)*b(2))-b(1)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
        x(3) = d_inv*(A(1,1)*(A(2,2)*b(3)-A(3,2)*b(2))+A(1,2)*(A(3,1)*b(2)-A(2,1)*b(3))+b(1)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    end function
    
    pure function mat4_det(A) result(d)
    implicit real(real64) (T)
    real(real64) :: d
    real(real64), intent(in) :: A(4,4)
    real(real64) :: t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13
    real(real64) :: t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27
      t2 = A(1,1)*A(2,2)*A(3,3)*A(4,4)
      t3 = A(1,1)*A(2,3)*A(3,4)*A(4,2)
      t4 = A(1,1)*A(2,4)*A(3,2)*A(4,3)
      t5 = A(1,2)*A(2,1)*A(3,4)*A(4,3)
      t6 = A(1,2)*A(2,3)*A(3,1)*A(4,4)
      t7 = A(1,2)*A(2,4)*A(3,3)*A(4,1)
      t8 = A(1,3)*A(2,1)*A(3,2)*A(4,4)
      t9 = A(1,3)*A(2,2)*A(3,4)*A(4,1)
      t10 = A(1,3)*A(2,4)*A(3,1)*A(4,2)
      t11 = A(1,4)*A(2,1)*A(3,3)*A(4,2)
      t12 = A(1,4)*A(2,2)*A(3,1)*A(4,3)
      t13 = A(1,4)*A(2,3)*A(3,2)*A(4,1)
      t16 = A(1,1)*A(2,2)*A(3,4)*A(4,3)
      t17 = A(1,1)*A(2,3)*A(3,2)*A(4,4)
      t18 = A(1,1)*A(2,4)*A(3,3)*A(4,2)
      t19 = A(1,2)*A(2,1)*A(3,3)*A(4,4)
      t20 = A(1,2)*A(2,3)*A(3,4)*A(4,1)
      t21 = A(1,2)*A(2,4)*A(3,1)*A(4,3)
      t22 = A(1,3)*A(2,1)*A(3,4)*A(4,2)
      t23 = A(1,3)*A(2,2)*A(3,1)*A(4,4)
      t24 = A(1,3)*A(2,4)*A(3,2)*A(4,1)
      t25 = A(1,4)*A(2,1)*A(3,2)*A(4,3)
      t26 = A(1,4)*A(2,2)*A(3,3)*A(4,1)
      t27 = A(1,4)*A(2,3)*A(3,1)*A(4,2)
      d = t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13-t16-t17-t18-t19-t20-t21-t22-t23-t24-t25-t26-t27
    end function
    
    
    pure function mat4_inv(A) result(B)
    real(real64) :: B(4,4), d_inv
    real(real64), intent(in) :: A(4,4)
    
      d_inv = 1/mat4_det(A)
      B(1,1) = d_inv*(A(2,2)*A(3,3)*A(4,4)-A(2,2)*A(3,4)*A(4,3)-A(2,3)*A(3,2)*A(4,4)+A(2,3)*A(3,4)*A(4,2)+A(2,4)*A(3,2)*A(4,3)-A(2,4)*A(3,3)*A(4,2))
      B(1,2) = -d_inv*(A(1,2)*A(3,3)*A(4,4)-A(1,2)*A(3,4)*A(4,3)-A(1,3)*A(3,2)*A(4,4)+A(1,3)*A(3,4)*A(4,2)+A(1,4)*A(3,2)*A(4,3)-A(1,4)*A(3,3)*A(4,2))
      B(1,3) = d_inv*(A(1,2)*A(2,3)*A(4,4)-A(1,2)*A(2,4)*A(4,3)-A(1,3)*A(2,2)*A(4,4)+A(1,3)*A(2,4)*A(4,2)+A(1,4)*A(2,2)*A(4,3)-A(1,4)*A(2,3)*A(4,2))
      B(1,4) = -d_inv*(A(1,2)*A(2,3)*A(3,4)-A(1,2)*A(2,4)*A(3,3)-A(1,3)*A(2,2)*A(3,4)+A(1,3)*A(2,4)*A(3,2)+A(1,4)*A(2,2)*A(3,3)-A(1,4)*A(2,3)*A(3,2))
      B(2,1) = -d_inv*(A(2,1)*A(3,3)*A(4,4)-A(2,1)*A(3,4)*A(4,3)-A(2,3)*A(3,1)*A(4,4)+A(2,3)*A(3,4)*A(4,1)+A(2,4)*A(3,1)*A(4,3)-A(2,4)*A(3,3)*A(4,1))
      B(2,2) = d_inv*(A(1,1)*A(3,3)*A(4,4)-A(1,1)*A(3,4)*A(4,3)-A(1,3)*A(3,1)*A(4,4)+A(1,3)*A(3,4)*A(4,1)+A(1,4)*A(3,1)*A(4,3)-A(1,4)*A(3,3)*A(4,1))
      B(2,3) = -d_inv*(A(1,1)*A(2,3)*A(4,4)-A(1,1)*A(2,4)*A(4,3)-A(1,3)*A(2,1)*A(4,4)+A(1,3)*A(2,4)*A(4,1)+A(1,4)*A(2,1)*A(4,3)-A(1,4)*A(2,3)*A(4,1))
      B(2,4) = d_inv*(A(1,1)*A(2,3)*A(3,4)-A(1,1)*A(2,4)*A(3,3)-A(1,3)*A(2,1)*A(3,4)+A(1,3)*A(2,4)*A(3,1)+A(1,4)*A(2,1)*A(3,3)-A(1,4)*A(2,3)*A(3,1))
      B(3,1) = d_inv*(A(2,1)*A(3,2)*A(4,4)-A(2,1)*A(3,4)*A(4,2)-A(2,2)*A(3,1)*A(4,4)+A(2,2)*A(3,4)*A(4,1)+A(2,4)*A(3,1)*A(4,2)-A(2,4)*A(3,2)*A(4,1))
      B(3,2) = -d_inv*(A(1,1)*A(3,2)*A(4,4)-A(1,1)*A(3,4)*A(4,2)-A(1,2)*A(3,1)*A(4,4)+A(1,2)*A(3,4)*A(4,1)+A(1,4)*A(3,1)*A(4,2)-A(1,4)*A(3,2)*A(4,1))
      B(3,3) = d_inv*(A(1,1)*A(2,2)*A(4,4)-A(1,1)*A(2,4)*A(4,2)-A(1,2)*A(2,1)*A(4,4)+A(1,2)*A(2,4)*A(4,1)+A(1,4)*A(2,1)*A(4,2)-A(1,4)*A(2,2)*A(4,1))
      B(3,4) = -d_inv*(A(1,1)*A(2,2)*A(3,4)-A(1,1)*A(2,4)*A(3,2)-A(1,2)*A(2,1)*A(3,4)+A(1,2)*A(2,4)*A(3,1)+A(1,4)*A(2,1)*A(3,2)-A(1,4)*A(2,2)*A(3,1))
      B(4,1) = -d_inv*(A(2,1)*A(3,2)*A(4,3)-A(2,1)*A(3,3)*A(4,2)-A(2,2)*A(3,1)*A(4,3)+A(2,2)*A(3,3)*A(4,1)+A(2,3)*A(3,1)*A(4,2)-A(2,3)*A(3,2)*A(4,1))
      B(4,2) = d_inv*(A(1,1)*A(3,2)*A(4,3)-A(1,1)*A(3,3)*A(4,2)-A(1,2)*A(3,1)*A(4,3)+A(1,2)*A(3,3)*A(4,1)+A(1,3)*A(3,1)*A(4,2)-A(1,3)*A(3,2)*A(4,1))
      B(4,3) = -d_inv*(A(1,1)*A(2,2)*A(4,3)-A(1,1)*A(2,3)*A(4,2)-A(1,2)*A(2,1)*A(4,3)+A(1,2)*A(2,3)*A(4,1)+A(1,3)*A(2,1)*A(4,2)-A(1,3)*A(2,2)*A(4,1))
      B(4,4) = d_inv*(A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(1,2)*A(2,1)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1))

    end function

    pure function mat4_solve_vec(A,b) result(x)
    real(real64) :: x(4)
    real(real64), intent(in) :: A(4,4), b(4)
        x = matmul(mat4_inv(A), b)
    end function
    
    function lu_mat_det(A) result(d)
    use mod_lu
    real(real64) :: d
    real(real64), intent(in) :: A(:,:)
    logical :: ok
    real(real64), allocatable :: temp(:), LU(:,:)
    integer, allocatable :: indx(:)

    integer :: i, rc, n, m
    
        n = size(A,1)
        m = size(A,2)
        if( n/= m) then
            error stop "Expecting a square matrix."
        end if

        ok = .false.
        allocate(LU(n,n))
        allocate(temp(n+1))
        allocate(INDX(n))
        LU = A
        !call LU decomposition routine
        call LUDCMP(LU,n,INDX,D,rc)
        do i=1, n
            d =d * LU(i,i)
        end do    
    end function
    
    pure function lu_mat_invert(A) result(A_inv)
    use mod_lu
    real(real64), allocatable :: A_inv(:,:)
    real(real64), intent(in) :: A(:,:)
    logical :: ok
    real(real64), allocatable :: temp(:), LU(:,:)
    integer, allocatable :: indx(:)
    real(real64) :: d
    integer :: i, j, rc, n, m
    
        n = size(A,1)
        m = size(A,2)
        if( n/= m) then
            error stop "Expecting a square matrix."
        end if
        allocate(A_inv(n,n))
        A_inv = 0.0_wp
        forall(i=1:n)
            A_inv(i,i) = 1.0_wp
        end forall
        
        ok = .false.
        allocate(LU(n,n))
        allocate(temp(n+1))
        allocate(INDX(n))
        LU = A
        !call LU decomposition routine
        call LUDCMP(LU,n,INDX,D,rc)

        !call appropriate solver if previous return code is ok
        if (rc == 0) then
            do j=1, n
                call LUBKSB(LU,n,INDX, A_inv(:,j))
            end do
            ok = .true.
        endif
        
    end function
    
    pure function lu_mat_solve_vec(A, b) result(x)
    use mod_lu
    real(real64), allocatable :: x(:)
    real(real64), intent(in) :: A(:,:), b(:)
    logical :: ok
    real(real64), allocatable :: temp(:), LU(:,:)
    integer, allocatable :: indx(:)
    real(real64) :: d
    integer :: rc, n, m
    
        n = size(A,1)
        m = size(A,2)
        if( n/= m) then
            error stop "Expecting a square matrix."
        end if

        ok = .false.
        allocate(LU(n,n))
        allocate(temp(n+1))
        allocate(INDX(n))
        LU = A
        x = b
        !call LU decomposition routine
        call LUDCMP(LU,n,INDX,D,rc)

        !call appropriate solver if previous return code is ok
        if (rc == 0) then
            call LUBKSB(LU,n,INDX,x)
            ok = .true.
        endif    
        
    end function
    
    pure recursive function mat_inv_reduce(M) result(W)
    real(real64), intent(in) :: M(:,:)
    real(real64), allocatable :: W(:,:)
    real(real64), allocatable :: A(:,:), b(:), c(:), d
    real(real64), allocatable :: A_inv(:,:), A_bc(:,:)
    integer :: n
        n = size(M,1)
        allocate(W(n,n))
        if( n>1) then
            A = M(1:n-1, 1:n-1)
            b = M(1:n-1, n)
            c = M(n, 1:n-1)
            d = M(n, n)
            A_bc = A - vec_outer(b, c)/d
            A_inv = mat_inv_reduce(A)
            W(1:n-1, 1:n-1) = mat_inv_reduce(A_bc)
            W(n,n) = 1/(d - dot_product(c, matmul(A_inv,b)))
            W(1:n-1, n) = matmul(A_inv, b)*W(n,n)
            W(n, 1:n-1) = matmul(c, W(1:n-1, 1:n-1))/d
        else
            W(1,1) = 1/M(1,1)
        end if
    end function
    
    pure recursive function mat_solve_vec_reduce(M,z) result(w)
    real(real64), intent(in) :: M(:,:), z(:)
    real(real64), allocatable :: w(:)
    real(real64), allocatable :: A(:,:), b(:), c(:), d, u(:), v
    real(real64), allocatable :: A_bc(:,:)
    integer :: n
        n = size(M,1)
        allocate(w(n))
        if( n>1) then
            A = M(1:n-1, 1:n-1)
            b = M(1:n-1, n)
            c = M(n, 1:n-1)
            d = M(n, n)
            u = z(1:n-1)
            v = z(n)
            A_bc = A - vec_outer(b,c)/d
            w(1:n-1) = mat_solve_vec_reduce(A_bc, u - b*(v/d))
            w(n) = (v-dot_product(c, w(1:n-1)))/d
        else
            w(1) = z(1)/M(1,1)
        end if
    end function
    
    pure function t_mat_solve(obj, b) result(x)
    class(t_matrix(*,*)), intent(in) :: obj
    real(kind=obj%knd), intent(in) :: b(:)
    real(kind=obj%knd), allocatable :: x(:)
        x = lu_mat_solve_vec(obj%data, b)
    end function
    
    pure function t_mat_invert(obj) result(A)
    class(t_matrix(*,*)), intent(in) :: obj
    real(kind=obj%knd), allocatable :: A(:,:)
        A = lu_mat_invert(obj%data)
    end function
    
    subroutine test_array_inv()    
    real(real64), allocatable :: A(:,:), b(:), x(:), A_inv(:,:)
    real(real64) :: d
    integer n
    
    do n=1, 4
        allocate(A(n,n))
        allocate(b(n))
        call RANDOM_NUMBER(A)
        call RANDOM_NUMBER(b)
        b = -1 + 2._wp*b
        
        print *,"Size = ", n
        
        print *, "A="
        call show(A)
                
        d = det(A)
        print *, "det(A) = ", d
        
        A_inv = inv(A)
        print *, "inv(A)="
        call show(A_inv)
        
        print *, "b="
        x = solve(A,b)
        print *, "x="
        call show(x)
        print *, "max(error)=", maxval(abs(b - matmul(A,x)))
        
        deallocate(A)
        deallocate(b)
    end do
    
    end subroutine
    
    
end module