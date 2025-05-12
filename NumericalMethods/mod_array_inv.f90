module mod_array_inv
use mod_common
use mod_show_matrix
implicit none

    !integer, parameter :: wp = real64
    
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
    
    interface   ! mat2
    module pure function mat2_det(A) result(d)
    real(real64) :: d    
    real(real64), intent(in) :: A(2,2)
    end function
    
    module pure function mat2_inv(A) result(B)
    real(real64) :: B(2,2)
    real(real64), intent(in) :: A(2,2)        
    end function
    
    module pure function mat2_solve_vec(A,b) result(x)
    real(real64) :: x(2)
    real(real64), intent(in) :: A(2,2), b(2)
    end function
    module function mat2_eigv(A,V) result(D)
    real(real64) :: D(2)
    real(real64), intent(in) :: A(2,2)
    real(real64), intent(out), optional :: V(2,2)
    end function
    end interface
    
    interface ! mat3
    module pure function mat3_det(A) result(d)
    real(real64) :: d
    real(real64), intent(in) :: A(3,3)
    end function
    
    module pure function mat3_inv(A) result(B)
    real(real64) :: B(3,3)
    real(real64), intent(in) :: A(3,3)            
    end function
    
    module pure function mat3_solve_vec(A,b) result(x)
    real(real64) :: x(3)
    real(real64), intent(in) :: A(3,3), b(3)
    end function
    module function mat3_eigv(A,V) result(D)
    real(real64) :: D(3)
    real(real64), intent(in) :: A(3,3)
    real(real64), intent(out), optional :: V(3,3)
    end function
    end interface
    
    interface ! mat4
    module pure function mat4_det(A) result(d)
    implicit real(real64) (T)
    real(real64) :: d
    real(real64), intent(in) :: A(4,4)
    end function
    
    
    module pure function mat4_inv(A) result(B)
    real(real64) :: B(4,4)
    real(real64), intent(in) :: A(4,4)
    end function

    module pure function mat4_solve_vec(A,b) result(x)
    real(real64) :: x(4)
    real(real64), intent(in) :: A(4,4), b(4)
    end function
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

    function mat_eigv(A,V) result(D)
    real(real64), intent(in) :: A(:,:)
    real(real64), intent(out), optional :: V(size(A,1),size(A,1))
    real(real64) :: d(size(A,1))
    integer :: n, m
    
        n = size(A, 1)
        m = size(A, 2)
        
        if(n /= m) then
            error stop "Expecting square matrix."
        end if
        
        select case(n)
        case(1)  
            d(1) = A(1,1)
            V(1,1) = 1.d0
        case(2)
            d = mat2_eigv(A,V)
        case(3)
            d = mat3_eigv(A,V)
        case default
            error stop 'size not supported in eigv'
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
        b = -1._real64 + 2._real64*b
        
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
    
submodule (mod_array_inv) mod_array_inv_mat2
contains
    module procedure mat2_det
    real(real64) :: t2, t5
        t2 = A(1,1)*A(2,2)
        t5 = A(1,2)*A(2,1)
        d = t2-t5
    end 
    
    module procedure mat2_inv
    real(real64) :: d_inv
        d_inv = 1/mat2_det(A)
        B(1,1) = A(2,2)*d_inv
        B(1,2) = -A(1,2)*d_inv
        B(2,1) = -A(2,1)*d_inv
        B(2,2) = A(1,1)*d_inv        
    end 
    
    module procedure mat2_solve_vec
    real(real64) :: d_inv
        d_inv = 1/mat2_det(A)
        x(1) = d_inv*(A(2,2)*b(1) - A(1,2)*b(2))
        x(2) = d_inv*(A(1,1)*b(2) - A(2,1)*b(1))    
    end 
    module procedure mat2_eigv
    implicit real(real64) (T)
    real(real64) :: t2,t3,t4,t5,t6,t7,t8,t9
    real(real64) :: t10,t11,t12,t13,t14,t15,t16,t17
    real(real64) :: & 
        A_11,A_12, &
        A_21,A_22
    
      A_11 = A(1,1)
      A_21 = A(2,1)
      A_12 = A(1,2)
      A_22 = A(2,2)
      
      t2 = A_11**2
      t3 = A_22**2
      t4 = 1.0D0/A_21
      t5 = A_11*A_22*2.0D0
      t6 = A_12*A_21*4.0D0
      t7 = A_11/2.0D0
      t8 = A_22/2.0D0
      t9 = -t5
      t10 = A_22*t4
      t11 = -t10
      t12 = t2+t3+t6+t9
      t13 = sqrt(t12)
      t14 = t13/2.0D0
      t15 = -t14
      t16 = t7+t8+t14
      t17 = t7+t8+t15
      
      if(present(V)) then
          V(1,1) = t11+t4*t17
          V(2,1) = 1.0D0
          
          V(1,2) = t11+t4*t16
          V(2,2) = 1.0D0
      end if
      
      D(1) = t17
      D(2) = t16
    
    end
end submodule

submodule (mod_array_inv) mod_array_inv_mat3
contains
    module procedure mat3_det
    real(real64) :: t2, t3, t4, t7, t8, t9
        t2 = A(1,1)*A(2,2)*A(3,3)
        t3 = A(1,2)*A(2,3)*A(3,1)
        t4 = A(1,3)*A(2,1)*A(3,2)
        t7 = A(1,1)*A(2,3)*A(3,2)
        t8 = A(1,2)*A(2,1)*A(3,3)
        t9 = A(1,3)*A(2,2)*A(3,1)
        d = t2+t3+t4-t7-t8-t9
    end 
    
    module procedure mat3_inv
    real(real64) :: d_inv
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
    end 
    
    module procedure mat3_solve_vec
    real(real64) :: d_inv
        d_inv = 1/mat3_det(A)
        
        x(1) = d_inv*(A(1,2)*(A(2,3)*b(3)-A(3,3)*b(2))+A(1,3)*(A(3,2)*b(2)-A(2,2)*b(3))+b(1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)))
        x(2) = d_inv*(A(1,1)*(A(3,3)*b(2)-A(2,3)*b(3))+A(1,3)*(A(2,1)*b(3)-A(3,1)*b(2))-b(1)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
        x(3) = d_inv*(A(1,1)*(A(2,2)*b(3)-A(3,2)*b(2))+A(1,2)*(A(3,1)*b(2)-A(2,1)*b(3))+b(1)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    end 
    
    module procedure mat3_eigv
    implicit real(real64) (T)
    real(real64) :: & 
        A_11,A_12,A_13, &
        A_21,A_22,A_23, &
        A_31,A_32,A_33
    
    real(real64) :: t2,t3,t4,t5,t6,t7,t8,t9
    real(real64) :: t10,t11,t12,t13,t14,t15,t16,t17,t18,t19
    real(real64) :: t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
    real(real64) :: t30,t31,t32,t33,t34,t35
    real(real64) :: t40,t41,t42,t43,t44,t45,t46,t47,t48,t49
    real(real64) :: t50,t51,t52,t53,t54,t56,t57,t58
    real(real64) :: t60,t61,t64,t65,t66,t67,t68
    real(real64) :: t73,t74,t75,t76,t78,t79
    real(real64) :: t80,t81,t82,t84,t85,t86,t87,t88,t89
    real(real64) :: t91,t92,t93,t94,t95,t97,t98,t99
    real(real64) :: t100,t101,t102,t103,t104
    
      A_11 = A(1,1)
      A_21 = A(2,1)
      A_31 = A(3,1)
      A_12 = A(1,2)
      A_22 = A(2,2)
      A_32 = A(3,2)
      A_13 = A(1,3)
      A_23 = A(2,3)
      A_33 = A(3,3)
    
      t2 = A_31**2
      t3 = A_32**2
      t4 = A_11*A_22
      t5 = A_12*A_21
      t6 = A_11*A_31
      t7 = A_12*A_31
      t8 = A_11*A_33
      t9 = A_13*A_31
      t10 = A_21*A_32
      t11 = A_22*A_32
      t12 = A_22*A_33
      t13 = A_23*A_32
      t14 = A_31*A_33
      t15 = A_32*A_33
      t24 = A_11+A_22+A_33
      t29 = sqrt(3.0D0)
      t30 = A_11/3.0D0
      t31 = A_22/3.0D0
      t32 = A_33/3.0D0
      t16 = A_32*t6
      t17 = A_33*t6
      t18 = A_33*t7
      t19 = A_32*t9
      t20 = A_31*t11
      t21 = A_33*t10
      t22 = A_31*t13
      t23 = A_33*t11
      t25 = A_12*t2
      t26 = A_13*t2
      t27 = A_21*t3
      t28 = A_23*t3
      t33 = -t5
      t34 = -t9
      t35 = -t13
      t41 = t24**2
      t42 = t24**3
      t43 = t4/3.0D0
      t44 = t5/3.0D0
      t45 = t8/3.0D0
      t46 = t9/3.0D0
      t47 = t12/3.0D0
      t48 = t13/3.0D0
      t49 = (A_33*t4)/2.0D0
      t50 = (A_11*t13)/2.0D0
      t51 = (A_33*t5)/2.0D0
      t52 = (A_23*t7)/2.0D0
      t53 = (A_13*t10)/2.0D0
      t54 = (A_22*t9)/2.0D0
      t65 = t6+t10+t14
      t66 = t7+t11+t15
      t56 = -t43
      t57 = -t45
      t58 = -t47
      t60 = -t50
      t61 = -t51
      t64 = -t54
      t67 = t41/9.0D0
      t68 = t42/2.7D+1
      t73 = t4+t8+t12+t33+t34+t35
      t74 = -1.0D0/(t16-t20-t25+t27)
      t78 = (t18-t19+t23-t28)/(t16-t20-t25+t27)
      t75 = (t24*t73)/6.0D0
      t79 = t44+t46+t48+t56+t57+t58+t67
      t80 = t74*(t17+t21-t22-t26)
      t76 = -t75
      t81 = t79**3
      t82 = -t81
      t84 = (t49+t52+t53+t60+t61+t64+t68+t76)**2
      t85 = t82+t84
      t86 = sqrt(t85)
      t87 = t49+t52+t53+t60+t61+t64+t68+t76+t86
      t88 = t87**(1.0D0/3.0D0)
      t89 = 1.0D0/t88
      t91 = t88/2.0D0
      t92 = -t91
      t93 = t79*t89
      t94 = t93/2.0D0
      t97 = t29*(t88-t93)*(0.0D0,-5.0D-1)
      t98 = t29*(t88-t93)*(0.0D0,5.0D-1)
      t99 = t30+t31+t32+t88+t93
      t95 = -t94
      t100 = t99**2
      t101 = t30+t31+t32+t92+t95+t97
      t102 = t30+t31+t32+t92+t95+t98
      t103 = t101**2
      t104 = t102**2
      
      if(present(V)) then
          V(1,1) = t78+(A_32*t100)/(t16-t20-t25+t27)+t66*t74*t99
          V(2,1) = t80+A_31*t74*t100+(t65*t99)/(t16-t20-t25+t27)
          V(3,1) = 1.0D0
          
          V(1,2) = t78+(A_32*t104)/(t16-t20-t25+t27)+t66*t74*t102
          V(2,2) = t80+A_31*t74*t104+(t65*t102)/(t16-t20-t25+t27)
          V(3,2) = 1.0D0
          
          V(1,3) = t78+(A_32*t103)/(t16-t20-t25+t27)+t66*t74*t101          
          V(2,3) = t80+A_31*t74*t103+(t65*t101)/(t16-t20-t25+t27)          
          V(3,3) = 1.0D0
      end if
            
      D(1) = t99
      D(2) = t102
      D(3) = t101
    end
end submodule

submodule (mod_array_inv) mod_array_inv_mat4
contains
    module procedure mat4_det
    !implicit real(real64) (T)
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
    end 
    
    
    module procedure mat4_inv
    real(real64) :: d_inv
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

    end 

    module procedure mat4_solve_vec
        x = matmul(mat4_inv(A), b)
    end 
end submodule