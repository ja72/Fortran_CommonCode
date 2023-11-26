program test
use mod_common
implicit none

    
    call RANDOM_SEED()
    
    !call test_mod_show()
    !call test_mod_vectors()
    !call test_array_inv()
    !call test_nasa_ode()
    !call test_nasa_quat()
    
    call test_rb()
    
    contains
    
    subroutine test_array_inv()
    use mod_array_inv
    use mod_show_matrix
    
    real(wp), allocatable :: A(:,:), b(:), x(:), A_inv(:,:)
    real(wp) :: d
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
    
    subroutine test_mod_show()
    use mod_show_matrix
    
    integer, parameter :: n = 12, m = 6
    
    integer :: iA(n,m), iV(n)
    real(real64) :: dA(n,m), dV(n)
    
    call RANDOM_NUMBER(dV)    
    call RANDOM_NUMBER(dA)
    
    dV = -1000.0d0 + 2000.0d0 * dV
    dA = 1 + 1.0d8 * dA
    
    iV = NINT(dV)
    iA = NINT(dA)
    
    print *, "Vectors"
    call show(iV)
    call show(dV)
    
    print *, "Matrices"
    call show(iA)
    call show(dA)
    
    end subroutine
    
    subroutine test_mod_vectors()
    use mod_spatial_vectors
    use mod_show_matrix
    integer :: n
    type(matrix) :: A
    type(vector) :: c, x, e
    n = 12
    
    A = 5._wp * ident(n) + rand(n,n)
    print *, "A="
    call show(A%data)
    c = rand(n, -1._wp, 1._wp)
    print *, "c="
    call show(c%data)
    x = solve(A,c)
    print *, "x="
    call show(x%data)
    e = c - A*x
    print *, "e="
    call show(e%data)
    end subroutine
    
    subroutine test_show_matrix
    use mod_show_matrix
    integer :: row(16), matrix(4,4)    
    real(real64) :: A(4,4)
    
    row = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]    
    matrix = reshape( row, [4, 4])
           
    call show(matrix)
    
    A = dble(matrix)
    
    A = sqrt( matmul( transpose(A), A) )
    
    call show(A, 8)
    call show(A, 12)
    call show(A, 16)
    
    end subroutine
    
    subroutine test_nasa_ode()
    use mod_nasa_ode    
    call srk4_test()    
    end subroutine
    
    subroutine test_nasa_quat()
    use mod_nasa_quat
    call quat_array_test()
    end subroutine
    
    subroutine test_rb()
    use mod_rigid_bodies
    call test_rb_sim()
    end subroutine
    
    
end program