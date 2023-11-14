program test
use mod_common
implicit none

    
    call RANDOM_SEED()
    
    call test_mod_show()
    !call test_mod_vectors()
    !call test_nasa_ode()
    !call test_nasa_quat()
    
    contains
    
    subroutine test_mod_show()
    use mod_show_matrix
    
    integer, parameter :: n = 12, m = 6
    
    integer :: iA(n,m), iV(n)
    real :: rA(n,m), rV(n)
    real(real64) :: dA(n,m), dV(n)
    
    call RANDOM_NUMBER(dV)    
    call RANDOM_NUMBER(dA)
    
    dV = -1000.0d0 + 2000.0d0 * dV
    dA = -1000.0d0 + 2000.0d0 * dA
    
    rV = REAL(dV)
    iV = NINT(dV)
    rA = REAL(dA)
    iA = NINT(dA)
    
    print *, "Vectors"
    call show(iV)
    call show(rV)
    call show(dV)
    
    print *, "Matrices"
    call show(iA)
    call show(rA)
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
    
    
end program