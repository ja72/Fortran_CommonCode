    module mod_show_matrix
    use mod_common
    implicit none
    
    interface show
        procedure show_vector_i, show_vector_r, show_vector_d
        procedure show_matrix_i, show_matrix_r, show_matrix_d
    end interface
    
    contains
    
        subroutine show_vector_i(v, w)
    ! Display the vector 'v' in a single column
    !   v : the array of real numbers
    !   w : the column width. default = 5
    !   s : sig. figures w-5 (calculated)
        integer, intent(in) :: v(:)
        integer, intent(in), optional :: w
        integer :: i,n,wt
        character(len=16) :: fmt
        if(present(w)) then
            wt = w
        else
            wt = 5
        end if
        n = size(v)
        write( fmt, "(a,g0,a)") "(*(1x,g",wt,".0))"
        write( * , fmt ) ( v(i), new_line("A"), i=1,n )    
    end subroutine
       
    subroutine show_vector_r(v, w)
    ! Display the vector 'v' in a single column
    !   v : the array of real numbers
    !   w : the column width. default = 12
    !   s : sig. figures w-5 (calculated)
        real(real32), intent(in) :: v(:)
        integer, intent(in), optional :: w
        integer :: i,n,dg,wt
        character(len=16) :: fmt
        if(present(w)) then
            wt = w
        else
            wt = 12
        end if
        dg = wt - 6
        n = size(v)
        write( fmt, "(a,g0,a,g0,a)") "(*(1x,g",wt,".",dg,"))"
        write( * , fmt ) ( v(i), new_line("A"), i=1,n )    
    end subroutine
    
    subroutine show_vector_d(v, w)
    ! Display the vector 'v' in a single column
    !   v : the array of real numbers
    !   w : the column width. default = 12
    !   s : sig. figures w-5 (calculated)
        real(real64), intent(in) :: v(:)
        real(real32), allocatable :: u(:)
        integer, intent(in), optional :: w
        u =real(v)
        where( abs(u)<1e-11 )
            u = 0.0
        end where
        call show_vector_r(u,w)
    end subroutine
    
    subroutine show_matrix_i(A, w)
    ! Display the matrix 'A' in columns
    !   A : the array of integers
    !   w : the column width. default = 5
        integer, intent(in) :: A(:,:)
        integer, intent(in), optional :: w
        integer :: i,j,n,m, wt
        character(len=16) :: fmt
        if(present(w)) then
            wt = w
        else
            wt = 5
        end if
        n = size(A,1)
        m = size(A,2)
        write( fmt, "(a,g0,a)") "(*(1x,g",wt,".0))"        
        write( * , fmt ) ( (A(i,j),j=1,m), new_line("A"), i=1,n )
    end subroutine
    
    subroutine show_matrix_r(A, w)
    ! Display the matrix 'A' in columns
    !   A : the array of real numbers
    !   w : the column width. default = 12
    !   s : sig. figures w-5 (calculated)
        real(real32), intent(in) :: A(:,:)
        integer, intent(in), optional :: w
        integer :: i,j,n,m,dg,wt
        character(len=16) :: fmt
        if(present(w)) then
            wt = w
        else
            wt = 12
        end if
        dg = wt - 6
        n = size(A,1)
        m = size(A,2)
        write( fmt, "(a,g0,a,g0,a)") "(*(1x,g",wt,".",dg,"))"
        write( * , fmt ) ( (A(i,j),j=1,m), new_line("A"), i=1,n )
    end subroutine
    
    subroutine show_matrix_d(A,w)
    ! Display the matrix 'A' in columns
    !   A : the array of dble numbers
    !   w : the column width. default = 12
    ! Converts 'A' into single precision and calls `show_matrix_r`
        real(real64), intent(in) :: A(:,:)
        real(real32), allocatable :: B(:,:)
        integer, intent(in), optional :: w
        B = real(A)
        where( abs(B)<1e-11 )
            B = 0.0
        end where
        call show_matrix_r(B,w)
    end subroutine
    
    subroutine test_show_matrix_random()

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
        
    subroutine test_show_matrix_large
    
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

    
    end module