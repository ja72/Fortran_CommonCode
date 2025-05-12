    module mod_polynomial
    use mod_common
    implicit none

    integer, parameter :: wp = real64

    type :: polynomial
        real(real64), allocatable :: coef(:)
    contains
    procedure :: order => p_order
    procedure :: eval => p_eval_x
    procedure :: show => p_to_string
    end type
    
    interface order
    module procedure :: p_order
    end interface
    
    interface operator (+)
    module procedure :: p_add
    end interface
    interface operator (-)
    module procedure :: p_sub
    end interface
    interface operator (*)
    module procedure :: p_scale1, p_scale2, p_mul
    end interface
    
    interface assignment (=)
    module procedure p_from_array, array_from_p
    end interface

    contains
    
    pure subroutine p_from_array(p,a)
    type(polynomial), intent(out) :: p
    real(real64), intent(in) :: a(:)
    integer :: order
        order = size(a)-1
        allocate(p%coef(0:order))
        p%coef(0:order) = a(:)
    end subroutine
    
    pure subroutine array_from_p(a,p)
    real(real64), intent(out), allocatable :: a(:)
    class(polynomial), intent(in) :: p
    integer :: o
        o = order(p)
        allocate(a(1:o+1))
        a(:) = p%coef(0:o)
    end subroutine
    
    pure function new_const(x) result(p)
    real(real64), intent(in) :: x
    type(polynomial) :: p
        allocate(p%coef(0:0))
        p%coef(0) = x
    end function
    
    pure function p_order(p) result(n)
    integer :: n
    class(polynomial), intent(in) :: p
        n = size(p%coef)-1
    end function

    pure function p_eval_x(p,x) result(y)
    class(polynomial), intent(in) :: p
    real(real64), intent(in) :: x
    real(real64):: y
    integer :: i, order
    order = size(p%coef)-1
    y = 0.0_wp
    do i=0, order
        y = x*y + p%coef(order-i)
    end do
    end function

    pure function p_to_string(p,fmt) result(polystr)
    class(polynomial), intent(in) :: p
    character(len=*), intent(in), optional :: fmt
    !real(real64), intent(in) :: p%coef(0:)
    character(len=:), allocatable :: polystr, buffer
    character(32) :: s
    integer :: i

    if(present(fmt)) then
        buffer = fmt
    else
        buffer = "g0.4"
    end if

    write(s,'("(",' // trim(buffer) // ',")")') p%coef(0)
    polystr = trim(s)

    do i = 1, ubound(p%coef,1)
        !write(s,'("(",g0.4,")")') p%coef(i)
        write(s,'("(",' // trim(buffer) // ',")")') p%coef(i)
        if (i == 1) then
            s = trim(s) // '*x'
        else if (i > 1) then
            s = trim(s) // '*x^' // achar(iachar('0') + i)
        end if
        polystr = polystr // ' + ' // trim(s)
    end do
    end function

    ! Algebra
    
    pure function p_add(a,b) result(c)
    type(polynomial), intent(in) :: a,b
    type(polynomial) :: c
    integer :: na, nb, nc
        na = order(a)
        nb = order(b)
        nc = max(na,nb)
        allocate(c%coef(0:nc))
        c%coef = 0.0_wp
        c%coef(0:na) = a%coef(0:na)
        c%coef(0:nb) = c%coef(0:nb) + b%coef(0:nb)
    end function
    pure function p_sub(a,b) result(c)
    type(polynomial), intent(in) :: a,b
    type(polynomial) :: c
    integer :: na, nb, nc
        na = order(a)
        nb = order(b)
        nc = max(na,nb)
        allocate(c%coef(0:nc))
        c%coef = 0.0_wp
        c%coef(0:na) = a%coef(0:na)        
        c%coef(0:nb) = c%coef(0:nb) - b%coef(0:nb)
    end function
    pure function p_scale1(a,b) result(c)
    real(real64), intent(in) :: a
    type(polynomial), intent(in) :: b
    type(polynomial) :: c
    integer :: nc
        nc = order(b)
        allocate(c%coef(0:nc))
        c%coef(0:nc) = a * b%coef(0:nc)
    end function
    pure function p_scale2(b,a) result(c)
    real(real64), intent(in) :: a
    type(polynomial), intent(in) :: b
    type(polynomial) :: c
    integer :: nc
        nc = order(b)
        allocate(c%coef(0:nc))
        c%coef(0:nc) = a * b%coef(0:nc)
    end function
    
    pure function p_mul(a,b) result(c)
    type(polynomial), intent(in) :: a, b
    type(polynomial) :: c
    integer :: na, nb, nc
    integer :: i,j,k
        na = order(a)
        nb = order(b)
        nc = na + nb
        allocate(c%coef(0:nc))
        c%coef = 0.0_wp
        do i=0, na
            do j=0, nb
                k = i + j
                c%coef(k) = c%coef(k) + a%coef(i)*b%coef(j)
            end do
        end do
    end function

    end module 