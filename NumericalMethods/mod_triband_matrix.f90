    module mod_triband_matrix
    use mod_common
    
    contains
    
    subroutine triband_decomposition(a_diag, A, D, L, D_inv, L_inv)
    ! Decompose a triband matrix into a lower triangular and a diagonal matrix, as well
    ! as their inverses for quick solution of linear equations
    ! The `A` matrix is specified by the diagonal elements `a_diag` and the upper
    !   and lower band is assumed to be filled with ones.
    ! The `D` matrix contains only diagonal elements.
    ! The `L` matrix has lower band elements equal to `L(i+1,i) = 1/d(i)`
    ! The `D_inv` matrix contains only diagonal elements.
    ! The `L_inv` matrix contains lower triangular elements
    !
    ! calls: triband_decomposition_core
    
    real(real64), intent(in) :: a_diag(:)
    real(real64), intent(out), allocatable :: A(:,:), D(:,:), L(:,:), D_inv(:,:), L_inv(:,:)
    real(real64), allocatable :: d_diag(:), l_lower_inv(:)
    integer :: n, i, j, k
    
    !tex: Performs the decompostion $\mathbf{A}=\mathbf{L}\,\mathbf{D}\,\mathbf{L}^{\top}$ for
    ! a matrix $A$ with the following specifc structure:
    ! $$
    ! \begin{bmatrix}a_{1} & 1\\
    ! 1 & a_{2} & 1\\
     ! & \ddots & \ddots & \ddots\\
     ! &  & 1 & a_{n-1} & 1\\
     ! &  &  & 1 & a_{n}
    ! \end{bmatrix}=\begin{bmatrix}1\\
    ! \tfrac{1}{d_{1}} & 1\\
    ! 0 & \tfrac{1}{d_{2}} & \ddots\\
    ! \ddots & 0 & \ddots & 1\\
    ! 0 & \ddots & 0 & \tfrac{1}{d_{n-1}} & 1
    ! \end{bmatrix}\begin{bmatrix}d_{1}\\
     ! & d_{2}\\
     ! &  & \ddots\\
     ! &  &  & d_{n-1}\\
     ! &  &  &  & d_{n}
    ! \end{bmatrix}\begin{bmatrix}1\\
    ! \tfrac{1}{d_{1}} & 1\\
    ! 0 & \tfrac{1}{d_{2}} & \ddots\\
    ! \ddots & 0 & \ddots & 1\\
    ! 0 & \ddots & 0 & \tfrac{1}{d_{n-1}} & 1
    ! \end{bmatrix}^{\top}!
    !$$
    ! in addtion to the elements of the lower triangular matrix inverse
    ! $$\mathbf{L}^{-1}=\begin{bmatrix}1\\
    !-\tfrac{1}{d_{1}} & 1\\
    !\tfrac{1}{d_{1}d_{2}} & -\tfrac{1}{d_{2}} & \ddots\\
    !\vdots & \ddots & \ddots & 1\\
    !\pm\tfrac{1}{d_{1}\cdots d_{k-2}d_{k-1}} & \cdots & \tfrac{1}{d_{k-2}d_{k-1}} & -\tfrac{1}{d_{k-1}} & 1
    !\end{bmatrix}$$
    
        n = size(a_diag)
        allocate(A(n,n), source = 0d0)
        allocate(D(n,n), source = 0d0)
        allocate(L(n,n), source = 0d0)
        allocate(D_inv(n,n), source = 0d0)
        allocate(L_inv(n,n), source = 0d0)
        
        forall(i=1:n)
            L(i,i) = 1d0
            L_inv(i,i) = 1d0
        end forall
        
        call triband_decomposition_core(a_diag, d_diag, l_lower_inv)
        
        ! Fill A matrix
        do i=1,n
            A(i,i) = a_diag(i)
            if(i>1) then
                A(i-1,i) = 1d0
                A(i,i-1) = 1d0
            end if
        end do
        
        ! Fill D & D_inv matrix
        do i=1, n
            D(i,i) = d_diag(i)
            D_inv(i,i) = 1/d_diag(i)
        end do
        
        ! Fill L & L_inv matrix
        k = 0
        do j=1, n-1
            L(j+1,j) = 1/d_diag(j)
            do i=j+1, n
                k = k + 1
                L_inv(i,j) = l_lower_inv(k)
            end do
        end do        
        
    end subroutine
    
    subroutine triband_decomposition_core(a_diag, d_diag, l_lower_inv)
    ! Decompose a triband matrix into a lower triangular and a diagonal matrix.
    ! The `A` matrix is specified by the diagonal elements `a_diag` and the upper
    !   and lower band is assumed to be filled with ones.
    ! The `D` matrix is specified by the diagonal elements `d_diag`
    ! The `L` matrix has lower band elements equal to `L(i+1,i) = 1/d(i)`
    ! The `inv(L)` matrix is specified by the `n*(n-1)/2` lower triangular elements
    !   arranged in columns (row major).
    real(real64), intent(in) :: a_diag(:)
    real(real64), intent(out), allocatable :: d_diag(:), l_lower_inv(:)
    integer i, j, n, k, r
    real(real64) :: p
    n = size(a_diag)
    
    allocate(d_diag(n))
    do i=1, n
        if (i>1) then
            d_diag(i) = a_diag(i) - 1/d_diag(i-1)
        else
            d_diag(i) = a_diag(i)
        end if
    end do
    
    allocate(l_lower_inv(n*(n-1)/2))
    k = 0
    do j=1, n-1
        do i=j+1,n
            k = k + 1
            p = dble(1 - 2*mod(i-j,2))
            do r=j, i-1
                p = p * d_diag(r)
            end do
            l_lower_inv(k) = 1/p
        end do
    end do
    
    end subroutine
    
    function solve_triband_system(a_diag, b_const) result(x)
    ! Solves the system of equations `A*x=b` where `A` is a triband matrix
    ! with `a_diag` elements on the diagoan, ones on the upper and lower band
    ! and zeros everywhere else.
    !
    ! calls: triband_decomposition
    real(real64), intent(in) :: a_diag(:), b_const(:)
    real(real64), allocatable :: x(:)
    real(real64), allocatable :: A(:,:), D(:,:), D_inv(:,:), L(:,:), L_inv(:,:)
    real(real64), allocatable :: u(:), v(:)
    integer :: n
        
        n = size(a_diag)
        call triband_decomposition(a_diag, A, D, L, D_inv, L_inv)
    
        !tex: If $A = L D L^\top$, to solve $Ax=b$ follow these steps: $\\$
        ! 1. Solve $\mathbf{L}\,\boldsymbol{u}=\boldsymbol{b}$ for $\boldsymbol{u}$. $\\$
        ! 2. Solve $\mathbf{D}\,\boldsymbol{v}=\boldsymbol{u}$ for $\boldsymbol{v}$. $\\$
        ! 3. Solve $\mathbf{L}^{\top}\boldsymbol{x}=\boldsymbol{v}$ for $\boldsymbol{x}$.    
    
        u = matmul(L_inv, b_const)
        v = matmul(D_inv, u)                ! opt. using `d_diag`
        x = matmul(transpose(L_inv), v)     ! opt. using `matmul(v, L_inv)`
        
    end function
    
    
end module