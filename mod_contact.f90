    module mod_contact
    use mod_common
    implicit none
    
    type t_contact        
        real(real64) :: dx, dy
        integer :: nx, ny
        real(real64) :: k       ! k = (1-ν_1)/(π*E_1)+(1-ν_2)/(π*E_2)
        real(real64), allocatable :: P(:), delta(:)
        real(real64) :: CP(4,4)    ! pressure coeffients inverted matrix
    contains
        procedure :: ctx_materials
        procedure :: ctx_patch_pressures
        procedure :: ctx_patch_deflection
    end type t_contact
    
    interface t_contact
        module procedure :: ctx_init
    end interface

    contains
    
    pure function ctx_init(lx,ly,nx,ny) result(obj)
    real(real64),intent(in) :: lx,ly
    integer, intent(in) :: nx,ny
    type(t_contact) :: obj
    real(real64) :: area, dx, dy, inv_area
    
        obj%dx = lx
        obj%dy = ly
        obj%nx = nx
        obj%ny = ny
        dx = lx/nx
        dy = ly/ny
        area = dx*dy
        inv_area = 1/area
        obj%CP(:,1) = [         1.0d0, -dx*inv_area     , -dy*inv_area      ,  inv_area     ]
        obj%CP(:,2) = [         0.0d0,  0.0d0           ,  dx*inv_area      , -inv_area     ]
        obj%CP(:,3) = [         0.0d0,  dy*inv_area     ,  0.0d0            , -inv_area     ]
        obj%CP(:,4) = [         0.0d0,  0.0d0           ,  0.0d0            ,  inv_area     ]
        obj%k = 1.4d-6  ! default is steel on steel
        allocate(obj%P(nx*ny))
        allocate(obj%delta(nx*ny))
                    
    end function
    
    pure subroutine ctx_materials(obj, E_1, v_1, E_2, v_2)
    class(t_contact), intent(inout) :: obj
    real(real64), intent(in) :: E_1, v_1
    real(real64), intent(in), optional :: E_2, v_2
        if(present(E_2) .and. present(v_2)) then
            obj%k = (1-v_1**2)/(pi*E_1)+(1-v_2**2)/(pi*E_2)
        else
            obj%k = (1-v_1**2)/(pi*E_1)
        end if
    end subroutine
    
    pure function ctx_patch_pressures(obj, i, j) result(P)
    class(t_contact), intent(in) :: obj
    integer, intent(in) :: i,j
    real(real64) :: P(4)
    integer :: k_11, k_12, k_21, k_22
    
        k_11 = ((i-1)*obj%ny+j)
        k_12 = ((i-1)*obj%ny+j+1)
        k_21 = ( (i) *obj%ny+j)
        k_22 = ( (i) *obj%ny+j+1)
        
        P = obj%P( [k_11, k_12, k_21, k_22] )
    
    end function
    
    pure function ctx_patch_deflection(obj, i, j) result(delta)
    class(t_contact), intent(in) :: obj
    integer, intent(in) :: i,j
    real(real64) :: delta, P(4), C(4), U(4)
    integer :: k_11, k_12, k_21, k_22
    real(real64) :: X_1, X_2, Y_1, Y_2, L_11, L_12, L_21, L_22
    
        k_11 = ((i-1)*obj%ny+j)
        k_12 = ((i-1)*obj%ny+j+1)
        k_21 = ( (i) *obj%ny+j)
        k_22 = ( (i) *obj%ny+j+1)
        
        P = obj%P( [k_11, k_12, k_21, k_22] )
        C = matmul( obj%CP, P)
        
        X_1 = (i-1)*obj%dx
        X_2 = X_1 + obj%dx
        Y_1 = (j-1)*obj%dy
        Y_2 = Y_1 + obj%dy
        
        L_11 = sqrt( X_1**2 + Y_1**2 )
        L_12 = sqrt( X_2**2 + Y_1**2 )
        L_21 = sqrt( X_1**2 + Y_2**2 )
        L_22 = sqrt( X_2**2 + Y_2**2 )
        
        U(1) = X_2*log((L_22+Y_2)/(L_21+Y_1))-X_1*log((L_12+Y_2)/(L_11+Y_1))+Y_2*log((L_22+X_2)/(L_12+X_2))-Y_1*log((L_21+X_2)/(L_11+X_2))
        U(2) = X_2**2/2*log((L_22+Y_2)/(L_21+Y_1))-X_1**2/2*log((L_12+Y_2)/(L_11+Y_1))+Y_2/2*(L_22-L_12)-Y_1/2*(L_21-L_11)
        U(3) = Y_2**2/2*log((L_22+X_2)/(L_12+X_2))-Y_1**2/2*log((L_21+X_2)/(L_11+X_1))+X_2/2*(L_22-L_21)-X_1/2*(L_12-L_11)
        U(4) = (L_22**3 - L_21**3 - L_12**2 + L_11**3)/3
        
        delta = obj%k * dot_product(U, C)
        
    end function
        
    end module