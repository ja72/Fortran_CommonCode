program test
use mod_common
implicit none

    
    call RANDOM_SEED()
    
    call test_mod_show()
    call test_mod_spatial_vectors()
    call test_mod_array_inv()
    call test_mod_nasa()
    call test_mod_rigid_bodies()
    call test_mod_linalg()
    
    contains
    
    subroutine test_mod_array_inv()
    use mod_array_inv
        call test_array_inv()    
    end subroutine
    
    subroutine test_mod_spatial_vectors()
    use mod_spatial_vectors
        call test_vectors()    
    end subroutine
        
    subroutine test_mod_show()
    use mod_show_matrix
        call test_show_matrix_random()
        call test_show_matrix_large()    
    end subroutine
    
    subroutine test_mod_nasa()
    use mod_nasa_quat
    use mod_nasa_ode    
    call quat_array_test()
    call srk4_test()    
    end subroutine    
    
    subroutine test_mod_rigid_bodies()
    use mod_rigid_bodies
    call test_rb_sim()
    end subroutine
    
    subroutine test_mod_linalg()
    use mod_linear_algebra
        call test_linear_algebra()        
    end subroutine
    
    
end program