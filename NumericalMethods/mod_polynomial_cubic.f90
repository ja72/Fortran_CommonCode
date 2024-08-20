module mod_polynomial_cubic
use mod_common
implicit none

    integer, parameter :: wp = real64

    contains
    
    function solve_cubic(Q_0,Q_1,Q_2,Q_3) result(T)
    implicit real(real64) (T)
    real(real64), intent(in) :: Q_0, Q_1, Q_2, Q_3
    real(real64) :: T(3)
    real(real64) :: t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t37, t38, t39
    real(real64) :: t36, t40, t41, t42, t43, t44, t45, t46
    
      t25 = 1.0D0/Q_3
      t26 = Q_2**2
      t27 = 1.0D0/Q_3**2
      t32 = Q_0*t25*(1.0D0/2.0D0)
      t33 = 1.0D0/Q_3**3
      t34 = Q_2*t26*t33*(1.0D0/2.7D1)
      t35 = Q_1*Q_2*t27*(1.0D0/6.0D0)
      t28 = t32+t34-t35
      t29 = Q_1*t25*(1.0D0/3.0D0)
      t37 = t26*t27*(1.0D0/9.0D0)
      t30 = t29-t37
      t31 = t30**2
      t36 = t28**2
      t38 = t36+t30*t31
      t39 = sqrt(t38)
      t40 = -t32-t34+t35+t39
      t41 = 1.0D0/t40**(1.0D0/3.0D0)
      t42 = t40**(1.0D0/3.0D0)
      t43 = sqrt(3.0D0)
      t44 = t30*t41
      t45 = t42+t44
      t46 = t30*t41*(1.0D0/2.0D0)
	  
      T(1) = t42-Q_2*t25*(1.0D0/3.0D0)-t30*t41
      T(2) = t42*(-1.0D0/2.0D0)+t46-Q_2*t25*(1.0D0/3.0D+0)-t43*t45*(0.0D0,5.0D-1)
      T(1) = t42*(-1.0D0/2.0D0)+t46-Q_2*t25*(1.0D0/3.0D+0)+t43*t45*(0.0D0,5.0D-1)

    end function
    
end module