module numerical_integration
  implicit none
contains
  subroutine runge_kutta_4(y, dydt, t, dt, y_next)
    real, dimension(:), intent(in) :: y
    real, dimension(:), intent(out) :: y_next
    real, intent(in) :: t, dt
    interface
      function dydt(t, y)
        real, intent(in) :: t
        real, dimension(:), intent(in) :: y
        real, dimension(size(y)) :: dydt
      end function dydt
    end interface
    
    real, dimension(size(y)) :: k1, k2, k3, k4
    
    k1 = dt * dydt(t, y)
    k2 = dt * dydt(t + 0.5*dt, y + 0.5*k1)
    k3 = dt * dydt(t + 0.5*dt, y + 0.5*k2)
    k4 = dt * dydt(t + dt, y + k3)
    
    y_next = y + (k1 + 2*k2 + 2*k3 + k4) / 6
  end subroutine runge_kutta_4
end module numerical_integration
