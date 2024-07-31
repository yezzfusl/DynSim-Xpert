! rigid_body_dynamics.f90

module rigid_body_dynamics
  use constants
  use vector_operations
  use numerical_integration
  implicit none

  type :: RigidBody
    real, dimension(3) :: position
    real, dimension(3) :: velocity
    real, dimension(3) :: acceleration
    real, dimension(3) :: angular_velocity
    real, dimension(3) :: angular_acceleration
    real, dimension(3,3) :: inertia_tensor
    real :: mass
  contains
    procedure :: update => update_rigid_body
  end type RigidBody

contains

  subroutine update_rigid_body(this, dt, external_force, external_torque)
    class(RigidBody), intent(inout) :: this
    real, intent(in) :: dt
    real, dimension(3), intent(in) :: external_force, external_torque
    real, dimension(3) :: linear_momentum, angular_momentum
    real, dimension(3,3) :: rotation_matrix

    ! Update linear motion using Euler integration
    this%acceleration = external_force / this%mass
    this%velocity = this%velocity + this%acceleration * dt
    this%position = this%position + this%velocity * dt

    ! Update angular motion
    this%angular_acceleration = matmul(inverse(this%inertia_tensor), &
                                       external_torque - cross_product(this%angular_velocity, &
                                       matmul(this%inertia_tensor, this%angular_velocity)))
    this%angular_velocity = this%angular_velocity + this%angular_acceleration * dt

    ! Update orientation (simplified, using small angle approximation)
    rotation_matrix = eye(3) + skew_symmetric(this%angular_velocity * dt)
    this%inertia_tensor = matmul(matmul(rotation_matrix, this%inertia_tensor), transpose(rotation_matrix))
  end subroutine update_rigid_body

  function inverse(A) result(Ainv)
    real, dimension(3,3), intent(in) :: A
    real, dimension(3,3) :: Ainv
    real :: det

    det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
        - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
        + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

    Ainv(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))/det
    Ainv(1,2) = (A(1,3)*A(3,2)-A(1,2)*A(3,3))/det
    Ainv(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(2,2))/det
    Ainv(2,1) = (A(2,3)*A(3,1)-A(2,1)*A(3,3))/det
    Ainv(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1))/det
    Ainv(2,3) = (A(1,3)*A(2,1)-A(1,1)*A(2,3))/det
    Ainv(3,1) = (A(2,1)*A(3,2)-A(2,2)*A(3,1))/det
    Ainv(3,2) = (A(1,2)*A(3,1)-A(1,1)*A(3,2))/det
    Ainv(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))/det
  end function inverse

  function eye(n) result(A)
    integer, intent(in) :: n
    real, dimension(n,n) :: A
    integer :: i
    A = 0.0
    do i = 1, n
      A(i,i) = 1.0
    end do
  end function eye

  function skew_symmetric(v) result(S)
    real, dimension(3), intent(in) :: v
    real, dimension(3,3) :: S
    S(1,:) = [0.0, -v(3), v(2)]
    S(2,:) = [v(3), 0.0, -v(1)]
    S(3,:) = [-v(2), v(1), 0.0]
  end function skew_symmetric

end module rigid_body_dynamics

! Update main program in dynamics_simulation.f90
program main
  use constants
  use vector_operations
  use numerical_integration
  use rigid_body_dynamics
  implicit none

  type(RigidBody) :: body
  real :: dt = 0.01
  real, dimension(3) :: external_force, external_torque
  integer :: i, num_steps = 1000

  ! Initialize rigid body
  body%mass = 1.0
  body%position = [0.0, 0.0, 0.0]
  body%velocity = [0.0, 0.0, 0.0]
  body%angular_velocity = [0.0, 0.0, 0.0]
  body%inertia_tensor = reshape([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], [3,3])

  external_force = [0.0, 0.0, -body%mass * G]  ! Gravity
  external_torque = [0.1, 0.0, 0.0]  ! Small torque around x-axis

  print *, "Simulating rigid body dynamics..."
  do i = 1, num_steps
    call body%update(dt, external_force, external_torque)
    if (mod(i, 100) == 0) then
      print *, "Time:", i*dt, "Position:", body%position, "Velocity:", body%velocity
    end if
  end do

  print *, "Simulation complete."

end program main
