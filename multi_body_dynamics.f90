! multi_body_dynamics.f90

module multi_body_dynamics
  use constants
  use vector_operations
  use numerical_integration
  use rigid_body_dynamics
  implicit none

  type :: MultiBodySystem
    type(RigidBody), dimension(:), allocatable :: bodies
    integer :: num_bodies
  contains
    procedure :: initialize => initialize_multi_body_system
    procedure :: update => update_multi_body_system
    procedure :: detect_collisions => detect_collisions_multi_body_system
    procedure :: resolve_collision => resolve_collision_multi_body_system
  end type MultiBodySystem

contains

  subroutine initialize_multi_body_system(this, num_bodies)
    class(MultiBodySystem), intent(inout) :: this
    integer, intent(in) :: num_bodies
    integer :: i

    this%num_bodies = num_bodies
    allocate(this%bodies(num_bodies))

    do i = 1, num_bodies
      this%bodies(i)%mass = 1.0
      this%bodies(i)%position = [real(i), 0.0, 0.0]
      this%bodies(i)%velocity = [0.0, 0.0, 0.0]
      this%bodies(i)%angular_velocity = [0.0, 0.0, 0.0]
      this%bodies(i)%inertia_tensor = reshape([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], [3,3])
    end do
  end subroutine initialize_multi_body_system

  subroutine update_multi_body_system(this, dt)
    class(MultiBodySystem), intent(inout) :: this
    real, intent(in) :: dt
    integer :: i
    real, dimension(3) :: external_force, external_torque
    logical :: collision_occurred
    
    do i = 1, this%num_bodies
      external_force = [0.0, 0.0, -this%bodies(i)%mass * G]  ! Gravity
      external_torque = [0.0, 0.0, 0.0]  ! No external torque
      
      call this%bodies(i)%update(dt, external_force, external_torque)
    end do

    ! Detect and resolve collisions
    call this%detect_collisions(collision_occurred)
    if (collision_occurred) then
      call this%resolve_collision()
    end if
  end subroutine update_multi_body_system

  subroutine detect_collisions_multi_body_system(this, collision_occurred)
    class(MultiBodySystem), intent(inout) :: this
    logical, intent(out) :: collision_occurred
    integer :: i, j
    real :: distance, collision_threshold

    collision_occurred = .false.
    collision_threshold = 1.0  ! Simplified collision detection using a distance threshold

    do i = 1, this%num_bodies - 1
      do j = i + 1, this%num_bodies
        distance = norm2(this%bodies(i)%position - this%bodies(j)%position)
        if (distance < collision_threshold) then
          collision_occurred = .true.
          return
        end if
      end do
    end do
  end subroutine detect_collisions_multi_body_system

  subroutine resolve_collision_multi_body_system(this)
    class(MultiBodySystem), intent(inout) :: this
    integer :: i, j
    real :: distance, collision_threshold
    real, dimension(3) :: normal, relative_velocity
    real :: impulse_magnitude
    real :: coefficient_of_restitution = 0.8

    collision_threshold = 1.0

    do i = 1, this%num_bodies - 1
      do j = i + 1, this%num_bodies
        distance = norm2(this%bodies(i)%position - this%bodies(j)%position)
        if (distance < collision_threshold) then
          ! Compute collision normal
          normal = (this%bodies(j)%position - this%bodies(i)%position) / distance

          ! Compute relative velocity
          relative_velocity = this%bodies(j)%velocity - this%bodies(i)%velocity

          ! Compute impulse magnitude (simplified - assuming direct central collision)
          impulse_magnitude = -(1 + coefficient_of_restitution) * dot_product(relative_velocity, normal) / &
                               (1/this%bodies(i)%mass + 1/this%bodies(j)%mass)

          ! Apply impulse to both bodies
          this%bodies(i)%velocity = this%bodies(i)%velocity - impulse_magnitude * normal / this%bodies(i)%mass
          this%bodies(j)%velocity = this%bodies(j)%velocity + impulse_magnitude * normal / this%bodies(j)%mass

          ! Separate overlapping bodies (simplified)
          this%bodies(i)%position = this%bodies(i)%position - 0.5 * (collision_threshold - distance) * normal
          this%bodies(j)%position = this%bodies(j)%position + 0.5 * (collision_threshold - distance) * normal
        end if
      end do
    end do
  end subroutine resolve_collision_multi_body_system

end module multi_body_dynamics

! Update main program in dynamics_simulation.f90
program main
  use constants
  use vector_operations
  use numerical_integration
  use rigid_body_dynamics
  use fluid_dynamics
  use multi_body_dynamics
  implicit none

  type(MultiBodySystem) :: system
  type(FluidGrid) :: fluid
  real :: dt = 0.01
  integer :: i, num_steps = 1000

  ! Initialize multi-body system
  call system%initialize(3)  ! Create a system with 3 bodies

  ! Initialize fluid grid
  call fluid%initialize(50, 50, 50, 0.1, 0.1, 0.1)

  print *, "Simulating multi-body dynamics with fluid interaction..."
  do i = 1, num_steps
    call system%update(dt)
    call fluid%update(dt)
    
    if (mod(i, 100) == 0) then
      print *, "Time:", i*dt
      print *, "Body 1 - Position:", system%bodies(1)%position, "Velocity:", system%bodies(1)%velocity
      print *, "Body 2 - Position:", system%bodies(2)%position, "Velocity:", system%bodies(2)%velocity
      print *, "Body 3 - Position:", system%bodies(3)%position, "Velocity:", system%bodies(3)%velocity
      print *, "Fluid - Center Cell Density:", fluid%cells(25,25,25)%density, &
               "Velocity:", fluid%cells(25,25,25)%velocity
    end if
  end do

  print *, "Simulation complete."

end program main
