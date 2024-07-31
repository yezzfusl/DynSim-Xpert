module visualization
  use fluid_dynamics
  use multi_body_dynamics
  implicit none

contains

  subroutine write_vtk_file(fluid, bodies, filename)
    type(FluidGrid), intent(in) :: fluid
    type(MultiBodySystem), intent(in) :: bodies
    character(len=*), intent(in) :: filename
    integer :: i, j, k, unit_num
    
    open(newunit=unit_num, file=filename, status='replace')
    
    ! Write VTK header
    write(unit_num, '(a)') '# vtk DataFile Version 3.0'
    write(unit_num, '(a)') 'Fluid and body simulation'
    write(unit_num, '(a)') 'ASCII'
    write(unit_num, '(a)') 'DATASET RECTILINEAR_GRID'
    write(unit_num, '(a,3i6)') 'DIMENSIONS', fluid%nx, fluid%ny, fluid%nz
    
    ! Write coordinates
    write(unit_num, '(a,i6)') 'X_COORDINATES', fluid%nx
    write(unit_num, '(f10.4)') (i*fluid%dx, i=0,fluid%nx-1)
    write(unit_num, '(a,i6)') 'Y_COORDINATES', fluid%ny
    write(unit_num, '(f10.4)') (j*fluid%dy, j=0,fluid%ny-1)
    write(unit_num, '(a,i6)') 'Z_COORDINATES', fluid%nz
    write(unit_num, '(f10.4)') (k*fluid%dz, k=0,fluid%nz-1)
    
    ! Write fluid data
    write(unit_num, '(a,i10)') 'POINT_DATA', fluid%nx*fluid%ny*fluid%nz
    
    ! Density
    write(unit_num, '(a)') 'SCALARS density float'
    write(unit_num, '(a)') 'LOOKUP_TABLE default'
    do k = 1, fluid%nz
      do j = 1, fluid%ny
        do i = 1, fluid%nx
          write(unit_num, '(f10.4)') fluid%cells(i,j,k)%density
        end do
      end do
    end do
    
    ! Velocity
    write(unit_num, '(a)') 'VECTORS velocity float'
    do k = 1, fluid%nz
      do j = 1, fluid%ny
        do i = 1, fluid%nx
          write(unit_num, '(3f10.4)') fluid%cells(i,j,k)%velocity
        end do
      end do
    end do
    
    ! Pressure
    write(unit_num, '(a)') 'SCALARS pressure float'
    write(unit_num, '(a)') 'LOOKUP_TABLE default'
    do k = 1, fluid%nz
      do j = 1, fluid%ny
        do i = 1, fluid%nx
          write(unit_num, '(f10.4)') fluid%cells(i,j,k)%pressure
        end do
      end do
    end do
    
    ! Temperature
    write(unit_num, '(a)') 'SCALARS temperature float'
    write(unit_num, '(a)') 'LOOKUP_TABLE default'
    do k = 1, fluid%nz
      do j = 1, fluid%ny
        do i = 1, fluid%nx
          write(unit_num, '(f10.4)') fluid%cells(i,j,k)%temperature
        end do
      end do
    end do
    
    ! Write body data
    write(unit_num, '(a,i6)') 'POINT_DATA', bodies%num_bodies
    write(unit_num, '(a)') 'VECTORS position float'
    do i = 1, bodies%num_bodies
      write(unit_num, '(3f10.4)') bodies%bodies(i)%position
    end do
    
    write(unit_num, '(a)') 'VECTORS velocity float'
    do i = 1, bodies%num_bodies
      write(unit_num, '(3f10.4)') bodies%bodies(i)%velocity
    end do
    
    write(unit_num, '(a)') 'VECTORS angular_velocity float'
    do i = 1, bodies%num_bodies
      write(unit_num, '(3f10.4)') bodies%bodies(i)%angular_velocity
    end do
    
    close(unit_num)
  end subroutine write_vtk_file

end module visualization
