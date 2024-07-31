module vector_operations
  implicit none
contains
  function cross_product(a, b)
    real, dimension(3), intent(in) :: a, b
    real, dimension(3) :: cross_product
    
    cross_product(1) = a(2)*b(3) - a(3)*b(2)
    cross_product(2) = a(3)*b(1) - a(1)*b(3)
    cross_product(3) = a(1)*b(2) - a(2)*b(1)
  end function cross_product
  
  function dot_product(a, b)
    real, dimension(:), intent(in) :: a, b
    real :: dot_product
    
    dot_product = sum(a * b)
  end function dot_product

  function norm2(v)
    real, dimension(:), intent(in) :: v
    real :: norm2
    
    norm2 = sqrt(sum(v**2))
  end function norm2
end module vector_operations
