program test_operators
  use operators
  use linfuncs
  implicit none
  type(operator_cont) :: op
  
  call op%initialize([2,2,2])
  call print_matrix(real(op%cc(:,:,2)))
end program test_operators
