program test_rngs
  use utils, only: fp,ip,dip
  use random
  implicit none
!#define SIMPLE_SPRNG
!#include "sprng_f.h"
  type(sprng_gen) :: sprng_obj
  integer(ip) :: n,ii
  real(fp),dimension(:),allocatable :: x
  !SPRNG_POINTER :: id
  integer(dip) :: id
  real(fp) :: w
  integer :: junk
  
  n = 500
  allocate(x(n))

  ! id = init_sprng(0,0,1,2,SPRNG_DEFAULT)
  ! print *, 'Print information about new stream:'
  ! junk = print_sprng(id)
  ! do ii=1,n
  !    w = sprng(id)
  !    write(*,*) w
  ! end do
  ! do ii=1,n
  !    w = sprng()
  !    write(*,*) w
  ! end do

  call sprng_obj%init(0,1,0)
  call sprng_obj%init(0,2,1)
  call sprng_obj%norm(x,n,0.0_fp,1.0_fp)
  write(*,*) xu
end program test_rngs
