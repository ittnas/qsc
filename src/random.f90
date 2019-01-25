module random
  use utils,only: fp,ip,dip
  implicit none

#include "sprng_f.h"

  interface generate_random_seed
     module procedure generate_random_seed_lin,generate_random_seed_par
  end interface generate_random_seed
  
  type,abstract :: rng
     integer(dip) :: min,max
   contains
     procedure(uniform_interface),deferred :: uniform_int
     procedure :: uniform=>uniform_rng
     !procedure(uniform_interface),deferred :: uniform_array
     procedure :: norm => norm_rng
     procedure :: get_parallel_seed => get_parallel_seed_rng
  end type rng
  
  abstract interface
     subroutine uniform_interface(this,x,N)
       use utils,only:ip,fp,dip
       import rng
       class(rng),intent(inout) :: this
       integer(ip),intent(in) :: N
       integer(dip),dimension(N),intent(out) :: x
     end subroutine uniform_interface
  end interface
  !-----------------------------------------------
  type,extends(rng) :: lcg
     integer(ip) :: a,c,seed
     integer(dip) :: m
   contains
     procedure :: uniform_int=>uniform_lcg
     !procedure :: uniform_array=>uniform_lcg,uniform_scalar
     procedure :: init => init_lcg
     procedure :: delete => delete_lcg
     !procedure :: norm
  end type lcg
  !-----------------------------------------------
  type,extends(rng) :: sprng_gen
     integer(ip) :: seed
     integer(ip) :: nstreams,streamnum
     integer(dip) :: sprng_pointer
   contains
     procedure :: uniform_int=>uniform_sprng
     procedure :: init => init_sprng_gen
     procedure :: delete => delete_sprng_gen
  end type sprng_gen
  
contains
  subroutine uniform_rng(this,x,N,lim_in)
    use utils,only:ip,fp,dip
    implicit none
    class(rng),intent(inout) :: this
    integer(ip),intent(in) :: N
    real(fp),dimension(N),intent(out) :: x
    !integer(dip),dimension(N) :: x_t
    integer(dip),allocatable,dimension(:) :: x_t ! Needs to be allocated from the heap
    real(fp),optional,intent(in),dimension(2) :: lim_in
    real(fp),dimension(2) :: lim
    !write(*,*) 'Here in uniform_rng.'
    if(present(lim_in)) then
       lim = lim_in
    else
       lim(1) = 0.0_fp
       lim(2) = 1.0_fp
    endif

    allocate(x_t(N))
    !x_t = this%uniform_int(N)
    call this%uniform_int(x_t,N)
    x = x_t*(1.0_fp/(this%max-this%min)*(lim(2)-lim(1)) + lim(1))
    deallocate(x_t)
  end subroutine uniform_rng
!------------------------------------------------------------------
  subroutine norm_rng(this,x,N,x0_in,sigma_in)
    use utils,only:ip,fp,dip
    use constants,only: pi
    implicit none
    class(rng),intent(inout) :: this
    integer(ip),intent(in) :: N
    real(fp),intent(in),optional :: x0_in,sigma_in
    real(fp),dimension(:),allocatable :: x_t
    real(fp),dimension(N),intent(out) :: x
    !real(fp),dimension(:) :: x(N)
    !real(fp),allocatable,dimension(:) :: x_t ! Eventhough we know this in advance, allocate it from the heap in order to save space
    real(fp) :: x0,sigma,r,theta
    integer(ip) :: ii
    !write(*,*) 'Here in norm_rng!'
    if(present(x0_in)) then
       x0 = x0_in
    else
       x0 = 0.0_fp
    end if
    
    if(present(sigma_in)) then
       sigma = sigma_in
    else
       sigma = 1.0_fp
    end if
    allocate(x_t(2*N))
    call this%uniform(x_t,2*N)
    do ii=1,N
       r=(-2.0_fp*log(x_t(ii*2-1)))**0.5_fp
       theta = 2.0_fp*pi*x_t(ii*2)
       x(ii) = sigma*r*sin(theta)+x0
    end do
    deallocate(x_t)
  end subroutine norm_rng
!------------------------------------------------------------------
  function get_parallel_seed_rng(this,seed_in,thread_num) result(seed_out)
    use utils,only:ip,fp,dip
    use constants,only: pi
    implicit none
    class(rng),intent(in) :: this
    integer(ip),intent(in) :: seed_in,thread_num
    integer(ip) :: seed_out
    seed_out = seed_in + thread_num
  end function get_parallel_seed_rng
!##############################################################
  subroutine init_lcg(this,seed_in,a_in,c_in,m_in)
    class(lcg),intent(inout) :: this
    integer(ip),intent(in),optional :: a_in,c_in,seed_in
    integer(dip),intent(in),optional :: m_in
    integer(ip),parameter :: a_d = 1103515245, c_d = 12345
    integer(dip),parameter :: m_d = 2147483648_dip
    if(present(seed_in)) then
       this%seed = seed_in
    else
       stop(1)
    end if
    if(present(a_in) .and. present(c_in) .and. present(m_in)) then
       this%a = a_in
       this%c = c_in
       this%m = m_in
    else
       this%a = a_d
       this%c = c_d
       this%m = m_d
    endif
    this%min = 0_dip
    this%max = this%m
  end subroutine init_lcg
!-----------------------------------------------------------------
  subroutine delete_lcg(this)
    class(lcg),intent(inout) :: this
  end subroutine delete_lcg
!-----------------------------------------------------------------
  subroutine uniform_lcg(this,x,N)
    use utils,only:ip,fp,dip
    implicit none
    class(lcg),intent(inout) :: this
    integer(ip),intent(in) :: N
    integer(dip),dimension(N),intent(out) :: x
    integer(dip) :: xt
    integer(ip) :: ii
    !write(*,*) 'Here in uniform_lcg!'
    xt = this%seed
    do ii=1,N
       xt = mod(this%a * xt + this%c,this%m)
       x(ii) = xt
    end do
    this%seed = xt
  end subroutine uniform_lcg
!#############################################################  
  function generate_random_seed_lin() result(seed)
    use utils,only:fp,ip
    implicit none
    integer(ip) :: seed
    integer(ip) :: un,istat,dt(8),pid,t(2),s
    integer(dip) :: count,tms
    
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
       if(seed < 0) then
          seed = -seed
       endif
       !write(*,*) 'Found generator. The seed is: ', seed
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(count)
       if (count /= 0) then
          t = transfer(count, t)
       else
          call date_and_time(values=dt)
          tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24 * 60 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
          t = transfer(tms, t)
       end if
       s = ieor(t(1), t(2))
       pid = getpid() + 1099279 ! Add a prime
       seed = ieor(s, pid)
       if(seed<0) then
          seed = -seed
       endif
       !write(*,*) 'No internal generator found. Using date and pid. The seed is: ', seed
    endif
  end function generate_random_seed_lin

  function generate_random_seed_par(thread_num) result(seed)
    use utils,only:fp,ip
    implicit none
    integer(ip),intent(in) :: thread_num
    integer(ip) :: seed
    seed = generate_random_seed_lin()+thread_num;
  end function generate_random_seed_par
!#############################################################
  subroutine init_sprng_gen(this,seed_in,nstreams_in,streamnum_in)
    class(sprng_gen),intent(inout) :: this
    integer(ip),intent(in),optional :: seed_in,nstreams_in,streamnum_in
    
    if(present(seed_in)) then
       this%seed = seed_in
    else
       this%seed = 0
    end if
    
    if(present(nstreams_in)) then
       this%nstreams = nstreams_in
    else
       this%nstreams = 1
    end if

    if(present(streamnum_in)) then
       this%streamnum = streamnum_in
    else
       this%streamnum = 0
    end if
    this%sprng_pointer = init_sprng(0,this%streamnum,this%nstreams,this%seed,SPRNG_DEFAULT)
    this%max = 2_dip**31_ip-1_dip
    this%min = 0_dip
  end subroutine init_sprng_gen
!-----------------------------------------------------------------
  subroutine delete_sprng_gen(this)
    class(sprng_gen),intent(inout) :: this
    integer(ip) :: junk
    junk = free_sprng(this%sprng_pointer)
  end subroutine delete_sprng_gen
!-----------------------------------------------------------------
  subroutine uniform_sprng(this,x,N)
    use utils,only:ip,fp,dip
    implicit none
    class(sprng_gen),intent(inout) :: this
    integer(ip),intent(in) :: N
    integer(dip),dimension(N),intent(out) :: x
    integer(ip) :: ii
    integer(ip) :: temp
    do ii=1,N
       x(ii) = isprng(this%sprng_pointer)
    end do
  end subroutine uniform_sprng
!#############################################################  
end module random
