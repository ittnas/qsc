module ode_solvers
  use utils, only: fp,ip
  abstract interface
     subroutine ode_solver(x,y0,f,y,dx_in,tol_in)
       use utils, only: fp,ip
       use fclass
       implicit none
       real(fp), dimension(:),intent(in) :: x
       complex(fp), dimension(:),intent(in) :: y0
       complex(fp), dimension(:,:),intent(out) :: y
       real(fp),optional :: dx_in,tol_in
       class(fd) :: f
     end subroutine ode_solver
  end interface
contains
  subroutine test_solver(x,y0,f,y,dx)
    use utils, only: fp,ip
    use fclass
    implicit none
    real(fp),dimension(:),intent(in) :: x
    complex(fp), dimension(:),intent(in) :: y0
    complex(fp), dimension(:,:),intent(out) :: y
    class(fd) :: f
    real(fp),optional :: dx
    call f%value_at(x(1),y0,y(1,:))
  end subroutine test_solver

  subroutine euler(x,y0,f,y,dx_in,tol_in)
    use utils, only: fp,ip
    use hamiltonians
    use fclass
    implicit none
    real(fp), dimension(:),intent(in) :: x
    complex(fp), dimension(:),intent(in) :: y0
    complex(fp), dimension(:,:),intent(out) :: y
    class(fd) :: f
    real(fp),optional :: dx_in,tol_in
    real(fp) :: dx,xc,xn,dxc
    integer(ip) :: ii,n,jj,dim
    complex(fp),dimension(:), allocatable :: fn,yn
    logical :: write_next

    n = size(x)
    dim = size(y0)
    if(n < 2) then
       write(*,*) 'At least 2 query points must be defined.'
       stop
    end if
    
    if(size(y,2) .NE. n) then
       write(*,*) 'Size of output array is not equal to query points.'
       stop
    end if

    if(present(dx_in)) then
       dx = dx_in
       if(dx > (x(2) - x(1))/2.0_fp) then
          dx = x(2) - x(1)
       end if
    else
       dx=x(2) - x(1)
    endif
    allocate(fn(dim))
    allocate(yn(dim))

    y(:,1) = y0
    yn = y0
    ii = 2
    write_next = .false.

    !jj = 0

    do while(ii <= n)
       write_next = .false.
       xc = x(ii-1)
       xn = x(ii)
       do while(.true.)
          if(xc + dx + dx/1000_fp >= xn) then
             dxc = xn - xc 
             write_next = .true.
          else
             dxc = dx
          end if
          !write(*,*) dx,dxc,xc,xn
          call f%value_at(xc,yn,fn)
          yn = yn + dxc*fn
          !jj = jj+1
          xc = xc + dxc
          if(write_next) then
             !write(*,*) jj
             y(:,ii) = yn
             exit
          end if
       end do
       ii = ii+1
    end do
    !write(*,*) ii-1,jj
    deallocate(fn,yn)
  end subroutine euler

  subroutine rk4(x,y0,f,y,dx_in,tol_in)
    use utils, only: fp,ip
    use fclass
    implicit none
    real(fp), dimension(:),intent(in) :: x
    complex(fp), dimension(:),intent(in) :: y0
    complex(fp), dimension(:,:),intent(out) :: y
    class(fd) :: f
    real(fp),optional :: dx_in,tol_in
    real(fp) :: dx,xc,xn,dxc
    integer(ip) :: ii,n,jj,dim
    complex(fp),dimension(:), allocatable :: fn,yn,k1,k2,k3,k4
    logical :: write_next

    n = size(x)
    dim = size(y0)
    if(n < 2) then
       write(*,*) 'At least 2 query points must be defined.'
       stop
    end if
    
    if(size(y,2) .NE. n) then
       write(*,*) 'Size of output array is not equal to query points.'
       stop
    end if

    if(present(dx_in)) then
       dx = dx_in
       if(dx > (x(2) - x(1))/2.0_fp) then
          dx = x(2) - x(1)
       end if
    else
       dx=x(2) - x(1)
    endif

    allocate(fn(dim))
    allocate(yn(dim))
    allocate(k1(dim),k2(dim),k3(dim),k4(dim))
    y(:,1) = y0
    yn = y0
    ii = 2

    do while(ii <= n)
       write_next = .false.
       xc = x(ii-1)
       xn = x(ii)
       do while(.true.)
          if(xc + dx + dx/1000_fp >= xn) then
             dxc = xn - xc 
             write_next = .true.
          else
             dxc = dx
          end if
          
          call f%value_at(xc,yn,k1)
          call f%value_at(xc+dxc/2.0_fp,yn+dxc*k1/2.0_fp,k2)
          call f%value_at(xc+dxc/2.0_fp,yn+dxc*k2/2.0_fp,k3)
          call f%value_at(xc+dxc,yn+dxc*k3,k4)
          yn = yn + dxc*(k1+2*k2+2*k3+k4)/6.0_fp
          xc = xc + dxc
          
          if(write_next) then
             !write(*,*) jj
             y(:,ii) = yn
             exit
          end if
       end do
       ii = ii+1
    end do
    deallocate(k1,k2,k3,k4,fn,yn)
  end subroutine rk4

! THIS IS DEPRACATED
  subroutine rk45_interp_end(x,y0,f,y,dx_in,tol_in)
    use utils, only: fp,ip
    use fclass
    use linfuncs
    implicit none
    real(fp), dimension(:),intent(in) :: x
    complex(fp), dimension(:),intent(in) :: y0
    complex(fp), dimension(:,:),intent(out) :: y
    class(fd) :: f
    real(fp),optional :: dx_in,tol_in
    real(fp) :: dx,xc,xn,dxc,s,tol_min,tol_max,eps
    integer(ip) :: ii,n,jj,dim,dim_yvs
    integer(ip),parameter :: growth_factor = 5
    complex(fp),dimension(:), allocatable :: yn,k1,k2,k3,k4,k5,k6,zn
    complex(fp),dimension(:,:), allocatable :: yvs,tmp_yvs
    real(fp),dimension(:),allocatable :: xvs,tmp_xvs
    logical :: write_next
    real(fp),parameter :: a1 = 0.0_fp,a2 = 0.25_fp,a3 = 3.0_fp/8.0_fp,a4 = 12.0_fp/13.0_fp,&
         a5 = 1.0_fp,a6 = 0.5_fp
    real(fp),parameter :: b11 = 0.0_fp,b21 = 0.25_fp,b31 = 3.0_fp/32.0_fp,b32=9.0_fp/32.0_fp,&
         b41=1932.0_fp/2197.0_fp,b42=-7200.0_fp/2197.0_fp,b43=7296.0_fp/2197.0_fp,&
         b51=439.0_fp/216.0_fp,b52 = -8.0_fp,b53=3680.0_fp/513.0_fp,b54 = -845.0_fp/4104.0_fp,&
         b61 = -8.0_fp/27.0_fp,b62 = 2.0_fp,b63=-3544.0_fp/2565.0_fp,b64=1859.0_fp/4104.0_fp,&
         b65=11.0_fp/40.0_fp
    real(fp),parameter :: c1 = 25.0_fp/216.0_fp,c2 = 0.0_fp,c3=1408.0_fp/2565.0_fp,&
         c4=2197.0_fp/4101.0_fp,c5=-0.2_fp
    real(fp),parameter :: d1 = 16.0_fp/135.0_fp,d2=0.0_fp,d3=6656.0_fp/12825.0_fp,&
         d4=28561.0_fp/56430.0_fp,d5=-9.0_fp/50.0_fp,d6=2.0_fp/55.0_fp

    n = size(x)
    dim = size(y0)
    if(n < 2) then
       write(*,*) 'At least 2 query points must be defined.'
       stop
    end if
    
    if(size(y,2) .NE. n) then
       write(*,*) 'Size of output array is not equal to query points.'
       stop
    end if

    if(present(dx_in)) then
       dx = dx_in
    else
       dx = (x(2) - x(1))/101
    endif

    if(present(tol_in)) then
       tol_max = tol_in
    else
       tol_max = 0.001_fp
    end if
    tol_min = tol_max/1.5_fp
    allocate(yn(dim))
    allocate(k1(dim),k2(dim),k3(dim),k4(dim),k5(dim),k6(dim))
    dim_yvs = max(2*size(x),100)
    allocate(yvs(dim,dim_yvs))
    allocate(xvs(dim_yvs))
    y(:,1) = y0
    yvs(:,1) = y0
    yn = y0
    ii=2

    dx = 0.001_fp
    xc = x(1)

    do while(xc + dx/10000.0_fp < x(n))
       !write(*,*) 'Entering the loop!'
       !write(*,*) 'ii=',ii
       if(xc+dx>x(n)) then
          dx = x(n) - xc
       endif
       !write(*,*) 'x=',xc,', dx=',dx
       if(ii>dim_yvs) then !need to reallocate the data arrays
          allocate(tmp_xvs(dim_yvs))
          allocate(tmp_yvs(dim,dim_yvs))
          tmp_xvs = xvs
          tmp_yvs = yvs
          deallocate(yvs,xvs)
          allocate(yvs(dim,dim_yvs*growth_factor),xvs(dim_yvs*growth_factor))
          yvs(:,1:dim_yvs) = tmp_yvs
          xvs(1:dim_yvs) = tmp_xvs
          deallocate(tmp_yvs,tmp_xvs)
          dim_yvs = growth_factor*dim_yvs
       endif
       call f%value_at(xc+a1,yn,k1)
       call f%value_at(xc+a2*dx,yn + dx*(b21*k1),k2)
       call f%value_at(xc+a3*dx,yn + dx*(b31*k1+b32*k2),k3)
       call f%value_at(xc+a4*dx,yn + dx*(b41*k1+b42*k2+b43*k3),k4)
       call f%value_at(xc+a5*dx,yn + dx*(b51*k1+b52*k2+b53*k3+b54*k4),k5)
       call f%value_at(xc+a6*dx,yn + dx*(b61*k1+b62*k2+b64*k3+b64*k4+b65*k5),k6)
       yn = yn + dx*(c1*k1+c3*k3+c4*k4+c5*k5)
       zn = yn + dx*(d1*k1+d3*k3+d4*k4+d5*k5+d6*k6)
       eps = norm2(abs(zn-yn))
       !write(*,*) 'eps=',eps,', tol_max=',tol_max
       if(eps > tol_max) then ! reject!
          s = 0.84_fp*(tol_min/eps)**0.25_fp ! and reduce step size
          dx = s*dx
       else ! accept
          xc = xc + dx
          xvs(ii) = xc
          yvs(:,ii) = yn
          !yn = zn
          ii = ii+1
          if(eps < tol_min) then ! but increase step
             s = 0.84_fp*(tol_min/eps)**0.25_fp
             dx = s*dx
          endif
       endif
    end do
    y = interp1(xvs(1:(ii-1)),yvs(:,1:(ii-1)),x)
    deallocate(k1,k2,k3,k4,k5,k6,yn,zn,yvs,xvs)
    write(*,*) 'dx:',dx,', ii',ii
  end subroutine rk45_interp_end
!------------------------------------------
  subroutine rk45(x,y0,f,y,dx_in,tol_in)
    use utils, only: fp,ip
    use fclass
    use linfuncs
    implicit none
    real(fp), dimension(:),intent(in) :: x
    complex(fp), dimension(:),intent(in) :: y0
    complex(fp), dimension(:,:),intent(out) :: y
    class(fd) :: f
    real(fp),optional :: dx_in,tol_in
    real(fp) :: dx,xc,xn,dxc,s,tol_min,tol_max,eps,dx_min,dx_max
    integer(ip) :: ii,n,jj,dim,dim_yvs,jj_t,kk
    integer(ip),parameter :: growth_factor = 5
    complex(fp),dimension(:), allocatable :: yn,k1,k2,k3,k4,k5,k6,zn,yn_t
    complex(fp),dimension(:,:), allocatable :: yvs,yvs_t
    real(fp),dimension(:),allocatable :: xvs,xvs_t
    real(fp),parameter :: a1 = 0.0_fp,a2 = 0.25_fp,a3 = 3.0_fp/8.0_fp,a4 = 12.0_fp/13.0_fp,&
         a5 = 1.0_fp,a6 = 0.5_fp
    real(fp),parameter :: b11 = 0.0_fp,b21 = 0.25_fp,b31 = 3.0_fp/32.0_fp,b32=9.0_fp/32.0_fp,&
         b41=1932.0_fp/2197.0_fp,b42=-7200.0_fp/2197.0_fp,b43=7296.0_fp/2197.0_fp,&
         b51=439.0_fp/216.0_fp,b52 = -8.0_fp,b53=3680.0_fp/513.0_fp,b54 = -845.0_fp/4104.0_fp,&
         b61 = -8.0_fp/27.0_fp,b62 = 2.0_fp,b63=-3544.0_fp/2565.0_fp,b64=1859.0_fp/4104.0_fp,&
         b65=11.0_fp/40.0_fp
    real(fp),parameter :: c1 = 25.0_fp/216.0_fp,c2 = 0.0_fp,c3=1408.0_fp/2565.0_fp,&
         c4=2197.0_fp/4101.0_fp,c5=-0.2_fp
    real(fp),parameter :: d1 = 16.0_fp/135.0_fp,d2=0.0_fp,d3=6656.0_fp/12825.0_fp,&
         d4=28561.0_fp/56430.0_fp,d5=-9.0_fp/50.0_fp,d6=2.0_fp/55.0_fp

    n = size(x)
    dim = size(y0)
    if(n < 2) then
       write(*,*) 'At least 2 query points must be defined.'
       stop
    end if
    
    if(size(y,2) .NE. n) then
       write(*,*) 'Size of output array is not equal to query points.'
       stop
    end if

    if(present(dx_in)) then
       dx = dx_in
    else
       dx = (x(2) - x(1))/101.0_fp
    endif

    dx_max = dx

    if(present(tol_in)) then
       tol_max = tol_in
    else
       tol_max = 0.001_fp
    end if
    tol_min = tol_max/1.5_fp
    allocate(yn(dim),yn_t(dim),zn(dim))
    allocate(k1(dim),k2(dim),k3(dim),k4(dim),k5(dim),k6(dim))
    allocate(yvs(dim,4),yvs_t(dim,4))
    allocate(xvs(4),xvs_t(4))
    xvs = 0.0_fp
    xvs_t = 0.0_fp
    y(:,1) = y0
    yvs(:,1) = y0
    xvs(1) = x(1)
    yn = y0
    ii=2
    jj=2
    kk=2

    !dx = 0.001_fp
    xc = x(1)

    !do while(xc + dx/10000.0_fp < x(n))
    do while(xc + tiny(1.0_fp) < x(n))
       !write(*,*) 'Entering the loop!'
       !write(*,*) 'ii=',ii
       if(xc+dx>x(n)) then
          dx = x(n) - xc
       endif
       !write(*,*) 'x=',xc,', dx=',dx
       call f%value_at(xc+a1,yn,k1)
       call f%value_at(xc+a2*dx,yn + dx*(b21*k1),k2)
       call f%value_at(xc+a3*dx,yn + dx*(b31*k1+b32*k2),k3)
       call f%value_at(xc+a4*dx,yn + dx*(b41*k1+b42*k2+b43*k3),k4)
       call f%value_at(xc+a5*dx,yn + dx*(b51*k1+b52*k2+b53*k3+b54*k4),k5)
       call f%value_at(xc+a6*dx,yn + dx*(b61*k1+b62*k2+b64*k3+b64*k4+b65*k5),k6)
       yn_t = yn + dx*(c1*k1+c3*k3+c4*k4+c5*k5)
       zn = yn + dx*(d1*k1+d3*k3+d4*k4+d5*k5+d6*k6)
       eps = norm2(abs(zn-yn_t))
       !write(*,*) 'eps=',eps,', tol_max=',tol_max
       if(eps > tol_max) then ! reject!
          s = 0.84_fp*(tol_min/eps)**0.25_fp ! and reduce step size
          dx = s*dx
          !write(*,*) 'Too small step at ',xc,'. New dx is :',dx
       else ! accept
          yn = yn_t
          
          xc = xc + dx
          !write(*,*) 'Next x is :',xc
          xvs(ii) = xc
          yvs(:,ii) = yn

          jj_t = jj
          ! Loops as long as the current element in xvs is larger than the highest
          ! unprocessed element in x
          
          do while(jj < size(x) .and. x(jj) <= xvs(ii))
             jj = jj+1
          end do
          ! In order to also get the last element
          if(jj .eq. size(x) .and. x(jj) <= xvs(ii)) then
             jj = jj+1 ! Add one, because below jj-1 is used.
          end if

          if(kk>4) then ! xvs is a circular buffer. Only order if xvs is filled.
             xvs_t = cshift(xvs,ii)
             yvs_t = cshift(yvs,ii,2)
          else
             xvs_t = xvs
             yvs_t = yvs
          endif
          ! jj tells which elements of x we need to process
          if(jj .NE. jj_t) then
             !write(*,*) 'dx:',dx,', ii:',ii,', jj:',jj
             !write(*,*) 'xvs:',xvs
             !write(*,*) 'xvs_t:',xvs_t(1:min(4,kk))
             !write(*,*) 'x(jj):',x(jj)
             !write(*,*) 'x:',x(jj_t:jj-1)
             y(:,jj_t:jj-1) = interp1(xvs_t(1:min(4,kk)),yvs_t(:,1:min(4,kk)),x(jj_t:jj-1))
          end if
          
          !yn = zn
          ii = mod(ii,4)+1
          kk = kk+1
          if(eps < tol_min) then ! but increase step
             s = 0.84_fp*(tol_min/eps)**0.25_fp
             dx = s*dx
             dx = min(dx,(x(n)-x(1))/10.0_fp,dx_max) ! Limit the maximum step size. Otherwise risk infinity
             !write(*,*) 'Too big step at ',xc,'. New dx is :',dx
          endif
       endif
    end do
    !y = interp1(xvs(1:(ii-1)),yvs(:,1:(ii-1)),x)
    deallocate(k1,k2,k3,k4,k5,k6,yn,zn,yvs,xvs)
    !write(*,*) 'RK45 finished. dx:',dx,', nbr of iterations',kk
  end subroutine rk45
end module ode_solvers

! program test_ode_solvers
!   use utils, only: fp,ip
!   use fclass
!   use ode_solvers
!   use linfuncs
!   use hamiltonians
!   use operators
!   use signals
!   use envelopes
!   use random
!   implicit none
!   real(fp),allocatable,dimension(:) :: x,wq,xquery
!   complex(fp),allocatable,dimension(:) :: y0s,y0m
!   complex(fp),allocatable,dimension(:,:) :: ys,ym,yquery
!   complex(fp) :: H(2,2)
!   type(fdstd) :: ftd
!   type(fdsti) :: fti
!   type(fdmtd) :: fmtd

!   type(fdstd) :: ftd_driven
!   integer(ip) :: ii,NN,c1,c2,nql,ncl,jj
!   type(simple_hamiltonian) :: Hs
!   type(jch) :: Hjc
!   type(djch) :: Hdjc
!   integer(ip) :: test_dimensions(3)
!   type(operator_cont) :: op
!   procedure(ode_solver),pointer :: solver => null()
!   type(signal) :: test_signal
!   type(constant_envelope) :: cenv
!   type(lcg) :: lcg_t
!   real(fp),dimension(:),allocatable :: rnumber
!   integer(dip),dimension(:),allocatable :: inumber
!   integer(ip) :: nr
!   type(fdlmtd) :: flmtd
  
!   test_dimensions = [2,1,1]
!   NN = 1001
!   nql = 5
!   ncl = 1

!   solver => rk45
!   !call op%initialize(test_dimensions)

!   allocate(y0s(nql*ncl))
!   allocate(y0m((nql*ncl)**2))
!   y0s = (0.0_fp,0.0_fp)
!   y0s(2) = 1.0_fp
!   y0m = reshape(matmul(reshape(y0s,[size(y0s),1]),reshape(y0s,[1,size(y0s)])),[size(y0m)])
!   !write(*,*) real(y0)
!   !y0 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
!   !allocate(y0(NH))
!   allocate(x(NN))
!   allocate(ys(size(y0s),NN))
!   allocate(ym(size(y0m),NN))
!   allocate(wq(nql))
!   allocate(xquery(NN/3),yquery(size(y0s),NN/3))
!   wq = 0.0_fp
!   do ii=1,size(wq)
!      wq(ii) = ii**0.7
!   end do

!   call linspace(0.0_fp,1000.0_fp,NN,x)
!   call linspace(0.0_fp,0.1_fp,NN/3,xquery)
!   !  f = fd(H)

!   !call Hs%initialize_simple_H(2,0.0_fp,1.0_fp)
!   call Hdjc%initialize_djch(nql,ncl,wq,0.5_fp,1.0_fp,hbar_in = 1.0_fp)
!   call test_signal%initialize(Hdjc%op%cc(:,:,1),1.0_fp)
!   call cenv%initialize_constant_envelope((0.01_fp,0.0_fp))
!   call test_signal%add_envelope(cenv)
!   call Hdjc%add_signal(test_signal)

!   call Hjc%initialize_jch(nql,ncl,wq,0.5_fp,1.0_fp,hbar_in = 1.0_fp)
!   call fti%fdsti_initialize(Hjc%value_at(0.0_fp))
!   call ftd%fdstd_initialize(Hjc)
!   call ftd_driven%fdstd_initialize(Hdjc)
!   call fmtd%fdmtd_initialize(Hjc)

!   call flmtd%initialize_fdlmtd(Hdjc)
!   !call flmtd%add_lindblad_operator_fdlmtd(Hjc%op%aa(:,:,1),(0.01_fp,0.0_fp))
!   !call flmtd%add_lindblad_operator_fdlmtd(cmplx(get_diag_matrix(wq),kind=fp),(0.5_fp,0.0_fp))
!   call system_clock(c1)
!   !call euler(x,y0m,fmtd,ym,0.00001_fp)
!   !call euler(x,y0m,fmtd,ym,0.001_fp)
!   call solver(x,y0m,flmtd,ym,0.01_fp,0.01_fp)
!   call system_clock(c2)
!   write(*,*) 'DM execution took', (c2-c1), 'milliseconds'
!   call system_clock(c1)
!   !call euler(x,y0s,ftd,ys,0.00001_fp)
!   !call euler(x,y0s,ftd,ys,0.001_fp)
!   call rk4(x,y0s,ftd,ys,0.001_fp,0.01_fp)
!   call system_clock(c2)
!   write(*,*) 'SE execution with Hjc took', (c2-c1), 'milliseconds'
!   call system_clock(c1)
!   call rk4(x,y0s,ftd_driven,ys,0.001_fp,0.01_fp)
!   call system_clock(c2)
!   write(*,*) 'SE execution with Hdjc took', (c2-c1), 'milliseconds'

  
!   !call print_matrix(Hjc%value_at(0.0_fp),.true.)
!   !call print_matrix(Hjc%qcca,.true.)
!   !call print_matrix(Hjc%qacc,.true.)
!   !write(*,*) real(CONJG(y(1,:))*y(1,:))
  
!   yquery = interp1(x,ys,xquery)
  
!   open(unit = 1, file = 'test_data_s.dat')
!   open(unit = 2, file = 'test_data_m.dat')
!   open(unit = 3, file = 'test_data_q.dat')
!   do ii=1,size(x)
!      !write(1,*) x(ii),real(CONJG(y(:,ii))*y(:,ii))
!      write(1,*) x(ii),real(conjg(ys(:,ii))*ys(:,ii))
!      !write(2,*) x(ii),real(ym(:,ii))
!      write(2,"(f8.3)",advance='no') x(ii)
!      do jj=1,nql*ncl
!         write(2,"(f8.3)",advance='no') ,real(ym(nql*ncl*(jj-1) + jj,ii))
!      end do
!      write(2,*) ''
!   end do
!   do ii=1,size(xquery)
!      write(3,*) xquery(ii),real(conjg(yquery(:,ii))*yquery(:,ii))
!   end do
!   close(1)
!   close(2)
!   close(3)
  
!   !write(*,*) Hdjc%value_at(1.0_fp)
!   call lcg_t%init(generate_random_seed())
!   !rnumber = lcg_t%uniform(1)
!   !write(*,*) rnumber
!   !nr = 20
!   !allocate(inumber(nr),rnumber(nr))
!   !inumber = lcg_t%uniform_int(nr)
!   !rnumber = lcg_t%norm(nr)
!   !write(*,*) rnumber

! end program test_ode_solvers
 
