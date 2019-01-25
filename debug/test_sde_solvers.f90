 program test_sde_solvers
   use utils,only: fp,ip
   use sclass
   use hamiltonians
   use random
   use linfuncs
   use ode_s_solvers
   use sde_solvers
   use envelopes
   use signals
   use fclass
   use equation_solvers
   use ode_solvers
   implicit none
   integer(ip) :: nql,ncl
   type(djch) :: Hdjc,Hsc
   type(lcg) :: lcg_t
   type(diff_schr) :: sfun
   real(fp),allocatable,dimension(:) :: wq
   complex(fp),allocatable,dimension(:) :: y0s
   integer(ip) :: ii,jj,kk,NN
   complex(fp),allocatable,dimension(:) :: d
   complex(fp),allocatable,dimension(:,:) :: s
   real(fp),allocatable,dimension(:) :: x,xquery,xsc
   complex(fp),allocatable,dimension(:,:),target :: y,yquery
   complex(fp),allocatable,dimension(:,:) :: y_step_split
   real(fp),allocatable,dimension(:,:) :: j
   complex(fp),allocatable,dimension(:,:,:) :: ysc,ysf
   complex(fp),allocatable,dimension(:,:) ::yso
   complex(fp),dimension(:,:,:),pointer :: rho
   type(signal) :: test_signal
   type(constant_envelope) :: cenv
   integer(ip) :: c1,c2
   type(fisrk15) :: f15
   real(fp) :: ts
   type(fdlmtd) :: ffun

   type(diff_sme) :: sme
   complex(fp),allocatable,dimension(:,:) :: y0m

#define TEST_SME

#ifdef TEST_SCHR

   nql = 75
   ncl = 1
   NN = 501

   allocate(x(NN))
   call linspace(0.0_fp,0.001_fp,NN,x)

   allocate(y0s(nql*ncl))
   y0s = (0.0_fp,0.0_fp)
   y0s(2) = 1.0_fp

   allocate(wq(nql))
   wq = 0.0_fp
   do ii=1,size(wq)
      wq(ii) = (ii-1)**0.7
   end do

   allocate(d(nql*ncl))
   allocate(s(nql*ncl,2))

   allocate(y(nql*ncl,NN))
   allocate(j(1,NN))

   call Hdjc%initialize(nql,ncl,wq,0.8_fp,0.01_fp,hbar_in = 1.0_fp)
   call test_signal%initialize(Hdjc%op%cc(:,:,1),1.0_fp)
   call cenv%initialize((0.05_fp,0.0_fp))
   call test_signal%add_envelope(cenv)
   call Hdjc%add_signal(test_signal)
   !write(*,*) Hdjc%value_at(0.0_fp)
   call sfun%initialize(Hdjc)
   call sfun%add_observable(Hdjc%op%aa(:,:,1),(0.0_fp,0.0_fp))
   call sfun%add_observable(Hdjc%op%aa(:,:,1),(0.0_fp,0.0_fp))
   !call sfun%value_at(0.0_fp,y0s,d,s)
   !call lcg_t%init(generate_random_seed())
   call lcg_t%init(0)
   !call test_sde_solver(x,y0s,sfun,y,j,lcg_t,0.01_fp)
   !call platen(x,y0s,sfun,y,j,lcg_t,0.2_fp)
   !call f15%initialize(x(1),y0s,reshape([0.0_fp,0.2_fp],[1,2]),sfun,0.01_fp,sfun%nobs)
   !call f15%value_at(x(1),y0s,y(:,1))
   !write(*,*) 'y:',y(:,1)
   !call picard_iteration(f15,y0s,y(:,1),1e-15_fp)
   !write(*,*) 'Picard iteration, y:'
   !write(*,*) y0s
   !write(*,*) y(:,1)
   !call isrk15(x,y0s,sfun,y,j,lcg_t,1.0_fp)
   !call euler_maryama(x,y0s,sfun,y,j,lcg_t,1.0_fp)
   !call platen(x,y0s,sfun,y,j,lcg_t,1.0_fp)

   
   !allocate(xquery(NN/3),yquery(size(y0s),NN/3))
   !call linspace(0.1_fp,290.5_fp,NN/3,xquery)
   !yquery = interp1(x,y,xquery)
  
   !do ii=1,size(x)
   !   y(:,ii) = y(:,ii)/l2norm(y(:,ii))
   !end do
  
   ! open(unit = 1, file = 'st_data_s.dat')
   ! open(unit = 2, file = 'st_data_s_interp.dat')
   ! do ii=1,size(x)
   !    write(1,*) x(ii), real(conjg(y(:,ii))*y(:,ii))
   ! end do
   ! do ii=1,size(xquery)
   !    write(2,*) xquery(ii),real(conjg(yquery(:,ii))*yquery(:,ii))
   ! end do
   ! close(1)
   ! close(2)

   ! call system_clock(c1)
   ! !call basic_ode_s_solver(x,y0s,sfun,y,euler_maryama,100,0.1_fp)
   ! !call basic_ode_s_solver(x,y0s,sfun,y,platen,1000,0.2_fp)
   ! call system_clock(c2)
   ! write(*,*) 'Linear', (c2-c1), 'milliseconds'
   ! call system_clock(c1)
   ! !call basic_ode_s_solver(x,y0s,sfun,y,euler_maryama,100,0.1_fp)
   ! !call basic_ode_s_solver(x,y0s,sfun,y,platen,100,2.0_fp)
   ! !call parallel_ode_s_solver(x,y0s,sfun,y,platen,40,0.2_fp)
   ! call system_clock(c2)
   ! write(*,*) 'Parallel', (c2-c1), 'milliseconds'


   ! open(unit = 1, file = 'st_data_ode_s.dat')
   ! do ii=1,size(x)
   !    write(1,*) x(ii), real(y(:,ii))
   ! end do
   ! close(1)

    ! allocate(y_step_split(NN,4))
    ! do ii=1,4
    !    call lcg_t%init(0)
    !    call isrk15(x,y0s,sfun,y,j,lcg_t,0.5_fp,2*(ii-1))
    !    !call platen(x,y0s,sfun,y,j,lcg_t,0.5_fp,2*(ii-1))
    !    y_step_split(:,ii) = y(1,:)
    ! end do
  
    ! open(unit = 1, file = 'step_split.dat')
    ! do ii=1,size(x)
    !    write(1,*) x(ii), real(conjg(y_step_split(ii,:))*y_step_split(ii,:))
    ! end do
    ! close(1)

   ! Comparison of different solvers

   ! ts = 0.1_fp
   ! allocate(ysc(nql*ncl,NN,4))
   ! ysc = (0.0_fp,0.0_fp)
   ! ! write(*,*) 'Solver comparison:'
   ! ! call lcg_t%init(0)
   ! ! call system_clock(c1)
   ! ! call euler_maryama(x,y0s,sfun,ysc(:,:,1),j,lcg_t,ts,2)
   ! ! call system_clock(c2)
   ! ! write(*,*) 'Euler-Maryama', (c2-c1), 'milliseconds'
   ! ! call lcg_t%init(0)
   ! ! call system_clock(c1)
   ! ! call platen(x,y0s,sfun,ysc(:,:,2),j,lcg_t,ts,2)
   ! ! call system_clock(c2)
   ! ! write(*,*) 'Platen', (c2-c1), 'milliseconds'
   ! call lcg_t%init(0)
   ! call system_clock(c1)
   ! call isrk15(x,y0s,sfun,ysc(:,:,3),j,lcg_t,ts,2)
   ! call system_clock(c2)
   ! write(*,*) 'ISRK15', (c2-c1), 'milliseconds'
   !  ! call lcg_t%init(0)
   !  ! call system_clock(c1)
   !  ! call isrk15(x,y0s,sfun,ysc(:,:,4),j,lcg_t,0.1*ts)
   !  ! call system_clock(c2)
   !  ! write(*,*) 'Reference finished in', (c2-c1), 'milliseconds'

   ! open(unit = 1, file = 'solver_comparison.dat')
   ! do ii=1,size(x)
   !    write(1,*) x(ii), real(conjg(ysc(1,ii,:))*ysc(1,ii,:))
   ! end do
   ! close(1)

   ! Stochastic vs ode
   allocate(yso((nql*ncl),NN))
   allocate(ysf((nql*ncl)**2,NN,2))
   yso = (0.0_fp,0.0_fp)

   call ffun%initialize_fdlmtd(Hdjc)
   call ffun%add_lindblad_operator_fdlmtd(Hdjc%op%aa(:,:,1),(0.0_fp,0.0_fp))
   call lcg_t%init(0)
   call system_clock(c1)
   call isrk15(x,y0s,sfun,yso,j,lcg_t,0.05_fp)
   call system_clock(c2)
   write(*,*) 'ISRK15', (c2-c1), 'milliseconds'
   ! call system_clock(c1)
   ! call rk4(x,reshape(matmul(conjg(reshape(y0s,[size(y0s),1])),reshape(y0s,[1,size(y0s)])),[size(yso,1)**2]),&
   !      ffun,ysf(:,:,2),0.05_fp,0.01_fp)
   ! call system_clock(c2)
   ! write(*,*) 'rk45', (c2-c1), 'milliseconds'

   do ii=1,NN
      ysf(:,ii,1) = reshape(matmul(conjg(reshape(yso(:,ii),[size(y0s),1])),reshape(yso(:,ii),[1,size(y0s)])),[size(y0s,1)**2])
   end do
   
   open(unit = 1, file = 'stode_comp.dat')
   do ii=1,size(x)
      write(1,*) x(ii), real(ysf(1,ii,1)),real(ysf(1,ii,2))
   end do
   close(1)

#endif

#ifdef TEST_SME

   nql = 2
   ncl = 1
   NN = 10001

   if(allocated(x)) then
      deallocate(x)
   else
      allocate(x(NN))
   end if
   call linspace(0.0_fp,50.0_fp,NN,x)

   allocate(y0m(nql*ncl,nql*ncl))
   allocate(y0s(nql*ncl))
   y0s = (0.0_fp,0.0_fp)
   y0s(1) = 0.0_fp
   y0s(2) = 1.0_fp
   y0s = y0s/l2norm(y0s)
   y0m = reshape(matmul(reshape(y0s,[size(y0s),1]),reshape(y0s,[1,size(y0s)])),[ncl*nql,ncl*nql])
   call print_matrix(y0m)
   !y0m = (0.0_fp,0.0_fp)
   !y0m(1,1) = 0.9_fp
   !y0m(2,2) = 0.1_fp

   allocate(wq(nql))
   wq = 0.0_fp
   do ii=1,size(wq)
      wq(ii) = (ii-1)**0.7
   end do

   !allocate(d(nql*ncl*nql*ncl))
   !allocate(s(nql*ncl*nql*ncl,1))
   allocate(y(nql*ncl*nql*ncl,NN))
   allocate(j(1,NN))

   y = 0.0_fp
   j = 0.0_fp
   call lcg_t%init(4)
   
   call Hdjc%initialize(nql,ncl,wq,0.8_fp,0.05_fp,hbar_in = 1.0_fp,initial_sg_cont_size = 5,wref_in=1.0_fp)
   call test_signal%initialize(Hdjc%op%cc(:,:,1),1.0_fp)
   call cenv%initialize((0.0_fp,0.0_fp))
   call test_signal%add_envelope(cenv)
   call Hdjc%add_signal(test_signal)
   call sme%initialize(Hdjc)
   !call sme%add_lindblad_operator(Hdjc%op%aa(:,:,1),(5.0_fp,0.0_fp))
   call sme%add_observable(Hdjc%op%aa(:,:,1),(0.2_fp,0.0_fp))
   !call platen(x,reshape(y0m,[nql**2*ncl**2]),sme,y,j,lcg_t,0.0001_fp,step_split_factor_in = 1)
   call isrk15(x,reshape(y0m,[nql**2*ncl**2]),sme,y,j,lcg_t,0.0001_fp,step_split_factor_in = 1)

   rho(1:(ncl*nql),1:(ncl*nql),1:NN) => y
   
   open(unit = 1, file = 'sme.dat')
   do ii=1,NN
      write(1,'(F7.4)',advance="no") x(ii)
      do jj=1,ncl*nql
         !write(*,*) (jj-1)*(nql*ncl)+jj
         write(1,'(F7.4)',advance="no") real(y((jj-1)*(nql*ncl)+jj,ii))
      end do
      write(1,*) ''
   end do
   close(1)
   open(unit = 1, file = 'sme_ex.dat')
   do ii=1,NN
      write(1,'(F8.4)',advance="no") x(ii)
      write(1,'(F8.4)',advance="no") real(trace(matmul(rho(:,:,ii),Hdjc%op%nn(:,:,1))))*2-1
      write(1,'(F8.4)',advance="no") real(trace(matmul(rho(:,:,ii),Hdjc%op%sx(:,:,1))))
      write(1,'(F8.4)',advance="no") real(trace(matmul(rho(:,:,ii),Hdjc%op%sy(:,:,1))))
      write(1,*) ''
   end do
   close(1)
#endif


 end program test_sde_solvers
