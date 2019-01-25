PROGRAM TEST_MAIN
   use utils, only: fp,ip
   use fclass
   use ode_solvers
   use linfuncs
   use hamiltonians
   use operators
   use signals
   use envelopes
   use random
   use constants, only: pi
   implicit none
   real(fp),allocatable,dimension(:) :: x,wq,xquery
   complex(fp),allocatable,dimension(:) :: y0s,y0m
   complex(fp),allocatable,dimension(:,:) :: ys,ym,yquery
   complex(fp) :: H(2,2)
   type(fdstd) :: ftd
   type(fdsti) :: fti
   type(fdsti_sparse) :: fti_sparse
   type(fdmtd) :: fmtd

   type(fdstd) :: ftd_driven
   integer(ip) :: ii,NN,c1,c2,nql,ncl,jj
   type(simple_hamiltonian) :: Hs
   type(jch) :: Hjc
   type(djch) :: Hdjc
   type(cpb) :: Hcpb
   integer(ip) :: test_dimensions(3)
   type(operator_cont) :: op
   procedure(ode_solver),pointer :: solver => null()
   type(signal) :: test_signal
   type(constant_envelope) :: cenv
   type(lcg) :: lcg_t
   real(fp),dimension(:),allocatable :: rnumber
   integer(dip),dimension(:),allocatable :: inumber
   integer(ip) :: nr
   type(fdlmtd) :: flmtd
   real,allocatable,dimension(:) :: Et

#define TEST_SPARSE
!#define TEST_SOLVERS
!#define TEST_RNG
  
   test_dimensions = [2,1,1]
   NN = 1001
   nql = 20
   ncl = 1

   solver => rk45
   !call op%initialize(test_dimensions)

   allocate(y0s(nql*ncl))
   allocate(y0m((nql*ncl)**2))
   y0s = (0.0_fp,0.0_fp)
   y0s(2) = 1.0_fp
   y0m = reshape(matmul(reshape(y0s,[size(y0s),1]),reshape(y0s,[1,size(y0s)])),[size(y0m)])
   !write(*,*) real(y0)
   !y0 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
   !allocate(y0(NH))
   allocate(x(NN))
   allocate(ys(size(y0s),NN))
   allocate(ym(size(y0m),NN))
   allocate(wq(nql))
   allocate(xquery(NN/3),yquery(size(y0s),NN/3))
   wq = 0.0_fp
   do ii=1,size(wq)
      wq(ii) = ii**0.7
   end do

   call linspace(0.0_fp,1000.0_fp,NN,x)
   call linspace(0.0_fp,0.1_fp,NN/3,xquery)
   !  f = fd(H)

   !call Hs%initialize_simple_H(2,0.0_fp,1.0_fp)
   call Hdjc%initialize(nql,ncl,wq,0.5_fp,1.0_fp,hbar_in = 1.0_fp,initial_sg_cont_size = 5)
   call test_signal%initialize(Hdjc%op%cc(:,:,1),1.0_fp)
   call cenv%initialize((0.01_fp,0.0_fp))
   call test_signal%add_envelope(cenv)
   call Hdjc%add_signal(test_signal)

   call Hjc%initialize(nql,ncl,wq,0.5_fp,1.0_fp,hbar_in = 1.0_fp)
   call fti%initialize(Hjc%value_at(0.0_fp))
   call fti_sparse%initialize(Hjc%value_at(0.0_fp))
   call ftd%initialize(Hjc)
   call ftd_driven%initialize(Hdjc)
   call fmtd%initialize(Hjc)

   call flmtd%initialize(Hdjc)
   !call flmtd%add_lindblad_operator_fdlmtd(Hjc%op%aa(:,:,1),(0.01_fp,0.0_fp))
   !call flmtd%add_lindblad_operator_fdlmtd(cmplx(get_diag_matrix(wq),kind=fp),(0.5_fp,0.0_fp))
#ifdef TEST_SOLVERS
   call system_clock(c1)
   !call euler(x,y0m,fmtd,ym,0.00001_fp)
   !call euler(x,y0m,fmtd,ym,0.001_fp)
   call solver(x,y0m,flmtd,ym,0.01_fp,0.01_fp)
   call system_clock(c2)
   write(*,*) 'DM execution took', (c2-c1), 'milliseconds'
   call system_clock(c1)
   !call euler(x,y0s,ftd,ys,0.00001_fp)
   !call euler(x,y0s,ftd,ys,0.001_fp)
   call rk4(x,y0s,ftd,ys,0.001_fp,0.01_fp)
   call system_clock(c2)
   write(*,*) 'SE execution with Hjc took', (c2-c1), 'milliseconds'
   call system_clock(c1)
   call rk4(x,y0s,ftd_driven,ys,0.001_fp,0.01_fp)
   call system_clock(c2)
   write(*,*) 'SE execution with Hdjc took', (c2-c1), 'milliseconds'

  
   !call print_matrix(Hjc%value_at(0.0_fp),.true.)
   !call print_matrix(Hjc%qcca,.true.)
   !call print_matrix(Hjc%qacc,.true.)
   !write(*,*) real(CONJG(y(1,:))*y(1,:))
  
   yquery = interp1(x,ys,xquery)
  
   open(unit = 1, file = 'test_data_s.dat')
   open(unit = 2, file = 'test_data_m.dat')
   open(unit = 3, file = 'test_data_q.dat')
   do ii=1,size(x)
      !write(1,*) x(ii),real(CONJG(y(:,ii))*y(:,ii))
      write(1,*) x(ii),real(conjg(ys(:,ii))*ys(:,ii))
      !write(2,*) x(ii),real(ym(:,ii))
      write(2,"(f8.3)",advance='no') x(ii)
      do jj=1,nql*ncl
         write(2,"(f8.3)",advance='no') ,real(ym(nql*ncl*(jj-1) + jj,ii))
      end do
      write(2,*) ''
   end do
   do ii=1,size(xquery)
      write(3,*) xquery(ii),real(conjg(yquery(:,ii))*yquery(:,ii))
   end do
   close(1)
   close(2)
   close(3)
#endif
#ifdef TEST_RNG
   !write(*,*) Hdjc%value_at(1.0_fp)
   call lcg_t%init(generate_random_seed())
   !rnumber = lcg_t%uniform(1)
   !write(*,*) rnumber
   !nr = 20
   !allocate(inumber(nr),rnumber(nr))
   !inumber = lcg_t%uniform_int(nr)
   !rnumber = lcg_t%norm(nr)
   !write(*,*) rnumber
   !call print_matrix(Hjc%value_at(0.0_fp))
   allocate(Et(5))
   Et = Hdjc%get_eigenmodes(0.0_fp,5)
   !write(*,*) 'Eigenenergies: ', Et
   !write(*,*) 'wq: ', wq

   call Hcpb%initialize(1.0_fp,50.0_fp,0.0_fp,PI*2,51,4) ! Eventhough in transmon regime, gives different values for ng = {even,odd}. Maybe this has something to do with ng/2 in the cpb Hamiltonian. CHECK!
   Et = Hcpb%get_eigenmodes(0.0_fp,5)
   write(*,*) 'Eigenmodes of cpb: ', Et
#endif
! ------------------------------------------ !
   ! Sparse vs non sparse
#ifdef TEST_SPARSE
   write(*,*) 'Comparing sparse and dense solvers.'
   ys = 0.0_fp
   call system_clock(c1)
   do ii=1,100000
      call fti%value_at(0.0_fp,y0s,ys(:,1))
   end do
   call system_clock(c2)
   write(*,*) 'Dense loop took', (c2-c1), 'milliseconds'
   ys = 0.0_fp
   call system_clock(c1)
   do ii=1,100000
      call fti_sparse%value_at(0.0_fp,y0s,ys(:,1))
   end do
   call system_clock(c2)
   write(*,*) 'Sparse loop took', (c2-c1), 'milliseconds'
!#ifdef DONT
   ys = 0.0_fp
   call system_clock(c1)
   call rk4(x,y0s,fti,ys,0.01_fp,0.01_fp)
   call system_clock(c2)
   write(*,*) 'Dense execution took', (c2-c1), 'milliseconds'
   ys = 0.0_fp
   call system_clock(c1)
   call rk4(x,y0s,fti_sparse,ys,0.01_fp,0.01_fp)
   call system_clock(c2)
   write(*,*) 'Sparse execution took', (c2-c1), 'milliseconds'
!#endif
#endif
END PROGRAM TEST_MAIN
