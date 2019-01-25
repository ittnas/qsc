program test_assignment
  use linfuncs
  use hamiltonians
  use fclass
  use signals
  use envelopes
  use equation_solvers
  use ode_solvers
  implicit none
  type(djch) :: H,Hg
  type(signal) ::sign,signg
  type(constant_envelope) :: cenv,cenvg
  type(fdlmtd) :: ffun,fung
  type(ftraceconst) :: ftc
  type(fggm2dens) :: fggm
  type(fdsti) :: fdsti1,fdsti2
  type(fdlmtd) :: fdlmtd1,fdlmtd2
  integer(ip) :: nl
  integer(ip) :: dim,ii
  complex(fp), dimension(:,:,:),allocatable :: ggm
  real(fp),dimension(:),allocatable :: ggmnorm
  complex(fp),dimension(:,:), allocatable :: rhog,rhor
  complex(fp),dimension(:,:),allocatable :: rhod
  real(fp),dimension(:), allocatable :: bg
  real(fp),dimension(:), allocatable :: wg
  complex(fp),dimension(:), allocatable :: b0,b1
  integer(ip) :: c1,c2
  real(fp),allocatable,dimension(:) :: t
  integer(ip) :: Nts
  complex(fp),dimension(:),allocatable :: rho0,rho

  dim = 2
  allocate(wg(dim))
  wg = 0.0_fp
  do ii=2,dim
     wg(ii) = wg(ii-1) + 0.7_fp**(ii-2)
  end do
  
  call Hg%initialize(dim,1,wg,0.0_fp,0.0_fp,hbar_in = 1.0_fp)
  call signg%initialize(Hg%op%cc(:,:,1),1.0_fp)
  call cenvg%initialize((1.0_fp,0.0_fp))
  call signg%add_envelope(cenvg)
  call Hg%add_signal(signg)

  call fdsti1%initialize(Hg%value_at(0.0_fp))
  call fdlmtd1%initialize(Hg)
  
  allocate(rho0(dim))
  allocate(rho(dim))
  rho0 = (0.0_fp,0.0_fp)
  rho0(2) = 1.0_fp

  fdsti2 = fdsti1
  call fdsti1%value_at(0.0_fp,rho0,rho)
  call fdsti2%value_at(0.0_fp,rho0,rho)
  !write(*,*) rho
  
  call fdsti1%delete()
  
  !write(*,*) fdlmtd2%H%value_at(0.0_fp)
  fdlmtd2 = fdlmtd1
  !call Hg%add_signal(signg)
  write(*,*) fdlmtd1%H%value_at(0.0_fp)
  write(*,*) fdlmtd2%H%value_at(0.0_fp)
  !call fdsti2%delete()

  ! call fung%initialize_fdlmtd(Hg)
  ! call fung%add_lindblad_operator(Hg%op%aa(:,:,1),(1.0_fp,0.0_fp))
  ! call fung%add_lindblad_operator(Hg%op%cc(:,:,1),(1.0_fp,0.0_fp))

  ! call fggm%initialize(fung,dim)

  ! allocate(b0(dim**2-1),b1(dim**2-1))
  ! allocate(ggm(dim,dim,dim**2-1))
  ! allocate(ggmnorm(dim**2-1))
  ! allocate(rhog(dim,dim))
  ! allocate(rhor(dim,dim))
  ! allocate(rhod(dim**2,Nts))
  ! allocate(t(Nts))
  ! allocate(bg(dim**2-1))

  ! call linspace(0.0_fp,100.0_fp,Nts,t)

  ! rhog = (0.0_fp,0.0_fp)
  ! !do ii=1,dim
  ! !   rhog(ii,ii) = 1.0_fp/dim
  ! !end do
  ! rhog(1,1) = 0.5_fp
  ! !rhog(1,2) = 0.5_fp
  ! !rhog(2,1) = 0.5_fp
  ! rhog(2,2) = 0.5_fp
  ! !rhog(10,10) = 1.0_fp

  ! call generate_ggm(dim,ggm,ggmnorm)
  ! call density_matrix_to_ggm(rhog,ggm,bg)
  ! b0 = bg
  ! !write(*,*) 'b0:',b0

  ! ! call fggm%value_at(0.0_fp,b0,b1)
  ! ! write(*,*) 'b1:',b1
  ! ! !b0(1) = b0(1) + 0.001_fp
  ! ! call fggm%value_at(0.0_fp,b0,b1)
  ! ! write(*,*) 'b1:',b1


  ! !call quasi_newton_real(fggm,b0,b1,tol_in = 1e-3_fp,dx_in = 0.1_fp)
  ! call system_clock(c1)
  ! !call broyden_real(fggm,b0,b1,tol_in = 1e-2_fp,dx_in = 0.1_fp,maxiter_in = 10000)
  ! call system_clock(c2)
  ! write(*,*) 'The duration:', c2-c1,'ms'
  ! call system_clock(c1)
  ! call spectral_real(fggm,b0,b1,tol_in = 1e-2_fp,dx_in = 0.1_fp,maxiter_in = 10000)
  ! call system_clock(c2)
  ! write(*,*) 'The duration:', c2-c1,'ms'
  ! !write(*,*) 'The result:', b1
  ! call ggm_to_density_matrix(rhor,ggm,ggmnorm,real(b1),dim)
  ! write(*,*) 'The purity of the system:',real(trace(matmul(rhor,rhor)))
  ! !write(*,*) 'The grand answer:'
  ! !write(*,*) real(diagonal(rhor))
  ! !call print_matrix(rhog)
  ! ! do ii=1,(dim**2-1)
  ! !    write(*,*) 'ii:',ii
  ! !    call print_matrix(ggm(:,:,ii))
  ! ! end do
  
  ! call rk45(t,reshape(rhog,[dim**2]),fung,rhod,0.001_fp,0.0001_fp)
  ! !write(*,*) 'The dynamical solver:'
  ! !write(*,*) real(diagonal(reshape(rhod(:,size(rhod,2)),[dim,dim])))

  ! write(*,*) 'The difference:', sum((real(diagonal(reshape(rhod(:,size(rhod,2)),[dim,dim]))) - real(diagonal(rhor)))**2)
  

end program test_assignment
