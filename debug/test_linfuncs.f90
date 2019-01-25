program test_linfuncs
  use linfuncs
  use hamiltonians
  use fclass
  use signals
  use envelopes
  use equation_solvers
  use ode_solvers
  use omp_lib
  implicit none
  type(djch) :: H,Hg
  type(signal) ::sign,signg
  type(constant_envelope) :: cenv,cenvg
  type(fdlmtd) :: ffun,fung
  type(ftraceconst) :: ftc
  type(fggm2dens) :: fggm
  type(fsuplin) :: suplin
  type(djch),dimension(:),allocatable,target :: Hg_array
  integer(ip) :: nl
  complex(fp),dimension(5) :: fun_value,rho
  integer(ip) :: dim,ii
  complex(fp), dimension(:,:,:),allocatable :: ggm
  real(fp),dimension(:),allocatable :: ggmnorm
  complex(fp),dimension(:,:), allocatable,target :: rhog,rhor
  complex(fp),dimension(:),allocatable :: rhog_lin,rhor_lin
  complex(fp),dimension(:,:),allocatable :: rhod
  real(fp),dimension(:), allocatable :: bg
  real(fp),dimension(:), allocatable :: wg
  complex(fp),dimension(:), allocatable :: b0,b1
  integer(ip) :: c1,c2
  real(fp),allocatable,dimension(:) :: t
  integer(ip) :: Nts

#ifdef _OPENMP
  write(*,*) omp_get_max_threads()
#else
  write(*,*) 'No openMP'
#endif
  ! call H%initialize(2,1,[-1.0_fp,1.0_fp],0.0_fp,0.0_fp,hbar_in = 1.0_fp)
  ! call sign%initialize(H%op%cc(:,:,1),1.0_fp)
  ! call cenv%initialize((1.0_fp,1.0_fp))
  ! call sign%add_envelope(cenv)
  ! call H%add_signal(sign)
  ! !call print_matrix(reshape(H%value_at(0.0_fp),[2,2]))

  ! rho = (0.0_fp,0.0_fp)
  ! rho(1) = 1.0_fp
  ! rho(4) = 1.0_fp
  ! rho(5) = (0.0_fp,0.0_fp)

  
  
  ! !rho(2) = (0.5_fp,0.5_fp)
  ! !rho(3) = (0.5_fp,-0.5_fp)
  ! call ffun%initialize_fdlmtd(H)
  ! call ffun%add_lindblad_operator_fdlmtd(H%op%aa(:,:,1),(1.0_fp,1.0_fp))

  ! call print_matrix(matmul(reshape(rho(1:4),[2,2]),H%value_at(0.0_fp)))
 
  
  ! !call ffun%value_at(0.0_fp,rho,fun_value)
  ! !call print_matrix(reshape(fun_value,[2,2]))
  ! !call print_matrix(jacobian_approximation_complex(ffun,rho,0.1_fp))

  ! call ftc%initialize(ffun,ffun%H%dim)

  ! call ftc%value_at(0.0_fp,rho,fun_value)
  ! write(*,*) fun_value
  ! !call print_matrix(jacobian_approximation_complex(ftc,rho,0.8_fp))

  
  ! call quasi_newton(ftc,rho,rho,tol_in = 1e-3_fp,dx_in = 0.0001_fp)
  !write(*,*) rho

  !dim = 2
  
  !allocate(ggm(dim,dim,dim**2-1))
  !allocate(rhog(dim,dim))
  !allocate(bg(dim**2-1))

  ! rhog = (0.0_fp,0.0_fp)
  ! rhog(1,1) = 0.5_fp
  ! rhog(1,2) = 0.5_fp
  ! rhog(2,1) = 0.5_fp
  ! rhog(2,2) = 0.5_fp

  ! call print_matrix(rhog)
  
  ! call generate_ggm(dim,ggm)
  ! call density_matrix_to_ggm(rhog,ggm,bg)
  ! write(*,*) 'bg:',bg
  ! call ggm_to_density_matrix(rhog,ggm,bg,dim)
  ! call print_matrix(rhog)

  dim = 2
  Nts = 101
  allocate(wg(dim))
  wg = 0.0_fp
  do ii=2,dim
     wg(ii) = wg(ii-1) + 0.7_fp**(ii-2)
  end do
  
  allocate(Hg_array(1))
  call Hg_array(1)%initialize(dim,1,wg*0,0.0_fp,0.0_fp,hbar_in = 1.0_fp)

  call Hg%initialize(dim,1,wg*0,0.0_fp,0.0_fp,hbar_in = 1.0_fp)
  call signg%initialize(Hg%op%cc(:,:,1),1.0_fp)
  call cenvg%initialize((0.0_fp,0.0_fp))
  call signg%add_envelope(cenvg)
  !call Hg%add_signal(signg)
  call Hg_array(1)%add_signal(signg)

  call fung%initialize(Hg_array(1))
  !call fung%initialize(Hg)
  call fung%add_lindblad_operator(Hg%op%aa(:,:,1),(1.0_fp,0.0_fp))
  call fung%add_lindblad_operator(Hg%op%cc(:,:,1),(1.0_fp,0.0_fp))

  call fggm%initialize(fung,dim)

  allocate(b0(dim**2),b1(dim**2))
  allocate(ggm(dim,dim,dim**2))
  allocate(ggmnorm(dim**2))
  allocate(rhog(dim,dim))
  allocate(rhor(dim,dim))
  allocate(rhod(dim**2,Nts))
  allocate(t(Nts))
  allocate(bg(dim**2))

  call linspace(0.0_fp,100.0_fp,Nts,t)

  rhog = (0.0_fp,0.0_fp)
  !do ii=1,dim
  !   rhog(ii,ii) = 1.0_fp/dim
  !end do
  rhog(1,1) = 0.5_fp
  rhog(1,2) = 0.5_fp
  rhog(2,1) = 0.5_fp
  rhog(2,2) = 0.5_fp
  !rhog(3,3) = 0.0_fp
  !rhog(4,4) = 0.0_fp
  !rhog(10,10) = 1.0_fp

  call generate_ggm(dim,ggm,ggmnorm)
  call density_matrix_to_ggm(rhog,ggm,bg)
  b0 = bg
  write(*,*) 'b0:',real(b0)

  ! call fggm%value_at(0.0_fp,b0,b1)
  ! write(*,*) 'b1:',b1
  ! !b0(1) = b0(1) + 0.001_fp
  ! call fggm%value_at(0.0_fp,b0,b1)
  ! write(*,*) 'b1:',b1


  !call quasi_newton_real(fggm,b0,b1,tol_in = 1e-3_fp,dx_in = 0.1_fp)
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

  call suplin%initialize(fung,ggm,ggmnorm,dim)
  call suplin%value_at(0.0_fp,b0,b1)
  !call ggm_to_density_matrix(rhor,ggm,ggmnorm,real(b1),dim)
  !call print_matrix(rhor)
  write(*,*) real(b1)

  allocate(rhog_lin(dim**2),rhor_lin(dim**2))
  rhog_lin = reshape(rhog,[dim**2])
  call fung%value_at(0.0_fp,rhog_lin,rhor_lin)
  rhor = reshape(rhor_lin,[dim,dim])
  call density_matrix_to_ggm(rhor,ggm,bg)
  call print_matrix(rhor)
  write(*,*) bg

  write(*,*) 'L:'
  call print_matrix(real(suplin%L))

end program test_linfuncs
