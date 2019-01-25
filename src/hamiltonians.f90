module hamiltonians
  use utils, only: fp,ip
  use operators
  use signals
  implicit none
  private

  public :: hamiltonian,simple_hamiltonian,jch,djch,cpb,djch_tlf,cpb_tl,disp_jch,jch_qubit,trq,djch_tlf6

  
  type, abstract :: hamiltonian
     integer(ip) :: dim
   contains
     procedure,private :: initialize_hamiltonian
     procedure(value_at_interface),deferred :: value_at
     procedure :: get_eigenmodes
     generic,public :: initialize => initialize_hamiltonian
  end type hamiltonian

  

  ! Returns the value of the Hamiltonian times -i/hbar. This convention tries to avoid multiplication and division by hbar, which is a very small number and might thus introduce a lot of numerical errors. The multiplication by -i is just for additional efficiency, when integrating the the master equation rho' = -i/hbar*[H,rho].
  abstract interface
     function value_at_interface(this,t) result(H_out)
       use utils, only:fp,ip
       import hamiltonian
       class(hamiltonian), intent(in) :: this
       real(fp), intent(in) :: t
       complex(fp),dimension(this%dim,this%dim) :: H_out
     end function value_at_interface
  end interface

  type :: signal_container
     class(signal),pointer :: obj
  end type signal_container

  type, extends(hamiltonian) :: simple_hamiltonian
     real(fp) :: E,g
   contains
     procedure :: value_at=>value_at_simple_h
     procedure,private ::initialize_simple_h
     generic,public ::  initialize => initialize_simple_h
  end type simple_hamiltonian

  type, extends(hamiltonian) :: driven_hamiltonian
     type(signal_container),dimension(:),allocatable :: signals
     integer(ip) :: nsg,sg_cont_size
     real(fp) :: wref
     real(fp), allocatable, dimension(:) :: wref_array
     logical :: use_wref_array
   contains
     procedure :: value_at=>value_at_driven_hamiltonian
     procedure,private :: initialize_driven_hamiltonian
     generic,public :: initialize=>initialize_driven_hamiltonian
     procedure :: add_signal => add_signal_driven_hamiltonian
     procedure :: set_wref => set_wref_driven_hamiltonian
     procedure :: set_wref_array => set_wref_array_driven_hamiltonian
  end type driven_hamiltonian

  type,extends(driven_hamiltonian) :: jch_qubit
     class(operator_cont),allocatable :: op
     complex(fp) :: alpha_e,alpha_g,beta
     real(fp) :: Gamma_d,Gamma_ci,wac
     real(fp) :: wq,eps_d,wc,wm,kappa,g,eta
   contains
     procedure :: value_at=>value_at_jch_qubit
     procedure,private :: initialize_jch_qubit
     generic,public :: initialize=>initialize_jch_qubit
     procedure :: set_eps_d => set_eps_d_jch_qubit
  end type jch_qubit

  type, extends(hamiltonian) :: jch
     !real(fp),dimension(:),allocatable :: wq
     real(fp) :: wc,hbar
     real(fp),dimension(:),allocatable :: g
     integer(ip) :: nql,ncl
     class(operator_cont),allocatable :: op
     real(fp),dimension(:,:),allocatable :: wq_mat,qcca,qacc
   contains
     procedure :: value_at=>value_at_jch
     procedure,private :: initialize_jch
     generic,public :: initialize => initialize_jch
     !procedure :: delete ! The resources must be freed!
  end type jch

  type,extends(jch) :: djch
     type(signal_container),dimension(:),allocatable :: signals
     integer(ip) :: nsg,sg_cont_size
     real(fp) :: wref
     real(fp), allocatable, dimension(:) :: wref_array
     logical :: use_wref_array
   contains
     procedure :: value_at=>value_at_djch
     procedure,private :: initialize_djch
     generic,public :: initialize => initialize_djch
     procedure :: add_signal => add_signal_djch
     procedure :: set_wref => set_wref_djch
     procedure :: set_wref_array => set_wref_array_djch
     !procedure :: delete ! The resources must be freed!
  end type djch

  type,extends(hamiltonian) :: disp_jch
     type(signal_container),dimension(:),allocatable :: signals
     integer(ip) :: nsg,sg_cont_size
     real(fp) :: wref
     real(fp), allocatable, dimension(:) :: wref_array
     logical :: use_wref_array
     real(fp) :: wc,hbar
     real(fp),dimension(:),allocatable :: g
     integer(ip) :: nql,ncl
     class(operator_cont),allocatable :: op
     real(fp),dimension(:,:),allocatable :: wq_mat,qcca,qacc
   contains
     procedure :: value_at=>value_at_disp_jch
     procedure,private :: initialize_disp_jch
     generic,public :: initialize => initialize_disp_jch
     procedure :: add_signal => add_signal_disp_jch
     procedure :: set_wref => set_wref_disp_jch
     procedure :: set_wref_array => set_wref_array_disp_jch
  end type disp_jch

  type,extends(djch) :: djch_tlf
     real(fp) :: gtlf
     integer(ip) :: ntlf
     real(fp),allocatable,dimension(:,:) :: wtlf_mat
   contains
     procedure :: value_at=>value_at_djch_tlf
     procedure,private :: initialize_djch_tlf
     generic,public :: initialize => initialize_djch_tlf
     !procedure :: delete ! The resources must be freed!
  end type djch_tlf

    type,extends(djch) :: djch_tlf6
     real(fp) :: v,Omega,eps,Delta
     real(fp),dimension(6,6) :: w_mat
   contains
     procedure :: value_at=>value_at_djch_tlf6
     procedure,private :: initialize_djch_tlf6
     generic,public :: initialize => initialize_djch_tlf6
     !procedure :: delete ! The resources must be freed!
  end type djch_tlf6

  type,extends(hamiltonian) :: cpb
     complex(fp),allocatable,dimension(:,:) :: H ! H/hbar
     integer(ip) :: ng
   contains
     procedure :: value_at=>value_at_cpb
     procedure,private ::initialize_cpb
     procedure,private :: delete_cpb
     generic,public ::  initialize => initialize_cpb
     generic,public :: delete => delete_cpb
  end type cpb

  type,extends(hamiltonian) :: cpb_tl
     complex(fp),allocatable,dimension(:,:) :: H ! H/hbar
     integer(ip) :: ng
   contains
     procedure :: value_at=>value_at_cpb_tl
     procedure,private ::initialize_cpb_tl
     generic,public ::  initialize => initialize_cpb_tl
     !procedure :: delete ! The resources must be freed!
  end type cpb_tl

  type,extends(hamiltonian) :: trq
     type(signal_container),dimension(:),allocatable :: signals
     integer(ip) :: nsg,sg_cont_size
     class(operator_cont),allocatable :: op
     real(fp) :: dq,da,db,hbar
     complex(fp) :: gab,gqa,gqb
     complex(fp) :: Aq,Aa,Ab
     integer(ip) :: na,nb
   contains
     procedure :: value_at=>value_at_trq
     procedure,private :: initialize_trq
     generic,public :: initialize => initialize_trq
     procedure :: add_signal => add_signal_trq
     procedure,public :: set_parameters_trq
  end type trq
  
contains
  subroutine initialize_hamiltonian(this)
    use utils, only:fp,ip
    class(hamiltonian) :: this
  end subroutine initialize_hamiltonian

  ! returns the eigenmodes of the Hamiltonian. You get the eigenenergies from E = hbar*w. The reason for calculating the eigenmodes (= the angular frequencies) is that hbar might not be known for all Hamiltonians. This way, there should not be any problems with classical or quantum Hamiltonians.
  function get_eigenmodes(this,t,n) result(w)
    use linfuncs
    implicit none
    class(hamiltonian),intent(in) :: this
    integer(ip),intent(in) :: n
    real(fp),intent(in) :: t
    real(fp),dimension(n) :: w
    real(fp),dimension(this%dim) :: wt

    complex(fp), allocatable,dimension(:) :: WORK
    real(fp), allocatable,dimension(:) :: RWORK
    integer(ip) :: info,ii

    w = 0.0_fp
    ! wt = 0.0_fp

    !call print_matrix(this%value_at(t))
    !write(*,*) 'Dim: ', this%dim
    !call print_matrix(this%value_at(t))
    
    allocate(WORK(this%dim*5)) ! Work array for the eigenvalue solver
    allocate(RWORK(max(1,3*this%dim-2)))
    call zheev('V','L',this%dim,(0.0_fp,1.0_fp)*this%value_at(t),this%dim,wt,WORK,this%dim*5,RWORK,INFO)
    if(INFO .ne. 0) then
       write(*,*) 'Eigenvalue solver failed with error code ', INFO
    end if
    do ii=1,min(n,this%dim)
       w(ii) = wt(ii) - wt(1)
    end do
  end function get_eigenmodes

  function value_at_simple_h(this,t) result(H_out)
    use utils, only:fp,ip
    class(simple_hamiltonian), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp),dimension(this%dim,this%dim) :: H_out

    H_out(1,1) = -this%E/2_fp
    H_out(2,2) = this%E/2_fp
    H_out(1,2) = this%g
    H_out(2,1) = this%g
  end function value_at_simple_h

  subroutine initialize_simple_h(this,n,E,g)
    use utils, only:fp,ip
    class(simple_hamiltonian) :: this
    real(fp), intent(in) :: E,g
    integer(ip), intent(in) :: n
    this%dim = n
    this%E = E
    this%g = g
  end subroutine initialize_simple_h

  function value_at_jch(this,t) result(H_out)
    use utils, only:fp,ip
    use linfuncs
    class(jch), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp),dimension(this%dim,this%dim) :: H_out
    !H_out = (0.0_fp,0.0_fp) ! slows down performance
    !H_out = -(0.0_fp,1.0_fp)*(this%hbar*(this%wq_mat + this%wc*this%op%nn(:,:,2)) + this%qcca + this%qacc)
    H_out = -(0.0_fp,1.0_fp)*((this%wq_mat + this%wc*this%op%nn(:,:,2)) + this%qcca + this%qacc)
  end function value_at_jch

  subroutine initialize_jch(this,nql_in,ncl_in,wq_in,wc_in,g_in,garray_in,hbar_in)
    use utils, only:fp,ip
    use constants, only: hbar
    use operators
    use linfuncs
    class(jch) :: this
    real(fp),dimension(:),intent(in) :: wq_in
    real(fp),intent(in) :: wc_in
    real(fp),intent(in),optional :: hbar_in,g_in
    real(fp),dimension(:),intent(in),optional :: garray_in
    integer(ip),intent(in) :: nql_in,ncl_in
    integer(ip) :: ii
    this%nql = nql_in
    this%ncl = ncl_in
    this%dim = this%nql*this%ncl
    if(allocated(this%op)) then
       deallocate(this%op)
    endif
    allocate(this%op)
    call this%op%initialize([this%nql,this%ncl])

    allocate(this%g(nql_in))
    this%g = 0.0_fp
    if(present(g_in)) then
       do ii=2,nql_in
          this%g(ii) = g_in
       end do
    elseif(present(garray_in)) then
       this%g(2:nql_in) = garray_in
    else
       do ii=2,nql_in
          this%g(ii) = 1.0_fp
       end do
    endif
    if(present(hbar_in)) then
       this%hbar = hbar_in
    else
       this%hbar = hbar
    endif
    !allocate(this%wq(nql_in))
    !this%wq = wq_in
    allocate(this%wq_mat(this%dim,this%dim))
    this%wq_mat = kron(get_diag_matrix(wq_in),get_eye(this%ncl))
    this%wc = wc_in
    allocate(this%qcca(this%dim,this%dim))
    this%qcca = matmul(kron(get_diag_matrix(this%g),get_eye(this%ncl)),matmul(this%op%cc(:,:,1),this%op%aa(:,:,2)))
    allocate(this%qacc(this%dim,this%dim))
    !call print_matrix(matmul(this%op%aa(:,:,1),this%op%cc(:,:,2)),.true.)
    !call print_matrix(kron(get_diag_matrix(this%g),get_eye(this%ncl)))
    this%qacc = matmul(matmul(this%op%aa(:,:,1),this%op%cc(:,:,2)),kron(get_diag_matrix(this%g),get_eye(this%ncl)))
  end subroutine initialize_jch
!##############################################################################
  function value_at_djch(this,t) result(H_out)
    use utils, only:fp,ip
    use linfuncs
    class(djch), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp),dimension(this%dim,this%dim) :: H_out
    integer(ip) :: ii
    !H_out = this%value_at_jch(t)
    if(this%use_wref_array) then
       H_out = -(0.0_fp,1.0_fp)*((this%wq_mat&
            + this%wc*this%op%nn(:,:,2) - get_diag_matrix(this%wref_array)) + this%qcca + this%qacc)
       !write(*,*) 'value_at_djch'
       !call print_matrix(H_out)
       !call print_matrix(get_diag_matrix(this%wref_array))
       !call print_matrix(this%wc*this%op%nn(:,:,2))
       !call print_matrix(this%wq_mat)

       do ii=1,this%nsg
          H_out = H_out - (0.0_fp,1.0_fp)*this%signals(ii)%obj%value_at(t,this%wref_array)
       end do
    else
       H_out = -(0.0_fp,1.0_fp)*((this%wq_mat - this%op%nn(:,:,1)*this%wref &
            + (this%wc-this%wref)*this%op%nn(:,:,2)) + this%qcca + this%qacc)
       do ii=1,this%nsg
          H_out = H_out - (0.0_fp,1.0_fp)*this%signals(ii)%obj%value_at(t,this%wref)
       end do
    end if
  end function value_at_djch
!--------------------------------------------------------------
  subroutine initialize_djch(this,nql_in,ncl_in,wq_in,wc_in,g_in,garray_in,hbar_in,wref_in,wref_array,initial_sg_cont_size)
    use utils, only:fp,ip
    use constants, only: hbar
    use operators
    use linfuncs
    implicit none
    class(djch) :: this
    real(fp),dimension(:),intent(in) :: wq_in
    real(fp),intent(in) :: wc_in
    real(fp),intent(in),optional :: hbar_in,g_in,wref_in
    real(fp),dimension(:),intent(in),optional :: garray_in
    real(fp),dimension(:),intent(in), optional :: wref_array
    integer(ip),intent(in) :: nql_in,ncl_in
    !integer(ip),parameter :: initial_nbr_sg = 5
    integer(ip),intent(in) :: initial_sg_cont_size
    call initialize_jch(this,nql_in,ncl_in,wq_in,wc_in,g_in,garray_in,hbar_in)
    allocate(this%signals(initial_sg_cont_size))
    this%sg_cont_size = initial_sg_cont_size
    this%nsg = 0
    if(present(wref_in)) then
       this%wref = wref_in
    else
       this%wref = 0.0_fp
    end if
    
    allocate(this%wref_array(this%dim))
    if(present(wref_array)) then
       this%wref_array = wref_array
       this%use_wref_array = .TRUE.
    else
       this%use_wref_array = .FAlSE.
    endif
  end subroutine initialize_djch
!---------------------------------------------------------------
  subroutine add_signal_djch(this,signal_in)
    use utils,only:fp,ip
    use signals
    class(djch),intent(inout) :: this
    class(signal),intent(in),target :: signal_in
    type(signal_container),allocatable,dimension(:) :: temp_cont
    if(this%nsg >= this%sg_cont_size) then
       allocate(temp_cont(this%nsg))
       temp_cont = this%signals(1:this%nsg)
       deallocate(this%signals)
       allocate(this%signals(this%nsg*2))
       this%signals(1:this%nsg) = temp_cont
       deallocate(temp_cont)
       this%sg_cont_size = size(this%signals)
    end if
    this%nsg = this%nsg + 1
    this%signals(this%nsg)%obj => signal_in
  end subroutine add_signal_djch
!--------------------------------------------------------------------!
  subroutine set_wref_djch(this,wref)
    implicit none
    class(djch),intent(inout) :: this
    real(fp),intent(in) :: wref
    this%wref = wref
    this%use_wref_array = .FALSE.
  end subroutine set_wref_djch
!--------------------------------------------------------------------!
  subroutine set_wref_array_djch(this,wref_array)
    implicit none
    class(djch),intent(inout) :: this
    real(fp),dimension(:),intent(in) :: wref_array
    if(size(wref_array) .ne. this%dim) then
       write(*,*) 'Size of wref_array is incorrect.'
       stop(1)
    endif
    this%wref_array = wref_array
    this%use_wref_array = .TRUE.
  end subroutine set_wref_array_djch
!#######################################################################
  subroutine initialize_djch_tlf(this,nql_in,ncl_in,ntlf_in,wq_in,wc_in,wtlf,gtlf,g_in,garray_in,hbar_in,wref_in,wref_array)
    use utils, only:fp,ip
    use constants, only: hbar
    use operators
    use linfuncs
    class(djch_tlf) :: this
    real(fp),dimension(:),intent(in) :: wq_in
    real(fp),intent(in) :: wc_in
    real(fp),intent(in),optional :: hbar_in,g_in,wref_in
    real(fp),dimension(:),intent(in),optional :: garray_in
    real(fp),dimension(:),intent(in), optional :: wref_array
    integer(ip),intent(in) :: nql_in,ncl_in,ntlf_in
    real(fp),intent(in) :: gtlf
    real(fp),dimension(:),intent(in) :: wtlf
    integer(ip),parameter :: initial_nbr_sg = 5
    integer(ip) :: ii
    !call initialize_jch(this,nql_in,ncl_in,wq_in,wc_in,g_in,garray_in,hbar_in)
    ! Initialization of jch (with small changes for the tlf).
    this%nql = nql_in
    this%ncl = ncl_in
    this%ntlf = ntlf_in
    this%dim = this%nql*this%ncl*this%ntlf ! *2 for tlf
    if(allocated(this%op)) then
       deallocate(this%op)
    endif
    allocate(this%op)
    call this%op%initialize([this%nql,this%ncl,this%ntlf])
    allocate(this%g(nql_in))
    this%g = 0.0_fp
    if(present(g_in)) then
       do ii=2,nql_in
          this%g(ii) = g_in
       end do
    elseif(present(garray_in)) then
       this%g(2:nql_in) = garray_in
    else
       do ii=2,nql_in
          this%g(ii) = 1.0_fp
       end do
    endif
    if(present(hbar_in)) then
       this%hbar = hbar_in
    else
       this%hbar = hbar
    endif
    allocate(this%wq_mat(this%dim,this%dim))
    this%wq_mat = kron(kron(get_diag_matrix(wq_in),get_eye(this%ncl)),get_eye(this%ntlf))
    allocate(this%wtlf_mat(this%dim,this%dim))
    this%wtlf_mat = kron(kron(get_eye(this%nql),get_eye(this%ncl)),get_diag_matrix(wtlf))
    this%wc = wc_in
    allocate(this%qcca(this%dim,this%dim))
    this%qcca = matmul(kron(kron(get_diag_matrix(this%g),get_eye(this%ncl)),get_eye(this%ntlf))&
         ,matmul(this%op%cc(:,:,1),this%op%aa(:,:,2)))
    allocate(this%qacc(this%dim,this%dim))
    this%qacc = matmul(matmul(this%op%aa(:,:,1),this%op%cc(:,:,2)),&
         kron(kron(get_diag_matrix(this%g),get_eye(this%ncl)),get_eye(this%ntlf)))
    !-----------------------------
    !djc
    allocate(this%signals(initial_nbr_sg))
    this%sg_cont_size = initial_nbr_sg
    this%nsg = 0
    if(present(wref_in)) then
       this%wref = wref_in
    else
       this%wref = 0.0_fp
    end if

    allocate(this%wref_array(this%dim))
    if(present(wref_array)) then
       this%wref_array = wref_array
       this%use_wref_array = .TRUE.
    else
       this%use_wref_array = .FAlSE.
    endif
    !-----------------------------
    !djc_tlf
    !this%w_tlf = w_tlf
    this%gtlf = gtlf
  end subroutine initialize_djch_tlf
!--------------------------------------------------------------
  function value_at_djch_tlf(this,t) result(H_out)
    use utils, only:fp,ip
    use linfuncs
    class(djch_tlf), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp),dimension(this%dim,this%dim) :: H_out
    integer(ip) :: ii
    !H_out = this%value_at_jch(t)


    if(this%use_wref_array) then
    H_out = -(0.0_fp,1.0_fp)*((this%wq_mat - get_diag_matrix(this%wref_array) &
         + this%wc*this%op%nn(:,:,2)) + this%qcca + this%qacc &
         + this%wtlf_mat &
         + this%gtlf*(matmul(this%op%cc(:,:,1),this%op%aa(:,:,3)) &
         + matmul(this%op%aa(:,:,1),this%op%cc(:,:,3))))
       do ii=1,this%nsg
          H_out = H_out - (0.0_fp,1.0_fp)*this%signals(ii)%obj%value_at(t,this%wref_array)
       end do

       !write(*,*) 'wq + wtlf'
       !call print_matrix(this%wq_mat)
       !call print_matrix(this%wtlf_mat)
    else
    H_out = -(0.0_fp,1.0_fp)*((this%wq_mat - this%op%nn(:,:,1)*this%wref &
         + (this%wc-this%wref)*this%op%nn(:,:,2)) + this%qcca + this%qacc &
         + this%wtlf_mat - this%wref*this%op%nn(:,:,3) &
         + this%gtlf*(matmul(this%op%cc(:,:,1),this%op%aa(:,:,3)) &
         + matmul(this%op%aa(:,:,1),this%op%cc(:,:,3))))
       do ii=1,this%nsg
          H_out = H_out - (0.0_fp,1.0_fp)*this%signals(ii)%obj%value_at(t,this%wref)
       end do
    end if
  end function value_at_djch_tlf
!--------------------------------------------------------------
  !#######################################################################
    subroutine initialize_djch_tlf6(this,wq_in,v,eps,Delta,wref_in,wref_array)
    use utils, only:fp,ip
    use constants, only: hbar
    use operators
    use linfuncs
    class(djch_tlf6) :: this
    real(fp),dimension(:),intent(in) :: wq_in
    real(fp),intent(in),optional :: wref_in
    real(fp),dimension(:),intent(in), optional :: wref_array
    real(fp),intent(in) :: v, eps, Delta
    integer(ip),parameter :: initial_nbr_sg = 5
    integer(ip) :: ii
    real(fp) :: Omega0,Omega1,Omega2,eps0,eps1,eps2,theta0,theta1,theta2
    real(fp),dimension(6) :: Omega_vec
    real(fp) :: w0,w1,w0d,w1d
    !call initialize_jch(this,nql_in,ncl_in,wq_in,wc_in,g_in,garray_in,hbar_in)
    ! Initialization of jch (with small changes for the tlf).
    this%dim = 6 ! *2 for tlf
    if(allocated(this%op)) then
       deallocate(this%op)
    endif
    allocate(this%op)
    call this%op%initialize([3,2])
    allocate(this%wq_mat(this%dim,this%dim))
    !-----------------------------
    !djc
    allocate(this%signals(initial_nbr_sg))
    this%sg_cont_size = initial_nbr_sg
    this%nsg = 0
    if(present(wref_in)) then
       this%wref = wref_in
    else
       this%wref = 0.0_fp
    end if

    allocate(this%wref_array(this%dim))
    if(present(wref_array)) then
       this%wref_array = wref_array
       this%use_wref_array = .TRUE.
    else
       this%use_wref_array = .FAlSE.
    endif
    !-----------------------------
    !djc_tlf
    this%v = v
    this%eps = eps
    this%Delta = Delta

    eps0 = eps - 0.5_fp*v
    eps1 = eps - 1.5_fp*v
    eps2 = eps - 2.5_fp*v

    Omega0 = sqrt(eps0**2 + Delta**2)
    Omega1 = sqrt(eps1**2 + Delta**2)
    Omega2 = sqrt(eps2**2 + Delta**2)

    Omega_vec = [-Omega0,Omega0,-Omega1,Omega1,-Omega2,Omega2]

    if(Delta .eq. 0.0_fp .and. eps0-v .eq. 0.0_fp) then
       theta0 = 0.0_fp
    else
       theta0 = atan(Delta/(eps0 - v))
    end if
    if(Delta .eq. 0.0_fp .and. eps1-v .eq. 0.0_fp) then
       theta1 = 0.0_fp
    else
       theta1 = atan(Delta/(eps1 - v))
    end if
    if(Delta .eq. 0.0_fp .and. eps2-v .eq. 0.0_fp) then
       theta2 = 0.0_fp
    else
       theta2 = atan(Delta/(eps2 - v))
    end if

    if(isnan(theta0)) then
       theta0 = 0.0_fp
    end if

    if(isnan(theta1)) then
       theta1 = 0.0_fp
    end if

    if(isnan(theta2 )) then
       theta2 = 0.0_fp
    end if

!    write(*,*) theta0,theta1,theta2

    w0 = cos((theta1-theta0)*0.5_fp)
    w1 = cos((theta2-theta1)*0.5_fp)
    w0d = sin((theta1-theta0)*0.5_fp)
    w1d = sin((theta2-theta1)*0.5_fp)

    this%wq_mat = kron(get_diag_matrix(wq_in),get_eye(2)) &
         + get_diag_matrix(Omega_vec)

    this%w_mat = 0.0_fp
    this%w_mat(3,1) = w0
    this%w_mat(4,1) = w0d
    this%w_mat(3,2) = -w0d
    this%w_mat(4,2) = w0

    this%w_mat(5,3) = w1
    this%w_mat(6,3) = w1d
    this%w_mat(5,4) = -w1d
    this%w_mat(6,4) = w1

    this%w_mat(1,3) = w0
    this%w_mat(1,4) = w0d
    this%w_mat(2,3) = -w0d
    this%w_mat(2,4) = w0

    this%w_mat(3,5) = w1
    this%w_mat(3,6) = w1d
    this%w_mat(4,5) = -w1d
    this%w_mat(4,6) = w1
    
  end subroutine initialize_djch_tlf6
!--------------------------------------------------------------
  function value_at_djch_tlf6(this,t) result(H_out)
    use utils, only:fp,ip
    use linfuncs
    class(djch_tlf6), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp),dimension(this%dim,this%dim) :: H_out
    integer(ip) :: ii
    !H_out = this%value_at_jch(t)


    if(this%use_wref_array) then
       H_out = -(0.0_fp,1.0_fp)*((this%wq_mat - get_diag_matrix(this%wref_array)))
       do ii=1,this%nsg
          H_out = H_out - this%w_mat*(0.0_fp,1.0_fp)*this%signals(ii)%obj%value_at(t,this%wref_array)
       end do

       !write(*,*) 'wq + wtlf'
       !call print_matrix(this%wq_mat)
       !call print_matrix(this%wtlf_mat)
    else
       write(*,*) 'value_at_djch_tlf6 cannot be used with wref (you must provide wref_array)'
       stop(1)
    end if
  end function value_at_djch_tlf6
!--------------------------------------------------------------
  !#######################################################################
  ! Initializes the Hamiltonian of a cooper pair box with asymmetric josephson junctions. Ec = charging energy/hb, wj1 = the energy/hb of the junction with higher energy, d = junction assymetry, phi = magnetic flux/flux quantum, nel = number of energy levels to include (you need enough of them to get correct eigenmodes), ng = gate charge.
  subroutine initialize_cpb(this,wc,wj1,d,phi,nel,ng)
    use utils, only:fp,ip
    use constants, only: hbar
    class(cpb),intent(inout) :: this
    real(fp),intent(in) :: d,wc,wj1,phi
    integer(ip),intent(in),optional :: nel,ng
    real(fp) :: wj2
    integer(ip) :: ii
    if(present(nel)) then
       this%dim = nel
    else
       this%dim = 21
    end if
    if(this%dim < 1) then 
       write(*,*) 'Number of levels must be larger than 0.'
       stop(1)
    end if
    
    if(present(ng)) then
       this%ng = ng
    else
       this%ng = this%dim-1
    end if

    if(d < 0 .or. d>1) then
       write(*,*) 'Error in initializing cpb. Junction assymetry d must be in interval [0,1].'
       stop(1)
    endif
    
    allocate(this%H(this%dim,this%dim))
    this%H = (0.0_fp,0.0_fp)
    
    do ii=1,this%dim
       this%H(ii,ii) = 4*wc*(ii - this%ng/2.0_fp)**2
    end do
    
    wj2 = wj1*(1-d)/(1+d)
    
    do ii=1,this%dim-1
       this%H(ii,ii+1) = (wj1 + wj2)/2.0_fp*cos(phi/2.0_fp) &
            + (wj1 - wj2)*(0.0_fp,1.0_fp)/2.0_fp*sin(phi/2.0_fp)
       this%H(ii+1,ii) = (wj1 + wj2)/2.0_fp*cos(phi/2.0_fp) &
            - (wj1 - wj2)*(0.0_fp,1.0_fp)/2.0_fp*sin(phi/2.0_fp)
    end do
  end subroutine initialize_cpb
!---------------------------------------------------------------
  function value_at_cpb(this,t) result(H_out)
    use utils, only:fp,ip
    class(cpb), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp),dimension(this%dim,this%dim) :: H_out
    H_out = -(0.0_fp,1.0_fp)*this%H
  end function value_at_cpb

  subroutine delete_cpb(this)
    implicit none
    class(cpb), intent(inout) :: this
    deallocate(this%H)
  end subroutine delete_cpb
!#############################################################!
  ! This is not checked and probably useless function. Tries to model two coupled cooper pair boxes having an equation H = [Hcpb_1, delta*sin(phi)+delta0;delta*sin(phi) + delta_0, Hcpb_2].
  subroutine initialize_cpb_tl(this,wc,wj1g,wj1e,d,phi,delta0,delta,nel,ng)
    use utils, only:fp,ip
    use constants, only: hbar
    class(cpb_tl),intent(inout) :: this
    real(fp),intent(in) :: d,wc,wj1g,wj1e,phi,delta0,delta
    integer(ip),intent(in),optional :: nel,ng
    real(fp) :: wj2g,wj2e
    integer(ip) :: ii,jj
    if(present(nel)) then
       this%dim = nel
    else
       this%dim = 21
    end if
    if(nel < 1) then 
       write(*,*) 'Number of levels must be larger than 0.'
       stop(1)
    end if
    
    if(present(ng)) then
       this%ng = ng
    else
       this%ng = this%dim-1
    end if

    if(d < 0 .or. d>1) then
       write(*,*) 'Error in initializing cpb. Junction assymetry d must be in interval [0,1].'
       stop(1)
    endif
    
    allocate(this%H(this%dim*2,this%dim*2))
    this%H = (0.0_fp,0.0_fp)
    
    do ii=1,this%dim
       do jj=0,1
          this%H(ii+this%dim*jj,ii+this%dim*jj) = 4*wc*(ii - this%ng/2.0_fp)**2
       end do
    end do
    
    wj2g = wj1g*(1-d)/(1+d)
    wj2e = wj1e*(1-d)/(1+d)
    
    do ii=1,this%dim-1
       this%H(ii,ii+1) = (wj1g + wj2g)/2.0_fp*cos(phi/2.0_fp) &
            + (wj1g - wj2g)*(0.0_fp,1.0_fp)/2.0_fp*sin(phi/2.0_fp)
       this%H(ii+1,ii) = (wj1g + wj2g)/2.0_fp*cos(phi/2.0_fp) &
            - (wj1g - wj2g)*(0.0_fp,1.0_fp)/2.0_fp*sin(phi/2.0_fp)
    end do
    do ii=1,this%dim-1
       this%H(ii+this%dim,ii+1+this%dim) = (wj1e + wj2e)/2.0_fp*cos(phi/2.0_fp) &
            + (wj1e - wj2e)*(0.0_fp,1.0_fp)/2.0_fp*sin(phi/2.0_fp)
       this%H(ii+1,ii) = (wj1e + wj2e)/2.0_fp*cos(phi/2.0_fp) &
            - (wj1e - wj2e)*(0.0_fp,1.0_fp)/2.0_fp*sin(phi/2.0_fp)
    end do
    
  end subroutine initialize_cpb_tl
!---------------------------------------------------------------
  function value_at_cpb_tl(this,t) result(H_out)
    use utils, only:fp,ip
    class(cpb_tl), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp),dimension(this%dim,this%dim) :: H_out
    H_out = -(0.0_fp,1.0_fp)*this%H
  end function value_at_cpb_tl
!################################################################
  function value_at_disp_jch(this,t) result(H_out)
    use utils, only:fp,ip
    use linfuncs
    class(disp_jch), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp),dimension(this%dim,this%dim) :: H_out
    integer(ip) :: ii
    !H_out = this%value_at_jch(t)
    if(this%use_wref_array) then
       H_out = -(0.0_fp,1.0_fp)*((this%wq_mat - get_diag_matrix(this%wref_array)&
            + this%wc*this%op%nn(:,:,2)) + this%qcca + this%qacc)
       do ii=1,this%nsg
          H_out = H_out - (0.0_fp,1.0_fp)*this%signals(ii)%obj%value_at(t,this%wref_array)
       end do
    else
       H_out = -(0.0_fp,1.0_fp)*((this%wq_mat - this%op%nn(:,:,1)*this%wref &
            + (this%wc-this%wref)*this%op%nn(:,:,2)) + this%qcca + this%qacc)
       do ii=1,this%nsg
          H_out = H_out - (0.0_fp,1.0_fp)*this%signals(ii)%obj%value_at(t,this%wref)
       end do
    end if
  end function value_at_disp_jch
!--------------------------------------------------------------
  subroutine initialize_disp_jch(this,nql_in,ncl_in,wq_in,wc_in,g_in,garray_in,hbar_in,wref_in,wref_array,initial_sg_cont_size)
    use utils, only:fp,ip
    use constants, only: hbar
    use operators
    use linfuncs
    implicit none
    class(disp_jch) :: this
    real(fp),dimension(:),intent(in) :: wq_in
    real(fp),intent(in) :: wc_in
    real(fp),intent(in),optional :: hbar_in,g_in,wref_in
    real(fp),dimension(:),intent(in),optional :: garray_in
    real(fp),dimension(:),intent(in), optional :: wref_array
    integer(ip),intent(in) :: nql_in,ncl_in
    !integer(ip),parameter :: initial_nbr_sg = 5
    integer(ip),intent(in) :: initial_sg_cont_size
    write(*,*) 'Warning NOT DONE YET'
    !call initialize_jch(this,nql_in,ncl_in,wq_in,wc_in,g_in,garray_in,hbar_in)
    allocate(this%signals(initial_sg_cont_size))
    this%sg_cont_size = initial_sg_cont_size
    this%nsg = 0
    if(present(wref_in)) then
       this%wref = wref_in
    else
       this%wref = 0.0_fp
    end if
    
    allocate(this%wref_array(this%dim))
    if(present(wref_array)) then
       this%wref_array = wref_array
       this%use_wref_array = .TRUE.
    else
       this%use_wref_array = .FAlSE.
    endif
  end subroutine initialize_disp_jch
!---------------------------------------------------------------
  subroutine add_signal_disp_jch(this,signal_in)
    use utils,only:fp,ip
    use signals
    class(disp_jch),intent(inout) :: this
    class(signal),intent(in),target :: signal_in
    type(signal_container),allocatable,dimension(:) :: temp_cont
    if(this%nsg >= this%sg_cont_size) then
       allocate(temp_cont(this%nsg))
       temp_cont = this%signals(1:this%nsg)
       deallocate(this%signals)
       allocate(this%signals(this%nsg*2))
       this%signals(1:this%nsg) = temp_cont
       deallocate(temp_cont)
       this%sg_cont_size = size(this%signals)
    end if
    this%nsg = this%nsg + 1
    this%signals(this%nsg)%obj => signal_in
  end subroutine add_signal_disp_jch
!--------------------------------------------------------------------!
  subroutine set_wref_disp_jch(this,wref)
    implicit none
    class(disp_jch),intent(inout) :: this
    real(fp),intent(in) :: wref
    this%wref = wref
    this%use_wref_array = .FALSE.
  end subroutine set_wref_disp_jch
!--------------------------------------------------------------------!
  subroutine set_wref_array_disp_jch(this,wref_array)
    implicit none
    class(disp_jch),intent(inout) :: this
    real(fp),dimension(:),intent(in) :: wref_array
    if(size(wref_array) .ne. this%dim) then
       write(*,*) 'Size of wref_array is incorrect.'
       stop(1)
    endif
    this%wref_array = wref_array
    this%use_wref_array = .TRUE.
  end subroutine set_wref_array_disp_jch
!##############################################################################
  function value_at_driven_hamiltonian(this,t) result(H_out)
    use utils, only:fp,ip
    use linfuncs
    class(driven_hamiltonian), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp),dimension(this%dim,this%dim) :: H_out
    integer(ip) :: ii
    if(this%use_wref_array) then
       do ii=1,this%nsg
          H_out = - (0.0_fp,1.0_fp)*this%signals(ii)%obj%value_at(t,this%wref_array)
       end do
    else
       do ii=1,this%nsg
          H_out = - (0.0_fp,1.0_fp)*this%signals(ii)%obj%value_at(t,this%wref)
       end do
    end if
  end function value_at_driven_hamiltonian
!--------------------------------------------------------------
  subroutine initialize_driven_hamiltonian(this,dim,initial_sg_cont_size,wref_in,wref_array)
    use utils, only:fp,ip
    use constants, only: hbar
    use operators
    use linfuncs
    implicit none
    class(driven_hamiltonian) :: this
    real(fp),intent(in),optional :: wref_in
    real(fp),dimension(:),intent(in), optional :: wref_array
    integer(ip),intent(in) :: initial_sg_cont_size,dim
    this%dim = dim
    allocate(this%signals(initial_sg_cont_size))
    this%sg_cont_size = initial_sg_cont_size
    this%nsg = 0
    if(present(wref_in)) then
       this%wref = wref_in
    else
       this%wref = 0.0_fp
    end if
    
    allocate(this%wref_array(this%dim))
    if(present(wref_array)) then
       this%wref_array = wref_array
       this%use_wref_array = .TRUE.
    else
       this%use_wref_array = .FAlSE.
    endif
  end subroutine initialize_driven_hamiltonian
!---------------------------------------------------------------
  subroutine add_signal_driven_hamiltonian(this,signal_in)
    use utils,only:fp,ip
    use signals
    class(driven_hamiltonian),intent(inout) :: this
    class(signal),intent(in),target :: signal_in
    type(signal_container),allocatable,dimension(:) :: temp_cont
    if(this%nsg >= this%sg_cont_size) then
       allocate(temp_cont(this%nsg))
       temp_cont = this%signals(1:this%nsg)
       deallocate(this%signals)
       allocate(this%signals(this%nsg*2))
       this%signals(1:this%nsg) = temp_cont
       deallocate(temp_cont)
       this%sg_cont_size = size(this%signals)
    end if
    this%nsg = this%nsg + 1
    this%signals(this%nsg)%obj => signal_in
  end subroutine add_signal_driven_hamiltonian
!--------------------------------------------------------------------!
  subroutine set_wref_driven_hamiltonian(this,wref)
    implicit none
    class(driven_hamiltonian),intent(inout) :: this
    real(fp),intent(in) :: wref
    this%wref = wref
    this%use_wref_array = .FALSE.
  end subroutine set_wref_driven_hamiltonian
!--------------------------------------------------------------------!
  subroutine set_wref_array_driven_hamiltonian(this,wref_array)
    implicit none
    class(driven_hamiltonian),intent(inout) :: this
    real(fp),dimension(:),intent(in) :: wref_array
    if(size(wref_array) .ne. this%dim) then
       write(*,*) 'Size of wref_array is incorrect.'
       stop(1)
    endif
    this%wref_array = wref_array
    this%use_wref_array = .TRUE.
  end subroutine set_wref_array_driven_hamiltonian
!#######################################################################
  function value_at_jch_qubit(this,t) result(H_out)
    use utils, only:fp,ip
    use linfuncs
    class(jch_qubit), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp),dimension(this%dim,this%dim) :: H_out
    integer(ip) :: ii
    complex(fp),dimension(2,2) :: sz
    sz = 0.0_fp
    sz(1,1) = -1.0_fp
    sz(2,2) = 1.0_fp

    if(this%use_wref_array) then
       H_out = this%driven_hamiltonian%value_at(t) - (0.0_fp,1.0_fp)*(this%wac/2.0_fp*sz - get_diag_matrix(this%wref_array))
    else
       H_out = this%driven_hamiltonian%value_at(t) - (0.0_fp,1.0_fp)*(this%wac/2.0_fp*sz - this%op%nn(:,:,1)*this%wref)
    end if
  end function value_at_jch_qubit
!--------------------------------------------------------------
  subroutine set_eps_d_jch_qubit(this,eps_d)
    use utils, only:fp,ip
    implicit none
    class(jch_qubit),intent(inout) :: this
    real(fp), intent(in) :: eps_d
    real(fp) :: B,delta,delta_c,ksi
    complex(fp) :: beta

    delta = this%wq-this%wc
    ksi = this%g**2/delta
    delta_c = this%wc-this%wm
    this%alpha_e = (0.0_fp,1.0_fp)*eps_d/(-(0.0_fp,1.0_fp)*(delta_c+ksi) - this%kappa/2.0_fp)
    this%alpha_g = (0.0_fp,1.0_fp)*eps_d/(-(0.0_fp,1.0_fp)*(delta_c-ksi) - this%kappa/2.0_fp)
    beta = this%alpha_e - this%alpha_g
    this%Gamma_d = 2*ksi*aimag(this%alpha_g*conjg(this%alpha_e))
    B = 2*ksi*real(this%alpha_g*conjg(this%alpha_e))
    this%wac = this%wq + B + ksi
    this%Gamma_ci = this%eta*this%kappa*beta*conjg(beta)
    this%eps_d = eps_d
  end subroutine set_eps_d_jch_qubit
!--------------------------------------------------------------
  subroutine initialize_jch_qubit(this,wq,eps_d,wc,wm,kappa,g,eta)
    use utils, only:fp,ip
    use constants, only: hbar
    use operators
    use linfuncs
    implicit none
    class(jch_qubit),intent(inout) :: this
    real(fp),intent(in) :: wq,eps_d,wc,wm,kappa,g,eta

    call this%driven_hamiltonian%initialize(2,2)
    if(allocated(this%op)) then
       deallocate(this%op)
    endif
    allocate(this%op)
    call this%op%initialize([2])

    this%wq = wq
    this%eps_d = eps_d
    this%wc = wc
    this%wm = wm
    this%kappa = kappa
    this%g = g
    this%eta = eta
    call this%set_eps_d(eps_d)
  end subroutine initialize_jch_qubit
!#######################################################################
  function value_at_trq(this,t) result(H_out)
    use utils, only:fp,ip
    use linfuncs
    class(trq), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp),dimension(this%dim,this%dim) :: H_out
    integer(ip) :: ii
    !H_out = this%value_at_jch(t)
    
    H_out = -(0.0_fp,1.0_fp)*(this%op%nn(:,:,1)*this%dq + this%op%nn(:,:,2)*this%da &
         + this%op%nn(:,:,3)*this%db + &
         this%gab*matmul(this%op%aa(:,:,2),this%op%cc(:,:,3))+conjg(this%gab)*matmul(this%op%cc(:,:,2),this%op%aa(:,:,3)) + &
         this%gqa*matmul(this%op%aa(:,:,1),this%op%cc(:,:,2))+conjg(this%gqa)*matmul(this%op%cc(:,:,1),this%op%aa(:,:,2)) + &
         this%gqb*matmul(this%op%aa(:,:,1),this%op%cc(:,:,3))+conjg(this%gqb)*matmul(this%op%cc(:,:,1),this%op%aa(:,:,3)) + &
         this%Aq*this%op%aa(:,:,1) + conjg(this%Aq)*this%op%cc(:,:,1) + &
         this%Aa*this%op%aa(:,:,2) + conjg(this%Aa)*this%op%cc(:,:,2) + &
         this%Ab*this%op%aa(:,:,3) + conjg(this%Ab)*this%op%cc(:,:,3))
    do ii=1,this%nsg
       H_out = H_out - (0.0_fp,1.0_fp)*this%signals(ii)%obj%value_at(t,0.0_fp)
    end do
  end function value_at_trq
!--------------------------------------------------------------
  subroutine initialize_trq(this,na,nb,dq,da,db,gab,gqa,gqb,Aq,Aa,Ab,hbar_in)
    use utils, only:fp,ip
    use constants, only: hbar
    use operators
    use linfuncs
    implicit none
    class(trq) :: this
    real(fp),intent(in) :: dq,da,db
    complex(fp),intent(in) :: gab,gqa,gqb
    complex(fp),intent(in) :: Aq,Aa,Ab
    integer(ip),intent(in) :: na,nb
    real(fp),intent(in),optional :: hbar_in

    integer(ip),parameter :: initial_sg_cont_size = 5
    
    allocate(this%signals(initial_sg_cont_size))
    this%sg_cont_size = initial_sg_cont_size
    this%nsg = 0

    this%na = na
    this%nb = nb
    this%dq = dq
    this%da = da
    this%db = db
    this%gab = gab
    this%gqa = gqa
    this%gqb = gqb
    this%Aq = Aq
    this%Aa = Aa
    this%Ab = Ab
    if(present(hbar_in)) then
       this%hbar = hbar_in
    else
       this%hbar = hbar
    end if
    


    this%dim = 2*this%na*this%nb

    if(allocated(this%op)) then
       deallocate(this%op)
    endif
    allocate(this%op)
    call this%op%initialize([2, this%na,this%nb])

    ! allocate(this%wq_mat(this%dim,this%dim))
    ! this%wq_mat = kron(get_diag_matrix(wq_in),get_eye(this%ncl))
    ! this%wc = wc_in
    ! allocate(this%qcca(this%dim,this%dim))
    ! this%qcca = matmul(kron(get_diag_matrix(this%g),get_eye(this%ncl)),matmul(this%op%cc(:,:,1),this%op%aa(:,:,2)))

        

  end subroutine initialize_trq
!---------------------------------------------------------------
  subroutine add_signal_trq(this,signal_in)
    use utils,only:fp,ip
    use signals
    class(trq),intent(inout) :: this
    class(signal),intent(in),target :: signal_in
    type(signal_container),allocatable,dimension(:) :: temp_cont
    if(this%nsg >= this%sg_cont_size) then
       allocate(temp_cont(this%nsg))
       temp_cont = this%signals(1:this%nsg)
       deallocate(this%signals)
       allocate(this%signals(this%nsg*2))
       this%signals(1:this%nsg) = temp_cont
       deallocate(temp_cont)
       this%sg_cont_size = size(this%signals)
    end if
    this%nsg = this%nsg + 1
    this%signals(this%nsg)%obj => signal_in
  end subroutine add_signal_trq

  subroutine set_parameters_trq(this,parameter_array)
    use utils,only:fp,ip
    class(trq), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    integer(ip) :: ii
    ii = 0
    if(size(parameter_array) > ii) then
       ii = ii+1
       this%dq = real(parameter_array(ii))
    end if
    if(size(parameter_array) > ii) then
       ii = ii+1
       this%da = real(parameter_array(ii))
    end if
    if(size(parameter_array) > ii) then
       ii = ii+1
       this%db = real(parameter_array(ii))
    end if
    if(size(parameter_array) > ii) then
       ii = ii+1
       this%gab = parameter_array(ii)
    end if
    if(size(parameter_array) > ii) then
       ii = ii+1
       this%gqa = parameter_array(ii)
    end if
    if(size(parameter_array) > ii) then
       ii = ii+1
       this%gqb = parameter_array(ii)
    end if
    if(size(parameter_array) > ii) then
       ii = ii+1
       this%Aq = parameter_array(ii)
    end if
    if(size(parameter_array) > ii) then
       ii = ii+1
       this%Aa = parameter_array(ii)
    end if
    if(size(parameter_array) > ii) then
       ii = ii+1
       this%Ab = parameter_array(ii)
    end if

  end subroutine set_parameters_trq

!#######################################################################
end module hamiltonians
