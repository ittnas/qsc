module fclass
  use utils, only: fp,ip
  use hamiltonians
  use sclass
  implicit none
  private
  !public :: fd,fdsti,fdstd,fdmtd,fdlmtd,fisrk15,fsin,ftraceconst,fggm2dens,fsuplin,linearize_superoperator,fdsti_sparse
  public :: fd,fdsti,fdstd,fdmtd,fdlmtd,fisrk15,fsin,ftraceconst,fggm2dens,fsuplin,linearize_superoperator

  ! This is commented because it causes problems with gcc 4.8.4. It is not quaranteed that the assignemt works properly now.
  !interface assignment(=)
  !   module procedure assign_fd
  !end interface assignment(=)

  type :: fd
   contains
     procedure :: value_at => fd_calculate_value_at
     procedure,private ::fd_initialize
     procedure :: assign_fd
     !procedure :: fd_initialize
     !procedure :: destructor => destructor_fd
     !procedure :: initialize => fd_initialize
     generic :: assignment(=) => assign_fd
     generic,public :: initialize => fd_initialize
     procedure :: delete => delete_fd
  end type fd

  type,extends(fd) :: fdsti
     complex(fp),allocatable,dimension(:,:) :: H
   contains
     procedure :: value_at => fdsti_calculate_value_at
     !procedure :: initialize => fdsti_initialize
     procedure,private :: fdsti_initialize
     procedure :: assign_fdsti
     !procedure :: destructor => destructor_fdsti
     procedure :: delete => delete_fdsti
     generic :: assignment(=) => assign_fdsti
     generic, public :: initialize => fdsti_initialize
  end type fdsti

  ! type,extends(fd) :: fdsti_sparse
  !    integer(ip) :: H
  !  contains
  !    procedure :: value_at => fdsti_sparse_calculate_value_at
  !    !procedure :: initialize => fdsti_initialize
  !    procedure,private :: fdsti_sparse_initialize
  !    procedure :: assign_fdsti_sparse
  !    !procedure :: destructor => destructor_fdsti
  !    procedure :: delete => delete_fdsti_sparse
  !    generic :: assignment(=) => assign_fdsti_sparse
  !    generic, public :: initialize => fdsti_sparse_initialize
  ! end type fdsti_sparse

  type,extends(fd) :: fdstd
     !procedure(hamiltonian),pointer,intent(in) :: H
     class(hamiltonian),pointer :: H
   contains
     procedure :: value_at => fdstd_calculate_value_at
     !procedure :: value_at => fdstd_lapack_calculate_value_at
     procedure,private :: fdstd_initialize
     generic,public :: initialize => fdstd_initialize
     procedure :: delete => delete_fdstd
     !final :: destructor_fdstd
  end type fdstd

  type,extends(fd) :: fdmtd
     class(hamiltonian),pointer :: H
   contains
     procedure :: value_at => fdmtd_calculate_value_at
     !procedure :: value_at => fdmtd_lapack_calculate_value_at
     procedure,private :: fdmtd_initialize
     generic,public :: initialize => fdmtd_initialize
     procedure :: delete => delete_fdmtd
     
  end type fdmtd

  type,extends(fd) :: fdlmtd
     class(hamiltonian),pointer :: H
     complex(fp),dimension(:,:,:),allocatable :: lindblad_operators
     complex(fp),dimension(:),allocatable :: lindblad_operator_couplings
     complex(fp),dimension(:,:,:),allocatable :: nd_lindblad_operators_left
     complex(fp),dimension(:,:,:),allocatable :: nd_lindblad_operators_right
     integer(ip) :: lo_cont_size,nlo
     integer(ip) :: nd_cont_size,nds
   contains
     procedure :: value_at => calculate_value_at_fdlmtd
     procedure,private :: initialize_fdlmtd
     generic,public :: initialize => initialize_fdlmtd
     procedure :: add_lindblad_operator => add_lindblad_operator_fdlmtd
     procedure :: add_non_diagonal_lindblad_operator => add_non_diagonal_lindblad_operator_fdlmtd
     procedure :: delete => delete_fdlmtd
     procedure :: set_lindblad_coupling => set_lindblad_coupling_fdlmtd
  end type fdlmtd

  type,extends(fd) :: fisrk15
     class(sd),pointer :: f
     integer(ip) :: nobs
     real(fp) :: x
     complex(fp),allocatable,dimension(:) :: y0,ytilde
     real(fp) :: dx
   contains
     procedure :: value_at => fisrk15_calculate_value_at
     procedure,private :: fisrk15_initialize
     generic,public :: initialize => fisrk15_initialize
     procedure :: delete => delete_fisrk15
  end type fisrk15

  type, extends(fd) :: fsin
     complex(fp) :: A
     complex(fp) :: B
   contains
     procedure :: value_at => fsin_value_at
     procedure,private :: fsin_initialize
     generic,public :: initialize => fsin_initialize
     procedure :: delete => delete_fsin
  end type fsin

  type, extends(fd) :: ftraceconst
     class(fd),pointer :: f
     integer(ip) :: dim
   contains
     procedure :: value_at => ftraceconst_value_at
     procedure,private :: ftraceconst_initialize
     generic,public :: initialize => ftraceconst_initialize
     procedure :: delete => delete_ftraceconst
  end type ftraceconst

  type, extends(fd) :: fggm2dens
     class(fd),pointer :: f
     integer(ip) :: dim
     complex(fp),allocatable,dimension(:) :: rho,rho_out
     complex(fp),allocatable,dimension(:,:,:) :: ggm
     real(fp),allocatable,dimension(:) :: ggmnorm
     real(fp),allocatable,dimension(:) :: b_out
   contains
     procedure :: value_at => fggm2dens_value_at
     procedure,private :: fggm2dens_initialize
     generic,public :: initialize => fggm2dens_initialize
     procedure :: delete => delete_fggm2dens
  end type fggm2dens
  
  type, extends(fd) :: fsuplin
     class(fd),pointer :: f
     complex(fp),allocatable,dimension(:,:) :: L
     !complex(fp),allocatable,dimension(:) :: c
     complex(fp),dimension(:,:,:),pointer :: M
     real(fp),dimension(:),pointer :: Mnorm
     integer(ip) :: dim
   contains
     procedure :: value_at => fsuplin_value_at
     procedure, private :: fsuplin_initialize
     generic,public :: initialize => fsuplin_initialize
     procedure :: delete => delete_fsuplin
  end type fsuplin
contains
  subroutine fd_calculate_value_at(this,x,y,y_out)
    implicit none
    class(fd) :: this
    real(fp),intent(in) :: x
    complex(fp), dimension(:),intent(in),target :: y
    complex(fp), dimension(:),intent(out) :: y_out
    y_out = y
  end subroutine fd_calculate_value_at
!------------------------------------------------------------!
  subroutine delete_fd(this)
    implicit none
    class(fd) :: this
  end subroutine delete_fd
!------------------------------------------------------------!
  subroutine assign_fd(to, from)
    class(fd), intent(out) :: to
    type(fd), intent(in) :: from
  end subroutine assign_fd
  !------------------------------------------------------------!
  subroutine assign(to, from)
    class(fd), intent(out) :: to
    type(fd), intent(in) :: from
  end subroutine assign
!------------------------------------------------------------!
  subroutine fd_initialize(this,n)
    implicit none
    class(fd) :: this
    integer(ip),intent(in) :: n
    integer(ip) :: ii
  end subroutine fd_initialize
!------------------------------------------------------------!
  subroutine fdsti_initialize(this,Hin)
    implicit none
    class(fdsti) :: this
    integer(ip) :: ii,dim
    complex(fp),dimension(:,:) :: Hin
    dim = size(Hin,1)
    if(allocated(this%H)) then
       deallocate(this%H)
    endif
    allocate(this%H(dim,dim))
    this%H = Hin
  end subroutine fdsti_initialize
!------------------------------------------------------------!
  subroutine assign_fdsti(to, from)
    class(fdsti), intent(out) :: to
    type(fdsti), intent(in) :: from
    if(allocated(from%H)) then
       allocate(to%H(size(from%H,1),size(from%H,2)))
       to%H = from%H
    end if
  end subroutine assign_fdsti
!------------------------------------------------------------!
  subroutine delete_fdsti(this)
    implicit none
    class(fdsti) :: this
    deallocate(this%H)
    call this%fd%delete()
  end subroutine delete_fdsti
!------------------------------------------------------------!
  subroutine fdsti_calculate_value_at(this,x,y,y_out)
    implicit none
    class(fdsti) :: this
    real(fp),intent(in) :: x
    complex(fp), dimension(:),intent(in),target :: y
    complex(fp), dimension(:),intent(out) :: y_out
    complex(fp), dimension(size(y),size(y)) :: test
    !y_out = -(0.0_fp,1.0_fp)*matmul(this%H,y)
    y_out = matmul(this%H,y)
  end subroutine fdsti_calculate_value_at
!------------------------------------------------------------!
  subroutine fdstd_initialize(this,Hin)
    use hamiltonians
    implicit none
    class(fdstd) :: this
    class(hamiltonian),target :: Hin
    this%H => Hin
  end subroutine fdstd_initialize
!------------------------------------------------------------!
  subroutine assign_fdstd(to, from)
    class(fdstd), intent(out) :: to
    type(fdstd), intent(in) :: from
    write(*,*) 'Warning: copying of Hamiltonian not impelmented in assign_fdstd'
    !to%H = from%H
  end subroutine assign_fdstd
!------------------------------------------------------------!
  subroutine fdstd_calculate_value_at(this,x,y,y_out)
    use hamiltonians
    implicit none
    class(fdstd) :: this
    real(fp),intent(in) :: x
    complex(fp), dimension(:),intent(in),target :: y
    complex(fp), dimension(:),intent(out) :: y_out
    !complex(fp), dimension(size(y),size(y)) :: temp
    !y_out = -(0.0_fp,1.0_fp)*matmul(this%H%value_at(x),y)
    if(this%H%dim > 5) then
       y_out = matmul(this%H%value_at(x),y)
    else
       call zgemv('n',this%H%dim,this%H%dim,(1.0_fp,0.0_fp),this%H%value_at(x),this%H%dim,y,1,(0.0_fp,0.0_fp),y_out,1_ip)
    endif
    !temp = this%H%value_at(x)
    !y_out = temp(:,1)
  end subroutine fdstd_calculate_value_at

  ! This subroutine is not required anymore. It is included in the fdstd_calculate_value_at.
  ! Remove when you think it is no longer required.
  subroutine fdstd_lapack_calculate_value_at(this,x,y,y_out)
    implicit none
    class(fdstd) :: this
    real(fp),intent(in) :: x
    complex(fp), dimension(:),intent(in),target :: y
    complex(fp), dimension(:),intent(out) :: y_out
    !complex(fp), dimension(size(y),size(y)) :: temp
    !y_out = -(0.0_fp,1.0_fp)*matmul(this%H%value_at(x),y)
    !y_out = matmul(this%H%value_at(x),y)
    call zgemv('n',this%H%dim,this%H%dim,(1.0_fp,0.0_fp),this%H%value_at(x),this%H%dim,y,1,(0.0_fp,0.0_fp),y_out,1_ip)
    !subroutine cgemv (TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
    !temp = this%H%value_at(x)
    !y_out = temp(:,1)
  end subroutine fdstd_lapack_calculate_value_at

!------------------------------------------------------------!
  subroutine delete_fdstd(this)
    implicit none
    class(fdstd) :: this
    write(*,*) 'Warning: Deleting hamiltonian not implemented in delete_fdstd.'
    call this%fd%delete()
  end subroutine delete_fdstd
!############################################################!  
  subroutine fdmtd_initialize(this,Hin)
    use hamiltonians
    implicit none
    class(fdmtd) :: this
    class(hamiltonian),target :: Hin
    this%H => Hin
  end subroutine fdmtd_initialize

  subroutine fdmtd_calculate_value_at(this,x,y,y_out)
    !use linfuncs
    implicit none
    class(fdmtd) :: this
    real(fp),intent(in) :: x
    complex(fp), dimension(:),intent(in),target :: y
    complex(fp), dimension(:),intent(out) :: y_out
    complex(fp), dimension(:,:),pointer :: y_reshaped
    complex(fp), dimension(this%H%dim,this%H%dim) :: y_out_temp
    complex(fp), dimension(this%H%dim,this%H%dim) :: H,lc,rc
    integer(ip) :: dim

    dim = this%H%dim
    H = this%H%value_at(x)
    y_reshaped(1:dim,1:dim) => y
    if(this%H%dim > 5) then
       lc = 0.0_fp
       rc = 0.0_fp
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),H,dim,&
            y_reshaped,dim,(0_fp,0_fp),lc,dim)
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),y_reshaped,dim,&
            H,dim,(0_fp,0_fp),rc,dim)
       y_out = reshape(lc - rc,[size(y)])
    else
       y_out = reshape(matmul(H,y_reshaped)-matmul(y_reshaped,H),[size(y)])
    end if
    !call print_matrix(H,.false.)
    !call print_matrix(y_reshaped,.false.)
    !call print_matrix(matmul(H,y_reshaped)-matmul(y_reshaped,H),.false.)
    !write(*,*) '############################'

  end subroutine fdmtd_calculate_value_at

  subroutine fdmtd_lapack_calculate_value_at(this,x,y,y_out)
    use linfuncs
    implicit none
    class(fdmtd) :: this
    real(fp),intent(in) :: x
    complex(fp), dimension(:),intent(in),target :: y
    complex(fp), dimension(:),intent(out) :: y_out
    complex(fp), dimension(:,:),pointer :: y_reshaped
    complex(fp), dimension(this%H%dim,this%H%dim) :: y_out_temp
    complex(fp), dimension(this%H%dim,this%H%dim) :: H,lc,rc
    H = this%H%value_at(x)
    y_reshaped(1:this%H%dim,1:this%H%dim) => y
    lc = 0.0_fp
    rc = 0.0_fp
    call zgemm('n','n',this%H%dim,this%H%dim,this%H%dim,(1.0_fp,0.0_fp),H,this%H%dim,&
         y_reshaped,this%H%dim,(0_fp,0_fp),lc,this%H%dim)
    call zgemm('n','n',this%H%dim,this%H%dim,this%H%dim,(1.0_fp,0.0_fp),y_reshaped,this%H%dim,&
         H,this%H%dim,(0_fp,0_fp),rc,this%H%dim)
    !call print_matrix(rc,.false.)
    !call print_matrix(lc,.false.)
    !call print_matrix(H,.false.)
    !call print_matrix(y_reshaped,.false.)
    !call print_matrix(lc,.false.)
    !call print_matrix(rc,.false.)
    !call print_matrix(lc - rc,.false.)
    !write(*,*) '############################'

    y_out = reshape(lc - rc,[size(y)])
    !y_out = reshape(matmul(H,y_reshaped)-matmul(y_reshaped,H),[size(y)])
    !SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA,
    !B, LDB, BETA, C, LDC )
  end subroutine fdmtd_lapack_calculate_value_at
!------------------------------------------------------------!
  subroutine delete_fdmtd(this)
    implicit none
    class(fdmtd) :: this
    write(*,*) 'Warning: Deleting hamiltonian not implemented in delete_fdmtd.'
    call this%fd%delete()
  end subroutine delete_fdmtd
!############################################################!
  subroutine calculate_value_at_fdlmtd(this,x,y,y_out)
    use linfuncs
    implicit none
    class(fdlmtd) :: this
    real(fp),intent(in) :: x
    complex(fp), dimension(:),intent(in),target :: y
    complex(fp), dimension(:),intent(out) :: y_out
    complex(fp), dimension(:,:),pointer :: y_reshaped
    complex(fp), dimension(this%H%dim,this%H%dim) :: y_out_temp
    complex(fp), dimension(this%H%dim,this%H%dim) :: H,lcom,rcom,orc,rco,cor,rc,or,c
    integer(ip) :: dim,ii
    dim = this%H%dim
    H = this%H%value_at(x)
    y_reshaped(1:dim,1:dim) => y
    if(this%H%dim > 5) then
       lcom = 0.0_fp
       rcom = 0.0_fp
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),H,dim,&
            y_reshaped,dim,(0_fp,0_fp),lcom,dim)
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),y_reshaped,dim,&
            H,dim,(0_fp,0_fp),rcom,dim)
       y_out = reshape(lcom - rcom,[size(y)])
    else
       y_out = reshape(matmul(H,y_reshaped)-matmul(y_reshaped,H),[size(y)])
    end if
    
    do ii=1,this%nlo
       c = conjg(transpose(this%lindblad_operators(:,:,ii)))
       rc = (0.0_fp,0.0_fp)
       or = (0.0_fp,0.0_fp)
       rco = (0.0_fp,0.0_fp)
       cor = (0.0_fp,0.0_fp)
       orc = (0.0_fp,0.0_fp)
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),y_reshaped,dim,&
            c,dim,(0_fp,0_fp),rc,dim)
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),this%lindblad_operators(:,:,ii),dim,&
            y_reshaped,dim,(0_fp,0_fp),or,dim)
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),this%lindblad_operators(:,:,ii),dim,&
            rc,dim,(0_fp,0_fp),orc,dim)
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),rc,dim,&
            this%lindblad_operators(:,:,ii),dim,(0_fp,0_fp),rco,dim)
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),c,dim,&
            or,dim,(0_fp,0_fp),cor,dim)
       !call print_matrix(c)
       y_out = y_out + this%lindblad_operator_couplings(ii)*reshape(orc - 0.5_fp*cor - 0.5_fp*rco,[size(y)])
    end do

    do ii=1,this%nds
       c = this%nd_lindblad_operators_right(:,:,ii)
       !call print_matrix(c)
       rc = (0.0_fp,0.0_fp)
       or = (0.0_fp,0.0_fp)
       rco = (0.0_fp,0.0_fp)
       cor = (0.0_fp,0.0_fp)
       orc = (0.0_fp,0.0_fp)
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),y_reshaped,dim,&
            c,dim,(0_fp,0_fp),rc,dim)
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),this%nd_lindblad_operators_left(:,:,ii),dim,&
            y_reshaped,dim,(0_fp,0_fp),or,dim)
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),this%nd_lindblad_operators_left(:,:,ii),dim,&
            rc,dim,(0_fp,0_fp),orc,dim)
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),rc,dim,&
            this%nd_lindblad_operators_left(:,:,ii),dim,(0_fp,0_fp),rco,dim)
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),c,dim,&
            or,dim,(0_fp,0_fp),cor,dim)
       !call print_matrix(orc - 0.5_fp*cor - 0.5_fp*rco)
       y_out = y_out + reshape(orc - 0.5_fp*cor - 0.5_fp*rco,[size(y)])
    end do
  end subroutine calculate_value_at_fdlmtd
!------------------------------------------------------------!
  subroutine initialize_fdlmtd(this,Hin)
    use hamiltonians
    implicit none
    class(fdlmtd) :: this
    class(hamiltonian),target :: Hin
    integer(ip),parameter :: initial_lo_cont_size = 2
    this%H => Hin
    allocate(this%lindblad_operators(this%H%dim,this%H%dim,initial_lo_cont_size))
    allocate(this%lindblad_operator_couplings(initial_lo_cont_size))
    allocate(this%nd_lindblad_operators_left(this%H%dim,this%H%dim,initial_lo_cont_size))
    allocate(this%nd_lindblad_operators_right(this%H%dim,this%H%dim,initial_lo_cont_size))
    this%nlo = 0
    this%lo_cont_size = initial_lo_cont_size
    this%nds = 0
    this%nd_cont_size = initial_lo_cont_size
  end subroutine initialize_fdlmtd
!------------------------------------------------------------!
  subroutine set_lindblad_coupling_fdlmtd(this,ind,coupling)
    implicit none
    class(fdlmtd),intent(inout) :: this
    complex(fp),intent(in) :: coupling
    integer(ip),intent(in) :: ind
    if(ind > this%nlo) then
       write(*,*) 'ind > nlo in set_lindblad_coupling_fdlmtd'
    else
       !this%lindblad_operators(:,:,ind) = this%lindblad_operators(:,:,ind)*sqrt(coupling)/&
       !     sqrt(this%lindblad_operator_couplings(ind))
    end if
    this%lindblad_operator_couplings(ind) = coupling
  end subroutine set_lindblad_coupling_fdlmtd
!------------------------------------------------------------!
  subroutine add_lindblad_operator_fdlmtd(this,op,coupling)
    implicit none
    class(fdlmtd),intent(inout) :: this
    complex(fp),intent(in) :: coupling
    complex(fp),allocatable,dimension(:,:,:) :: temp_cont
    complex(fp),allocatable,dimension(:) :: temp_cont_couplings
    complex(fp),dimension(this%H%dim,this%H%dim),intent(in) :: op

!    if(size(op,1) .ne. this%H%dim .or. size(op,2) .ne. size(op,1)) then
!       write(*,*) 'Trying to add Lindblad operator with incorrect size.'
!       stop(1)
!    endif
    
    if(this%nlo>=this%lo_cont_size) then
       allocate(temp_cont(this%H%dim,this%H%dim,this%nlo))
       allocate(temp_cont_couplings(this%nlo))
       temp_cont = this%lindblad_operators(:,:,1:this%nlo)
       temp_cont_couplings = this%lindblad_operator_couplings(1:this%nlo)
       deallocate(this%lindblad_operators)
       deallocate(this%lindblad_operator_couplings)
       allocate(this%lindblad_operators(this%H%dim,this%H%dim,this%nlo*2))
       allocate(this%lindblad_operator_couplings(this%nlo*2))
       this%lindblad_operators(:,:,1:this%nlo) = temp_cont
       this%lindblad_operator_couplings(1:this%nlo) = temp_cont_couplings
       deallocate(temp_cont)
       deallocate(temp_cont_couplings)
       this%lo_cont_size = size(this%lindblad_operators,3)
    end if
    this%nlo = this%nlo+1
    !this%lindblad_operators(:,:,this%nlo) = sqrt(coupling)*op
    this%lindblad_operators(:,:,this%nlo) = op
    this%lindblad_operator_couplings(this%nlo) = coupling
  end subroutine add_lindblad_operator_fdlmtd
!------------------------------------------------------------!
  ! Non-diagonal operators store coupling and the operator to the same matrix. This should be changed.
  subroutine add_non_diagonal_lindblad_operator_fdlmtd(this,left_op,right_op,coupling)
    implicit none
    class(fdlmtd),intent(inout) :: this
    complex(fp),allocatable,dimension(:,:,:) :: temp_cont_left,temp_cont_right
    complex(fp),dimension(this%H%dim,this%H%dim),intent(in) :: left_op,right_op
    complex(fp),intent(in) :: coupling

!    if(size(op,1) .ne. this%H%dim .or. size(op,2) .ne. size(op,1)) then
!       write(*,*) 'Trying to add Lindblad operator with incorrect size.'
!       stop(1)
!    endif
    
    if(this%nds>=this%nd_cont_size) then
       allocate(temp_cont_left(this%H%dim,this%H%dim,this%nds))
       allocate(temp_cont_right(this%H%dim,this%H%dim,this%nds))
       temp_cont_right = this%nd_lindblad_operators_right(:,:,1:this%nds)
       temp_cont_left = this%nd_lindblad_operators_left(:,:,1:this%nds)
       deallocate(this%nd_lindblad_operators_left)
       deallocate(this%nd_lindblad_operators_right)
       allocate(this%nd_lindblad_operators_left(this%H%dim,this%H%dim,this%nds*2))
       allocate(this%nd_lindblad_operators_left(this%H%dim,this%H%dim,this%nds*2))
       this%nd_lindblad_operators_left(:,:,1:this%nds) = temp_cont_left
       this%nd_lindblad_operators_right(:,:,1:this%nds) = temp_cont_right
       deallocate(temp_cont_left)
       deallocate(temp_cont_right)
       this%nd_cont_size = size(this%nd_lindblad_operators_left,3)
    end if
    this%nds = this%nds+1
    this%nd_lindblad_operators_left(:,:,this%nds) = sqrt(coupling)*left_op
    this%nd_lindblad_operators_right(:,:,this%nds) = sqrt(coupling)*right_op
  end subroutine add_non_diagonal_lindblad_operator_fdlmtd
!------------------------------------------------------------!
  subroutine delete_fdlmtd(this)
    implicit none
    class(fdlmtd) :: this
    write(*,*) 'Warning: Deleting hamiltonian not implemented in delete_fdlmtd.'
    deallocate(this%lindblad_operators)
    deallocate(this%nd_lindblad_operators_left)
    deallocate(this%nd_lindblad_operators_right)
    call this%fd%delete()
  end subroutine delete_fdlmtd
!############################################################!
  subroutine fisrk15_initialize(this,x,y0,ng,f,dx,nobs,j_avg,dW_out)
    use hamiltonians
    use sclass
    implicit none
    class(fisrk15),intent(inout) :: this
    class(sd),target,intent(in) :: f
    integer(ip),intent(in) :: nobs
    real(fp),intent(in) :: x
    complex(fp),dimension(:),intent(in) :: y0
    real(fp),dimension(:,:),intent(in) :: ng
    real(fp),intent(in) :: dx
    real(fp),dimension(nobs),intent(out),optional :: j_avg
    real(fp),dimension(nobs),intent(out),optional :: dW_out

    complex(fp),dimension(size(y0)) :: gp,gm,pp,pm,a,agp,agm,ad
    complex(fp),dimension(size(y0),nobs) :: b,bgp,bgm,bpp,bpm
    real(fp),dimension(nobs) :: dW,dZ
    integer(ip) :: jj
    real(fp) :: sqrtdx
    integer(ip) :: c1,c2
    
    if(.not. allocated(this%ytilde)) then
       allocate(this%y0(size(y0)))
       !allocate(this%ng(size(ng,1),size(ng,2)))
       allocate(this%ytilde(size(y0)))
    endif

    this%f => f
    this%nobs = nobs
    this%x = x
    this%y0 = y0
    !call system_clock(c1)
    !this%ng = ng
    !call system_clock(c2)
    !write(*,*) 'Allocation: ', c2-c1
    !write(*,*) size(this%ng),size(ng)
    this%dx = dx
    
    sqrtdx = sqrt(dx)

    !call system_clock(c1)
    if(present(j_avg)) then
       call f%value_at(x,this%y0,a,b,j_avg)
    else
       call f%value_at(x,this%y0,a,b)
    end if
    gp = this%y0+a*dx
    gm = this%y0+a*dx

    do jj=1,nobs
       gp = gp + sqrtdx*b(:,jj)
       gm = gm - sqrtdx*b(:,jj)
    end do

    call f%value_at(x,gp,agp,bgp)
    call f%value_at(x,gm,agm,bgm)

    pp = gp
    pm = gm

    do jj=1,nobs
       pp = pp + bgp(:,jj)*sqrtdx
       pm = pm - bgp(:,jj)*sqrtdx
    end do

    call f%value_at(x,pp,ad,bpp)
    call f%value_at(x,pm,ad,bpm)
    !call system_clock(c2)
    !write(*,*) 'precalc: ', c2-c1

    !call system_clock(c1)
    ! The stochastic part
    dW = ng(:,1)*sqrtdx
    dZ = 0.5_fp*dx**(1.5_fp)*(ng(:,1) + 1.0_fp/sqrt(3.0_fp)*ng(:,2))
    if(present(dW_out)) then
       dW_out = dW
    end if

    this%ytilde = this%y0 + 0.5_fp*dx*a
    do jj=1,nobs
       this%ytilde = this%ytilde + b(:,jj)*dW(jj) &
            + 0.25_fp/sqrtdx*(bgp(:,jj) - bgm(:,jj)) &
            * (dw(jj)**2 - dx) &
            + 0.5_fp/dx*(bgp(:,jj)-2*b(:,jj) + bgm(:,jj))&
            * (dw(jj)*dx - dZ(jj))&
            + 0.5_fp/sqrtdx*(agp - agm)*(dZ(jj) - 0.5_fp*dW(jj)*dx)&
            + 0.25_fp/dx*(bpp(:,jj) - bpm(:,jj) - bgp(:,jj) + bgm(:,jj))&
            * (1.0_fp/3.0_fp*dW(jj)**2 - dx)*dW(jj)
    end do
    !call system_clock(c2)
    !write(*,*) 'ytilde: ', c2-c1

    !write(*,*) 'this%ytilde',this%ytilde
    !write(*,*) 'b',b
    !write(*,*) 'bgp',bgp
    !write(*,*) 'bgm',bgm
    !write(*,*) 'agp',agp
    !write(*,*) 'agm',agm
    !write(*,*) 'bpp',bpp
    !write(*,*) 'bpm',bpm
    
  end subroutine fisrk15_initialize
  !------------------------------------------------------------!
  subroutine fisrk15_calculate_value_at(this,x,y,y_out)
    implicit none
    class(fisrk15) :: this
    real(fp),intent(in) :: x
    complex(fp), dimension(:),intent(in),target :: y
    complex(fp), dimension(:),intent(out) :: y_out
    complex(fp), dimension(size(y),this%nobs) :: bd
    call this%f%value_at(this%x,y,y_out,bd) ! Depends only at time at x0
    y_out = (0.5_fp*this%dx)*y_out+this%ytilde
  end subroutine fisrk15_calculate_value_at
  !------------------------------------------------------------!
  subroutine delete_fisrk15(this)
    implicit none
    class(fisrk15) :: this
    if(allocated(this%ytilde)) then ! In practice it should always be allocated
       deallocate(this%y0,this%ytilde)
    endif
    call this%fd%delete()
  end subroutine delete_fisrk15
!############################################################!
  subroutine fsin_value_at(this,x,y,y_out)
    implicit none
    class(fsin) :: this
    real(fp),intent(in) :: x
    complex(fp), dimension(:),intent(in),target :: y
    complex(fp), dimension(:),intent(out) :: y_out
    y_out = this%A*sin(this%B*y + x)
  end subroutine fsin_value_at
  !------------------------------------------------------------!
  subroutine fsin_initialize(this,A,B)
    implicit none
    class(fsin) :: this
    complex(fp),intent(in) :: A
    complex(fp),intent(in) :: B
    this%A = A
    this%B = B
  end subroutine fsin_initialize
!------------------------------------------------------------!
  subroutine delete_fsin(this)
    implicit none
    class(fsin) :: this
    call this%fd%delete()
  end subroutine delete_fsin
!############################################################!
  subroutine ftraceconst_initialize(this,f,dim)
    implicit none
    class(ftraceconst) :: this
    class(fd),target,intent(in) :: f
    integer(ip),intent(in) :: dim
    this%f => f
    this%dim = dim
  end subroutine ftraceconst_initialize
  !------------------------------------------------------------!
  subroutine ftraceconst_value_at(this,x,y,y_out)
    implicit none
    class(ftraceconst) :: this
    real(fp),intent(in) :: x
    complex(fp), dimension(:),intent(in),target :: y
    complex(fp), dimension(:),intent(out) :: y_out
    !complex(fp),dimension(this%dim) :: yf
    integer(ip) :: ii
    call this%f%value_at(x,y(1:(this%dim**2)),&
         y_out(1:(this%dim**2)))
    !y_out(this%dim**2+1) = -(1.0_fp,0.0_fp)
    y_out(this%dim**2+1) = -(1.0_fp,0.0_fp)!+y(this%dim**2+1)
    do ii=1,this%dim
       y_out(this%dim**2+1) = y_out(this%dim**2+1) + y(this%dim*(ii-1)+ii)
    end do
    y_out(1:this%dim**2) = y_out(1:this%dim**2) + y(this%dim**2+1)
    !write(*,*) 'y:',y(this%dim**2+1)
  end subroutine ftraceconst_value_at
!------------------------------------------------------------!
  subroutine delete_ftraceconst(this)
    implicit none
    class(ftraceconst) :: this
    call this%fd%delete()
  end subroutine delete_ftraceconst
!############################################################!
  subroutine fggm2dens_initialize(this,f,dim)
    use linfuncs
    implicit none
    class(fggm2dens) :: this
    class(fd),target,intent(in) :: f
    integer(ip),intent(in) :: dim
    integer(ip) :: c1,c2
    this%f => f
    this%dim = dim
    allocate(this%rho(dim**2),this%rho_out(dim**2))
    allocate(this%ggm(dim,dim,dim**2))
    allocate(this%b_out(dim**2))
    allocate(this%ggmnorm(dim**2))
    !call system_clock(c1)
    call generate_ggm(dim,this%ggm,this%ggmnorm)
    !call system_clock(c2)
    !write(*,*) 'Generate_ggm. Duration: ', c2-c1,'ms.'
  end subroutine fggm2dens_initialize
  !------------------------------------------------------------!
  subroutine fggm2dens_value_at(this,x,y,y_out)
    use linfuncs
    implicit none
    class(fggm2dens) :: this
    real(fp),intent(in) :: x
    complex(fp), dimension(:),intent(in),target :: y
    complex(fp), dimension(:),intent(out) :: y_out
    call ggm_to_density_matrix(this%rho,this%ggm,this%ggmnorm,real(y),this%dim)
    !write(*,*) 'rho:'
    !call print_matrix(reshape(this%rho,[this%dim,this%dim]))
    call this%f%value_at(x,this%rho,this%rho_out)
    call density_matrix_to_ggm(reshape(this%rho_out,[this%dim,this%dim]),&
         this%ggm,this%b_out)
    y_out = this%b_out
  end subroutine fggm2dens_value_at
!---------------------------------------------------------------!
  subroutine delete_fggm2dens(this)
    implicit none
    class(fggm2dens) :: this
    deallocate(this%b_out,this%ggm,this%rho,this%rho_out,this%ggmnorm)
    call this%fd%delete()
  end subroutine delete_fggm2dens
!############################################################!
  subroutine fsuplin_initialize(this,f,M,Mnorm,dim)
    use linfuncs
    implicit none
    class(fsuplin) :: this
    class(fd),intent(in),target :: f
    integer(ip),intent(in) :: dim
    complex(fp),dimension(:,:),pointer :: Lrho,Li
    complex(fp),dimension(dim,dim) :: Lrho_out
    complex(fp),dimension(dim**2),target :: Lrho_lin,Li_lin
    complex(fp),intent(in),dimension(dim,dim,dim**2),target :: M
    real(fp),intent(in),dimension(dim**2),target :: Mnorm
    complex(fp),dimension(dim**2) :: M_rs
    integer(ip) :: ii,jj
    this%M => M
    this%Mnorm => Mnorm
    this%dim = dim
    this%f => f
    if(.not. allocated(this%L)) then
       allocate(this%L(dim**2,dim**2))
    end if
    call linearize_superoperator(f,M,Mnorm,dim,this%L)
  end subroutine fsuplin_initialize
!---------------------------------------------------------------!
  ! this is not used at the moment
  subroutine fsuplin_reinitialize(this)
    implicit none
    class(fsuplin) :: this
    call linearize_superoperator(this%f,this%M,this%Mnorm,this%dim,this%L)
  end subroutine fsuplin_reinitialize
!---------------------------------------------------------------!
  subroutine fsuplin_value_at(this,x,y,y_out)
    use linfuncs
    implicit none
    class(fsuplin) :: this
    real(fp),intent(in) :: x
    complex(fp), dimension(:),intent(in),target :: y
    complex(fp), dimension(:),intent(out) :: y_out
    integer(ip) :: dim
    dim = size(y)
    !y_out = this%c
    call zgemv('n',dim,dim,(1.0_fp,0.0_fp),this%L,dim,&
         y,1,(0.0_fp,0.0_fp),y_out,1)
    !write(*,*) 'y_out before:',real(y_out)
    !y_out = y_out + this%c
    !write(*,*) 'y_out after:',real(y_out)
    !y_out = matmul(this%L,y)
    !write(*,*) y_out
  end subroutine fsuplin_value_at
!---------------------------------------------------------------!
  subroutine delete_fsuplin(this)
    implicit none
    class(fsuplin) :: this
    deallocate(this%L)
  end subroutine delete_fsuplin
  !---------------------------------------------------------------!
  ! Creates a linear presentation of a time independent superoperator. Linearized in the basis given by M. Lmn = (Mm,L(Mn)) = Tr(Mm^d*L(Mn)), from http://physics.stackexchange.com/questions/163546/finding-the-matrix-representation-of-a-superoperator
  subroutine linearize_superoperator(f,M,Mnorm,dim,L)
    use linfuncs
    implicit none
    class(fd),intent(in),target :: f
    complex(fp),intent(in),dimension(dim,dim,dim**2) :: M
    real(fp),intent(in),dimension(dim**2) :: Mnorm
    integer(ip),intent(in) :: dim
     complex(fp),dimension(dim**2,dim**2),intent(out) :: L
    complex(fp),dimension(:,:),pointer :: Lrho
    complex(fp),dimension(dim**2),target :: Lrho_lin
    complex(fp),dimension(dim,dim) :: Lrho_out
    integer(ip) :: ii,jj,c1,c2
    
    Lrho(1:dim,1:dim) => Lrho_lin
    L = 0.0_fp

    do ii=1,dim**2
       call f%value_at(0.0_fp,reshape(M(:,:,ii),[dim**2]),Lrho_lin)
       !call system_clock(c1)
       do jj=1,dim**2
          call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),Lrho,dim,&
               conjg(transpose(M(:,:,jj))),dim,(0.0_fp,0.0_fp),Lrho_out,dim)
          L(jj,ii) = trace(Lrho_out)/Mnorm(jj)
       end do
       !call system_clock(c2)
       !write(*,*) 'One round',c2-c1,'ms.'
    end do
    L(:,1) = L(:,1)*2.0_fp/dim ! MAGIC
    write(*,*) 'Warnining, this is a weird hack linearize_superoperator'
  end subroutine linearize_superoperator
!#########################################################!
!   subroutine fdsti_sparse_initialize(this,Hin)
!     use blas_sparse
!     implicit none
!     class(fdsti_sparse) :: this
!     integer(ip) :: ii,dim,jj,istat,H,property
!     complex(fp),dimension(:,:) :: Hin
!     dim = size(Hin,1)
!     !if(allocated(this%H)) then
!     !   deallocate(this%H)
!     !endif
!     property= blas_general + blas_col_major
!     call zuscr_begin(dim,dim,H,istat)
!     call ussp(H,property,istat)
!     write(*,*) 'Initializing sparse matrix.'
!     do ii=1,dim
!        do jj=1,dim
!           if(abs(Hin(ii,jj)) > 0.0000001_fp) then
!              call uscr_insert_entry(H,Hin(ii,jj),ii,jj,istat)
!              write(*,*) 'ii=',ii,', jj=',jj,', H(ii,jj)=',Hin(ii,jj)
!           end if
!        end do
!     end do
!     call uscr_end(H,istat)
!     this%H = H
!   end subroutine fdsti_sparse_initialize
! !------------------------------------------------------------!
!   subroutine assign_fdsti_sparse(to, from)
!     class(fdsti_sparse), intent(out) :: to
!     type(fdsti_sparse), intent(in) :: from
!     write(*,*) 'Not implemented!'
!     call exit(1)
!   end subroutine assign_fdsti_sparse
! !------------------------------------------------------------!
!   subroutine delete_fdsti_sparse(this)
!     use blas_sparse
!     implicit none
!     class(fdsti_sparse) :: this
!     integer(ip) :: istat
!     call usds(this%H,istat)
!     call this%fd%delete()
!   end subroutine delete_fdsti_sparse
! !------------------------------------------------------------!
!   subroutine fdsti_sparse_calculate_value_at(this,x,y,y_out)
!     use blas_sparse
!     implicit none
!     class(fdsti_sparse) :: this
!     real(fp),intent(in) :: x
!     complex(fp), dimension(:),intent(in),target :: y
!     complex(fp), dimension(:),intent(out) :: y_out
!     complex(fp), dimension(size(y),size(y)) :: test
!     integer(ip) :: istat
!     !y_out = -(0.0_fp,1.0_fp)*matmul(this%H,y)
!     !y_out = matmul(this%H,y)
!     call usmv(this%H,y,y_out,istat)
!   end subroutine fdsti_sparse_calculate_value_at
! !------------------------------------------------------------!
end module fclass
