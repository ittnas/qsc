module operators
  use utils, only: fp, ip
  use linfuncs
  implicit none
  private

  public :: operator_cont

  type :: operator_cont
     integer(ip) :: nbr_elements
     integer(ip) :: dim
     integer(ip), allocatable,dimension(:) :: el_dim
     complex(fp), allocatable,dimension(:,:,:) :: cc,aa,nn,Iphi,Qphi,sx,sy,sz
     complex(fp), allocatable,dimension(:,:) :: eye
     contains
       procedure :: initialize
    end type operator_cont
contains
  
  subroutine initialize(this,el_dim_in)
    use linfuncs
    implicit none
    class(operator_cont) :: this
    integer(ip),dimension(:),intent(in) :: el_dim_in
    integer(ip) :: ii,jj,csize
    integer(ip), dimension(size(el_dim_in)+1) :: edw

    !write(*,*) 'Initializing the operators now.'
    this%nbr_elements = size(el_dim_in)
    allocate(this%el_dim(this%nbr_elements))
    this%el_dim = el_dim_in
    !write(*,*) 'el_dim:', this%el_dim
    edw(1) = 1
    edw(2:this%nbr_elements+1) = this%el_dim
    this%dim = product(this%el_dim)
    !write(*,*) 'dim:', this%dim
    !write(*,*) 'nbr_elements:', this%nbr_elements

    allocate(this%cc(this%dim,this%dim,this%nbr_elements))
    allocate(this%aa(this%dim,this%dim,this%nbr_elements))
    allocate(this%nn(this%dim,this%dim,this%nbr_elements))
    allocate(this%Iphi(this%dim,this%dim,this%nbr_elements))
    allocate(this%Qphi(this%dim,this%dim,this%nbr_elements))
    allocate(this%eye(this%dim,this%dim))
    allocate(this%sx(this%dim,this%dim,this%nbr_elements))
    allocate(this%sy(this%dim,this%dim,this%nbr_elements))
    allocate(this%sz(this%dim,this%dim,this%nbr_elements))
    
    !write(*,*) 'Allocated them all.'

    this%cc = (0.0_fp,0.0_fp)
    this%aa = (0.0_fp,0.0_fp)
    this%nn = (0.0_fp,0.0_fp)
    this%Iphi = (0.0_fp,0.0_fp)
    this%Qphi = (0.0_fp,0.0_fp)
    this%eye = get_eye(this%dim)
    this%sx = (0.0_fp,0.0_fp)
    this%sy = (0.0_fp,0.0_fp)

    do ii=1,this%nbr_elements
       csize = 1
       this%cc(:,:,ii) = get_eye(this%dim)
       this%aa(:,:,ii) = get_eye(this%dim)
       do jj=1,this%nbr_elements
          if(jj==ii) then
             this%cc(1:(csize*edw(jj+1)),1:(csize*edw(jj+1)),ii)&
                  =kron(this%cc(1:csize,1:csize,ii),get_creation_op(edw(jj+1)))
             this%aa(1:(csize*edw(jj+1)),1:(csize*edw(jj+1)),ii)&
                  =kron(this%aa(1:csize,1:csize,ii),get_annihilation_op(edw(jj+1)))
             !write(*,*) 'The same!'
          else
             this%cc(1:(csize*edw(jj+1)),1:(csize*edw(jj+1)),ii)&
                  =kron(this%cc(1:csize,1:csize,ii),get_eye(edw(jj+1)))
             this%aa(1:(csize*edw(jj+1)),1:(csize*edw(jj+1)),ii)&
                  =kron(this%aa(1:csize,1:csize,ii),get_eye(edw(jj+1)))
             !write(*,*) 'Different!'
          endif
          !call print_matrix(cc(1:(csize*edw(jj+1)),1:(csize*edw(jj+1)),ii),.true.)
          
          csize = csize*this%el_dim(jj)
       end do
    end do
    
    do ii=1,this%nbr_elements
       this%nn(:,:,ii) = matmul(this%cc(:,:,ii),this%aa(:,:,ii))
       this%Iphi(:,:,ii) = (0.5,0.0)*(this%cc(:,:,ii)*EXP((0.0,1.0)*0.0)&
            +this%aa(:,:,ii)*EXP(-(0.0,1.0)*0.0))
       this%Qphi(:,:,ii) = (0.0,0.5)*(this%cc(:,:,ii)*EXP((0.0,1.0)*0.0)&
            -this%aa(:,:,ii)*EXP(-(0.0,1.0)*0.0))
       !Note the added factor of 2. It seems to produce correct results, but I don't know why.
       ! The correct results to where???? I removed the factor of two, because it definitely should not be there
       this%sx(:,:,ii) = 1.0*(this%cc(:,:,ii) + this%aa(:,:,ii))
       this%sy(:,:,ii) = -(0.0_fp,1.0_fp)*(this%cc(:,:,ii) - this%aa(:,:,ii))
       this%sz(:,:,ii) = -this%nn(:,:,ii)*2+this%eye ! sz now corresponds to the standard definition
    end do
    !write(*,*) 'The final results.'
    !call print_matrix(cc(:,:,1),.true.)
    !call print_matrix(aa(:,:,1),.true.)
  end subroutine initialize

end module operators
