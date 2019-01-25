module sclass
  use utils, only: fp, ip
  use hamiltonians
  implicit none
  private
  public :: sd,diff_schr,diff_sme
  
  type :: sd
     complex(fp),dimension(:,:,:),allocatable :: observables
     complex(fp),dimension(:),allocatable :: couplings
     integer(ip) :: nobs,obs_cont_size
   contains
     !procedure(value_at_interface),deferred :: value_at
     procedure :: value_at => value_at_sd
     procedure :: normalize => normalize_sd
     procedure :: add_observable => add_observable_sd
     procedure :: set_observable_coupling => set_observable_coupling_sd
  end type sd
  
  abstract interface
     subroutine value_at_interface(this,x,y,d,s,j_avg)
       use utils,only:ip,fp
       import sd
       implicit none
       class(sd),intent(in) :: this
       real(fp),intent(in) :: x
       complex(fp),dimension(:),intent(in),target :: y
       complex(fp), dimension(size(y)),intent(out) :: d
       complex(fp), dimension(size(y),this%nobs),intent(out) :: s
       real(fp),dimension(this%nobs),intent(out),optional :: j_avg
     end subroutine value_at_interface
  end interface

  type,extends(sd) :: diff_schr
     class(hamiltonian),pointer :: H
   contains
     procedure :: value_at => value_at_diff_schr
     procedure :: initialize => init_diff_schr
     procedure :: normalize => normalize_diff_schr
  end type diff_schr

  type,extends(sd) :: diff_sme
     class(hamiltonian),pointer :: H
     complex(fp),dimension(:,:,:),allocatable :: lindblad_operators
     complex(fp),dimension(:),allocatable :: lindblad_operator_couplings
     integer(ip) :: lo_cont_size,nlo
     complex(fp),dimension(:,:),allocatable :: lo_comb,loclo_comb
   contains
     procedure :: value_at => value_at_diff_sme
     procedure :: initialize => init_diff_sme
     procedure :: normalize => normalize_diff_sme
     procedure :: add_lindblad_operator => add_lindblad_operator_diff_sme
     procedure :: add_observable => add_observable_diff_sme
     procedure :: set_lindblad_coupling => set_lindblad_coupling_diff_sme

  end type diff_sme
contains
!##############################################################
  ! subroutine init_sd(this)
  !   use utils,only:fp,ip
  !   use hamiltonians
  !   implicit none
  !   class(diff_schr),intent(inout) :: this
  !   this%nobs = 0
  !   this%obs_cont_size = 5
  ! end subroutine init_sd

  subroutine value_at_sd(this,x,y,d,s,j_avg)
    use utils,only:ip,fp
    use linfuncs
    implicit none
    class(sd),intent(in) :: this
    real(fp),intent(in) :: x
    complex(fp),dimension(:),intent(in),target :: y
    complex(fp), dimension(size(y)),intent(out) :: d
    complex(fp), dimension(size(y),this%nobs),intent(out) :: s
    real(fp),dimension(this%nobs),intent(out),optional :: j_avg
    
  end subroutine value_at_sd


  subroutine add_observable_sd(this,A,coupling,efficiency_in)
    implicit none
    class(sd) :: this
    complex(fp),dimension(:,:),intent(in) :: A
    complex(fp),intent(in) :: coupling
    integer(ip),parameter :: obs_cont_initial_size = 5
    complex(fp),allocatable,dimension(:,:,:) :: observables_t
    complex(fp),allocatable,dimension(:) :: couplings_t
    real(fp),intent(in),optional :: efficiency_in
    real(fp) :: efficiency

    if(present(efficiency_in)) then
       efficiency = efficiency_in
    else
       efficiency = 1.0_fp
    end if
    !integer(ip) :: ii
    if(allocated(this%observables) .eqv. .false.) then
       allocate(this%observables(size(A,1),size(A,2),obs_cont_initial_size))
       allocate(this%couplings(obs_cont_initial_size))
       this%nobs = 0
       this%obs_cont_size = obs_cont_initial_size
    endif
    
    if(this%nobs >= this%obs_cont_size) then !make more space
       allocate(observables_t(size(A,1),size(A,2),this%nobs))
       allocate(couplings_t(this%nobs))

       observables_t = this%observables
       couplings_t = this%couplings
       deallocate(this%observables)
       deallocate(this%couplings)
       this%obs_cont_size = this%obs_cont_size*2
       allocate(this%observables(size(A,1),size(A,2),this%obs_cont_size))
       allocate(this%couplings(this%obs_cont_size))
       this%observables(:,:,1:this%nobs) = observables_t
       this%couplings(1:this%nobs) = couplings_t
       deallocate(observables_t)
       deallocate(couplings_t)
    endif
    this%nobs = this%nobs + 1
    this%observables(:,:,this%nobs) = A
    this%couplings(this%nobs) = coupling*efficiency
  end subroutine add_observable_sd

  subroutine set_observable_coupling_sd(this,coupling,n,efficiency)
    implicit none
    class(sd),intent(inout) :: this
    complex(fp),intent(in) :: coupling
    integer(ip),intent(in) :: n
    real(fp),intent(in),optional :: efficiency

    if(this%nobs < n) then
       write(*,*) 'Trying to set the coupling of nonexisting operator.'
       return
    end if
    !write(*,*) 'Setting observable coupling'
    if(present(efficiency)) then
       this%couplings(n) = coupling*efficiency
    else 
       this%couplings(n) = coupling*1.0_fp
    end if

  end subroutine set_observable_coupling_sd

  subroutine normalize_sd(this,y)
    use utils,only:fp,ip
    implicit none
    class(sd),intent(in) :: this
    complex(fp),dimension(:),intent(inout),target :: y
  end subroutine normalize_sd
!##############################################################
  subroutine init_diff_schr(this,H)
    use utils,only:fp,ip
    use hamiltonians
    implicit none
    class(diff_schr),intent(inout) :: this
    class(hamiltonian), target :: H
    this%H => H
  end subroutine init_diff_schr

  subroutine value_at_diff_schr(this,x,y,d,s,j_avg)
    use utils,only:ip,fp
    use linfuncs
    implicit none
    class(diff_schr),intent(in) :: this
    real(fp),intent(in) :: x
    complex(fp),dimension(:),intent(in),target :: y
    complex(fp), dimension(size(y)),intent(out) :: d
    complex(fp), dimension(size(y),this%nobs),intent(out) :: s
    real(fp),dimension(this%nobs),intent(out),optional :: j_avg
    complex(fp),dimension(this%nobs) :: avg_a
    integer(ip) :: ii
    complex(fp),dimension(this%H%dim,this%H%dim,this%nobs) :: Act
    complex(fp),dimension(this%H%dim,this%H%dim) :: temp
    complex(fp),dimension(size(y)) :: y_temp

    ! Can be found from Breur-Petruccione pp. 341

    ! if(allocated(this%observables) .eqv. .false.) then
    !    this%nobs = 0
    ! endif
    
    !d = (0.0_fp,0.0_fp)
    !s = (0.0_fp,0.0_fp)

    do ii=1,this%nobs
       Act(:,:,ii) = conjg(transpose(this%observables(:,:,ii)))
       !avg_a(ii) = dot_product(y,matmul(this%observables(:,:,ii) + Act(:,:,ii),y))
       call zgemv('n',this%H%dim,this%H%dim,(1.0_fp,0.0_fp),this%observables(:,:,ii) + Act(:,:,ii),&
            this%H%dim,y,1,(0.0_fp,0.0_fp),y_temp,1)
       avg_a(ii) = dot_product(y,y_temp)
       !write(*,*) 'avg_a ',ii,':',avg_a(ii)
    end do
    temp = this%H%value_at(x)
    !write(*,*) 'H:'
    !call print_matrix(temp)
    ! (A + a)|psi> = A|psi> + a|psi> = (A + aI)|psi>
    do ii=1,this%nobs
       !temp = temp + this%couplings(ii)*0.5_fp*avg_a(ii)*(this%observables(:,:,ii) - 0.25_fp*avg_a(ii))
       temp = temp + this%couplings(ii)*0.5_fp*avg_a(ii)*(this%observables(:,:,ii))
    end do
    !d = matmul(temp,y)
    call zgemv('n',this%H%dim,this%H%dim,(1.0_fp,0.0_fp),temp,this%H%dim&
         ,y,1,(0.0_fp,0.0_fp),d,1)

    do ii=1,this%nobs
       d = d - this%couplings(ii)*0.125_fp*avg_a(ii)**2*y
    end do

    !write(*,*) 'temp:'
    !call print_matrix(temp)
    !write(*,*) 'avg_a', avg_a
    do ii=1,this%nobs
       !d = d - this%couplings(ii)*0.5_fp*matmul(Act(:,:,ii),matmul(this%observables(:,:,ii),y))
       call zgemv('n',this%H%dim,this%H%dim,(1.0_fp,0.0_fp),this%observables(:,:,ii),this%H%dim,y,1,(0.0_fp,0.0_fp),y_temp,1)
       call zgemv('n',this%H%dim,this%H%dim,this%couplings(ii)*0.5_fp,Act(:,:,ii),this%H%dim,y_temp,1,(0.0_fp,0.0_fp),y_temp,1)
       d = d - y_temp
    end do
    !s = (0.0_fp,0.0_fp)
    do ii=1,this%nobs
       !s(:,ii) = sqrt(this%couplings(ii))*(matmul(this%observables(:,:,ii),y) - 0.5_fp*avg_a(ii)*y)
       s(:,ii) = y
       call zgemv('n',this%H%dim,this%H%dim,sqrt(this%couplings(ii)),&
            this%observables(:,:,ii),this%H%dim,y,1,-sqrt(this%couplings(ii))*0.5_fp*avg_a(ii),s(:,ii),1)
       !s(:,ii) = y     ! y CANNOT BE USED ANYMORE, ITS MODIFIED IN ZGEMV!
    end do
    
    ! if(this%H%dim > 5) then
    !    d = matmul(this%H%value_at(x),y)
    ! else
    !    call zgemv('n',this%H%dim,this%H%dim,(1.0_fp,0.0_fp),this%H%value_at(x),this%H%dim,y,1,(0.0_fp,0.0_fp),d,1_ip)
    ! endif
    if(present(j_avg)) then
       write(*,*) 'j_avg not implemented in value_at_diff_schr!'
    end if
  end subroutine value_at_diff_schr

  subroutine normalize_diff_schr(this,y)
    use utils,only:fp,ip
    use linfuncs
    implicit none
    class(diff_schr),intent(in) :: this
    complex(fp),dimension(:),intent(inout),target :: y

    y = y/l2norm(y)
  end subroutine normalize_diff_schr
!##############################################################
  subroutine init_diff_sme(this,H)
    use utils,only:fp,ip
    use hamiltonians
    implicit none
    class(diff_sme),intent(inout) :: this
    class(hamiltonian), target :: H
    integer(ip),parameter :: initial_lo_cont_size = 2

    this%H => H

    allocate(this%lindblad_operators(this%H%dim,this%H%dim,initial_lo_cont_size))
    allocate(this%lindblad_operator_couplings(initial_lo_cont_size))
    allocate(this%lo_comb(this%H%dim,this%H%dim),this%loclo_comb(this%H%dim,this%H%dim))
    this%nlo = 0
    this%lo_cont_size = initial_lo_cont_size

  end subroutine init_diff_sme
!---------------------------------------------------------------!
  subroutine value_at_diff_sme(this,x,y,d,s,j_avg)
    use utils,only:ip,fp
    use linfuncs
    implicit none
    class(diff_sme),intent(in) :: this
    integer(ip) :: dim
    real(fp),intent(in) :: x
    complex(fp),dimension(:),target,intent(in) :: y
    complex(fp),dimension(:,:),pointer :: rho
    complex(fp), dimension(size(y)),intent(out) :: d
    complex(fp), dimension(size(y),this%nobs),intent(out) :: s
    real(fp),dimension(this%nobs),intent(out),optional :: j_avg
    real(fp),dimension(this%nobs) :: j_avg_temp
    complex(fp),dimension(this%nobs) :: avg_a
    integer(ip) :: ii
    !complex(fp),dimension(this%H%dim,this%H%dim,this%nobs) :: Act
    complex(fp),dimension(this%H%dim,this%H%dim) :: temp
    complex(fp),dimension(this%H%dim,this%H%dim) :: rho_temp
    complex(fp),dimension(this%H%dim,this%H%dim) :: Arho,rhoA,Lrho,rhoL,Ht
    complex(fp),dimension(this%H%dim,this%H%dim) :: orc
    complex(fp),dimension(this%H%dim,this%H%dim) :: rco 
    complex(fp),dimension(this%H%dim,this%H%dim) ::cor
    complex(fp),dimension(this%H%dim,this%H%dim) ::rc
    complex(fp),dimension(this%H%dim,this%H%dim) ::or

    dim = this%H%dim
    rho(1:dim,1:dim) => y
    ! Source: http://math.univ-lyon1.fr/~attal/Mesarticles/Stoheat.pdf
    ! dp_t = L(p_t) dt + (C p_t + p_tC* - Tr[p_t(C+C*)]p_t)*dW_t

    Ht = this%H%value_at(x)
    ! Calculates only Lrho - rhoL
    Lrho = 0.0_fp
    rhoL = 0.0_fp
    d = 0.0_fp
    s = 0.0_fp
    Arho = 0.0_fp
    rhoA = 0.0_fp
    orc = 0.0_fp
    rco = 0.0_fp
    cor = 0.0_fp
    rc = 0.0_fp
    or = 0.0_fp
    call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),Ht,&
         dim,rho,dim,(0.0_fp,0.0_fp),Lrho,dim)
    call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),rho,&
         dim,Ht,dim,(0.0_fp,0.0_fp),rhoL,dim)
    d = reshape(Lrho - rhoL,[size(y)])

    call zgemm('n','c',dim,dim,dim,(1.0_fp,0.0_fp),rho,dim,&
         this%lo_comb,dim,(0_fp,0_fp),rc,dim)
    call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),this%lo_comb,dim,&
         rc,dim,(0_fp,0_fp),orc,dim)
    call zgemm('n','n',dim,dim,dim,(0.5_fp,0.0_fp),this%loclo_comb,dim,&
         rho,dim,(0_fp,0_fp),cor,dim)
    call zgemm('n','n',dim,dim,dim,(0.5_fp,0.0_fp),rho,dim,&
         this%loclo_comb,dim,(0_fp,0_fp),rco,dim)
    d = d + reshape(orc - cor - rco,[size(y)])
    ! do ii=1,this%nlo
    !    c = conjg(transpose(this%lindblad_operators(:,:,ii)))
    !    rc = (0.0_fp,0.0_fp)
    !    or = (0.0_fp,0.0_fp)
    !    rco = (0.0_fp,0.0_fp)
    !    cor = (0.0_fp,0.0_fp)
    !    orc = (0.0_fp,0.0_fp)
    !    call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),rho,dim,&
    !         c,dim,(0_fp,0_fp),rc,dim)
    !    call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),this%lindblad_operators(:,:,ii),dim,&
    !         rho,dim,(0_fp,0_fp),or,dim)
    !    call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),this%lindblad_operators(:,:,ii),dim,&
    !         rc,dim,(0_fp,0_fp),orc,dim)
    !    call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),rc,dim,&
    !         this%lindblad_operators(:,:,ii),dim,(0_fp,0_fp),rco,dim)
    !    call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),c,dim,&
    !         or,dim,(0_fp,0_fp),cor,dim)
    !    !call print_matrix(c)
    !    d  = d + this%lindblad_operator_couplings(ii)*reshape(orc - 0.5_fp*cor - 0.5_fp*rco,[size(y)])
    ! end do
    do ii=1,this%nobs
       call zgemm('n','n',dim,dim,dim,(1.0_fp,0.0_fp),this%observables(:,:,ii),&
            dim,rho,dim,(0.0_fp,0.0_fp),Arho,dim)
       call zgemm('n','c',dim,dim,dim,(1.0_fp,0.0_fp),rho,&
            dim,this%observables(:,:,ii),dim,(0.0_fp,0.0_fp),rhoA,dim)
       !call print_matrix(Arho)
       !call print_matrix(rhoA)
       j_avg_temp(ii) = trace(Arho+rhoA)
       s(:,ii) = sqrt(this%couplings(ii))*reshape(Arho + rhoA - j_avg_temp(ii)*rho,[size(y)])
       if(present(j_avg)) then
             j_avg(ii) = sqrt(this%couplings(ii))*j_avg_temp(ii)
       end if
    end do
    !write(*,*) 'd real', real(d)
    !write(*,*) 'd imag', aimag(d)
    !write(*,*) 's real', real(s)
    !write(*,*) 's imag', aimag(s)
    !write(*,*) ''
  end subroutine value_at_diff_sme
!---------------------------------------------------------------!
  subroutine normalize_diff_sme(this,y)
    use utils,only:fp,ip
    use linfuncs
    implicit none
    class(diff_sme),intent(in) :: this
    complex(fp),dimension(:),target,intent(inout) :: y
    complex(fp),dimension(:,:),pointer :: rho
    rho(1:this%H%dim,1:this%H%dim) => y
    rho = rho/trace_abs(rho)
    !    call print_matrix(rho)
    ! if(maxval(abs(y))>1.001) then
    !    write(*,*) 'Trace[rho]=',trace(rho),', Maxval=',maxval(abs(y))
    ! end if
  end subroutine normalize_diff_sme

  subroutine add_lindblad_operator_diff_sme(this,op,coupling)
    implicit none
    class(diff_sme),intent(inout) :: this
    complex(fp),intent(in) :: coupling
    complex(fp),allocatable,dimension(:,:,:) :: temp_cont
    complex(fp),allocatable,dimension(:) :: temp_cont_couplings
    complex(fp),dimension(this%H%dim,this%H%dim),intent(in) :: op
    integer(ip) :: dim,ii
    dim = this%H%dim

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

    this%lo_comb = sqrt(this%lindblad_operator_couplings(1))*this%lindblad_operators(:,:,1)
    do ii=2,this%nlo
       this%lo_comb = this%lo_comb + sqrt(this%lindblad_operator_couplings(ii))*&
            this%lindblad_operators(:,:,ii)
    end do
    call zgemm('c','n',dim,dim,dim,(1.0_fp,0.0_fp),this%lo_comb,dim,this%lo_comb,dim,(0.0_fp,0.0_fp),this%loclo_comb,dim)
  end subroutine add_lindblad_operator_diff_sme

  subroutine add_observable_diff_sme(this,A,coupling,efficiency_in)
    implicit none
    class(diff_sme) :: this
    complex(fp),dimension(:,:),intent(in) :: A
    complex(fp),intent(in) :: coupling
    integer(ip),parameter :: obs_cont_initial_size = 5
    complex(fp),allocatable,dimension(:,:,:) :: observables_t
    complex(fp),allocatable,dimension(:) :: couplings_t
    real(fp),intent(in),optional :: efficiency_in
    real(fp) :: efficiency

    if(present(efficiency_in)) then
       efficiency = efficiency_in
    else
       efficiency = 1.0_fp
    end if
    call this%add_lindblad_operator(A,coupling)
    call this%sd%add_observable(A,coupling,efficiency)
  end subroutine add_observable_diff_sme

  subroutine set_lindblad_coupling_diff_sme(this,coupling,n)
    implicit none
    class(diff_sme) :: this
    complex(fp),intent(in) :: coupling
    integer(ip),intent(in) :: n
    integer(ip) :: dim,ii
    if(this%nlo < n) then
       write(*,*) 'Trying to set the coupling of nonexisting operator.'
       return
    end if
    dim = this%H%dim    
    this%lindblad_operator_couplings(n) = coupling

    this%lo_comb = sqrt(this%lindblad_operator_couplings(1))*this%lindblad_operators(:,:,1)
    do ii=2,this%nlo
       this%lo_comb = this%lo_comb + sqrt(this%lindblad_operator_couplings(ii))*&
            this%lindblad_operators(:,:,ii)
    end do
    call zgemm('c','n',dim,dim,dim,(1.0_fp,0.0_fp),this%lo_comb,dim,this%lo_comb,dim,(0.0_fp,0.0_fp),this%loclo_comb,dim)
  end subroutine set_lindblad_coupling_diff_sme
!################################################################
end module sclass
