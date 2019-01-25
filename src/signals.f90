module signals
  use utils, only: fp,ip
  use envelopes
  implicit none
  private
  public :: signal
  integer(ip),parameter :: initial_nbr_env = 2

  !type :: container
  !   class(envelope),pointer :: obj
  !end type container

  type :: signal
     real(fp) :: w
     complex(fp),dimension(:,:),allocatable :: opc
     complex(fp),dimension(:,:),allocatable :: opa
     integer(ip) :: dim,env_cont_size
     !integer(ip),pointer :: nenv
     integer(ip) :: nenv
     type(envelope_container), dimension(:), allocatable :: envelopes
   contains
     procedure :: initialize => initialize_signal
     procedure,private :: value_at_signal_wref
     procedure,private :: value_at_signal_wref_array
     generic,public :: value_at => value_at_signal_wref,value_at_signal_wref_array
     procedure :: add_envelope => add_envelope_signal
     procedure :: set_w => set_w_signal
     procedure :: get_nbr_of_envelopes => get_nbr_of_envelopes_signal
     !procedure :: get_nbr_of_envelopes_ptr => get_nbr_of_envelopes_ptr_signal
     procedure :: get_envelope => get_envelope_signal
     procedure :: get_envelopes => get_envelopes_signal
  end type signal

contains
  subroutine initialize_signal(this,op_in,w)
    use utils, only: fp,ip
    implicit none
    class(signal) :: this
    complex(fp),dimension(:,:),intent(in) :: op_in
    real(fp),intent(in) :: w
    this%dim = size(op_in,1)
    allocate(this%opc(this%dim,this%dim),this%opa(this%dim,this%dim))
    this%opc = op_in
    this%opa = conjg(transpose(op_in))
    this%nenv = 0
    allocate(this%envelopes(initial_nbr_env))
    this%env_cont_size = initial_nbr_env
    this%w = w
  end subroutine initialize_signal
  !-------------------------------------------------------------
  function value_at_signal_wref(this,t,wref_in) result(y)
    use utils, only:fp,ip
    class(signal), intent(in) :: this
    real(fp),intent(in) :: t
    real(fp),optional,intent(in) :: wref_in
    complex(fp),dimension(this%dim,this%dim) :: y
    integer(ip) :: ii
    real(fp) :: wref
    complex(fp) :: g
    if(present(wref_in)) then
       wref = wref_in
    else
       wref = 0.0_fp
    endif
    !y = (0.0_fp,0.0_fp)
    g = (0.0_fp,0.0_fp)
    do ii=1,this%nenv
       g = g + this%envelopes(ii)%obj%value_at(t)
    end do
    ! There is a mistake here: Usually the convention is to multiply conjg(g) with opc and not vice versa. The net effect is that all the phases are reversed (positives become negative).
    !y = (g*exp(-(0.0_fp,1.0_fp)*t*(this%w-wref))*this%opc + &
    !     conjg(g)*exp((0.0_fp,1.0_fp)*t*(this%w-wref))*this%opa)
    ! fixed version
    y = (conjg(g)*exp(-(0.0_fp,1.0_fp)*t*(this%w-wref))*this%opc + &
         g*exp((0.0_fp,1.0_fp)*t*(this%w-wref))*this%opa)
    !write(*,*) this%w - wref
  end function value_at_signal_wref
  !-------------------------------------------------------------
  function get_nbr_of_envelopes_signal(this) result(n)
    use utils,only:fp,ip
    class(signal),intent(in) :: this
    integer(ip) :: n
    n = this%nenv
  end function get_nbr_of_envelopes_signal
  !-------------------------------------------------------------
  ! function get_nbr_of_envelopes_ptr_signal(this) result(n)
  !   use utils,only:fp,ip
  !   class(signal),intent(in) :: this
  !   integer(ip),pointer :: n
  !   n => this%nenv
  ! end function get_nbr_of_envelopes_ptr_signal
  !------------------------------------------------------------
  function value_at_signal_wref_array(this,t,wref_array) result(y)
    use utils, only:fp,ip
    use linfuncs
    class(signal), intent(in) :: this
    real(fp),intent(in) :: t
    real(fp),dimension(:),intent(in) :: wref_array
    complex(fp),dimension(this%dim,this%dim) :: y
    integer(ip) :: ii,jj
    complex(fp) :: g
    real(fp) :: w
    !y = (0.0_fp,0.0_fp)
    g = (0.0_fp,0.0_fp)
    do ii=1,this%nenv
       g = g + this%envelopes(ii)%obj%value_at(t)
    end do
    !write(*,*) wref_array
    do ii=1,this%dim
       do jj=1,this%dim
          if(ii>jj) then
             !w = this%w + (wref_array(jj) - wref_array(ii))
             !fixed version
             w = this%w + (wref_array(jj) - wref_array(ii))
          else
             !w = this%w - (wref_array(jj) - wref_array(ii))
             !fixed version
             w = this%w - (wref_array(jj) - wref_array(ii))
          endif
    ! There is a mistake here: Usually the convention is to multiply conjg(g) with opc and not vice versa. Also the order (this%w -wref) is non-convential. The net effect is that all the phases are reversed (positives become negative). If you fix this, remember to also change value_at_signal_wref_array. The problem was actually only in conjg(g). The order  (this%w -wref) is correct.
          !y(jj,ii) = (g*exp(-(0.0_fp,1.0_fp)*t*w)*this%opc(jj,ii) + conjg(g)*exp((0.0_fp,1.0_fp)*t*w)*this%opa(jj,ii))
          !fixed version
          y(jj,ii) = (conjg(g)*exp(-(0.0_fp,1.0_fp)*t*w)*this%opc(jj,ii) + g*exp((0.0_fp,1.0_fp)*t*w)*this%opa(jj,ii))
          !write(*,*) 'w', ii,jj,':',w
          !y = g*(exp(-(0.0_fp,1.0_fp)*t*(this%w-wref))*this%opc + &
          !     exp((0.0_fp,1.0_fp)*t*(this%w-wref))*this%opa)
       end do
    end do
    !call print_matrix(real(this%opc))
    !write(*,*) this%w - wref,y
  end function value_at_signal_wref_array
  !------------------------------------------------------------
  subroutine add_envelope_signal(this,envelope_in)
    use utils,only:fp,ip
    class(signal),intent(inout) :: this
    class(envelope),intent(in),target :: envelope_in
    type(envelope_container),allocatable,dimension(:) :: temp_cont
    if(this%nenv >= this%env_cont_size) then
       allocate(temp_cont(this%nenv))
       temp_cont = this%envelopes(1:this%nenv)
       deallocate(this%envelopes)
       allocate(this%envelopes(this%nenv*2))
       this%envelopes(1:this%nenv) = temp_cont
       deallocate(temp_cont)
       this%env_cont_size = size(this%envelopes)
       !write(*,*) 'Extending space for envelopes. New size is', this%env_cont_size,'.'
    end if
    this%nenv = this%nenv + 1
    this%envelopes(this%nenv)%obj => envelope_in
  end subroutine add_envelope_signal

  subroutine set_w_signal(this,w_in)
    use utils,only:fp,ip
    class(signal),intent(inout) :: this
    real(fp),intent(in) :: w_in
    this%w = w_in
  end subroutine set_w_signal

  function get_envelope_signal(this,n) result(env)
    use utils,only:fp,ip
    class(signal),intent(in) :: this
    integer(ip),intent(in) :: n
    class(envelope),pointer :: env
    env => this%envelopes(n)%obj
  end function get_envelope_signal

  subroutine get_envelopes_signal(this,env)
    use utils,only:fp,ip
    class(signal),intent(in) :: this
    type(envelope_container), dimension(:),intent(out) :: env
    env = this%envelopes
  end subroutine get_envelopes_signal

end module signals
