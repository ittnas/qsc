module envelopes
  use utils,only: fp,ip
  use constants,only: pi
  
  implicit none

  type, abstract :: envelope

   contains
     procedure(value_at_interface),deferred :: value_at
     procedure(set_parameters_interface),deferred :: set_parameters
  end type envelope

  abstract interface
     function value_at_interface(this,t) result(y)
       use utils,only:fp,ip
       import envelope
       class(envelope), intent(in) :: this
       real(fp), intent(in) :: t
       complex(fp) :: y
     end function value_at_interface
  end interface

  abstract interface
     subroutine set_parameters_interface(this,parameter_array)
       use utils,only:fp,ip
       import envelope
       class(envelope), intent(inout) :: this
       complex(fp),dimension(:),intent(in) :: parameter_array
     end subroutine set_parameters_interface
  end interface

  type :: envelope_container
     class(envelope),pointer :: obj
  end type envelope_container

  type,extends(envelope) :: constant_envelope
  complex(fp) :: A
contains
  procedure :: value_at=>value_at_constant_envelope
  procedure :: initialize=>initialize_constant_envelope
  procedure :: set_parameters=>set_parameters_constant_envelope
end type constant_envelope

type,extends(envelope) :: gaussian_envelope
   complex(fp) :: A
   real(fp) :: sigma,t0
   real(fp) :: cutoff
   real(fp) :: trnc_point,trnc_point_s
 contains
   procedure :: value_at=>value_at_gaussian_envelope
   procedure :: initialize=>initialize_gaussian_envelope
   procedure :: set_t0=>set_t0_gaussian_envelope
   procedure :: set_A=>set_A_gaussian_envelope
   procedure :: set_parameters=>set_parameters_gaussian_envelope
end type gaussian_envelope

type,extends(envelope) :: rectangular_envelope
   complex(fp) :: A
   real(fp) :: duration,t0
 contains
   procedure :: value_at=>value_at_rectangular_envelope
   procedure :: initialize=>initialize_rectangular_envelope
   procedure :: set_parameters=>set_parameters_rectangular_envelope
end type rectangular_envelope

type,extends(envelope) :: dummy_envelope
 contains
   procedure :: value_at=>value_at_dummy_envelope
   procedure :: initialize=>initialize_dummy_envelope
   procedure :: set_parameters=>set_parameters_dummy_envelope
end type dummy_envelope


type,extends(envelope) :: lorenzian_envelope
   complex(fp) :: A
   real(fp) :: Gamma,t0
 contains
   procedure :: value_at=>value_at_lorenzian_envelope
   procedure :: initialize=>initialize_lorenzian_envelope
   procedure :: set_parameters=>set_parameters_lorenzian_envelope
end type lorenzian_envelope

type,extends(envelope) :: flat_gaussian_envelope
   complex(fp) :: A
   real(fp) :: sigma,t0,saturation_start,saturation_end
 contains
   procedure :: value_at=>value_at_flat_gaussian_envelope
   procedure :: initialize=>initialize_flat_gaussian_envelope
   procedure :: set_parameters=>set_parameters_flat_gaussian_envelope
end type flat_gaussian_envelope

type,extends(envelope) :: inverse_cosh_envelope
   complex(fp) :: A
   real(fp) :: ts,t0,sigma,cutoff,trnc_point
 contains
   procedure :: value_at=>value_at_inverse_cosh_envelope
   procedure :: initialize=>initialize_inverse_cosh_envelope
   procedure :: set_parameters=>set_parameters_inverse_cosh_envelope
end type inverse_cosh_envelope

type,extends(inverse_cosh_envelope) :: inverse_sqrt_cosh_envelope
 contains
   procedure :: value_at=>value_at_inverse_sqrt_cosh_envelope
   procedure :: initialize=>initialize_inverse_sqrt_cosh_envelope
   procedure :: initialize_inverse_sqrt_cosh_envelope
   procedure :: set_parameters=>set_parameters_inverse_sqrt_cosh_envelope
end type inverse_sqrt_cosh_envelope

type,extends(inverse_sqrt_cosh_envelope) :: inverse_sqrt_cosh_envelope_dfa
   class(envelope),pointer :: target_envelope
   real(fp) :: eps,delta,tis,dt,g
 contains
   procedure :: value_at=>value_at_inverse_sqrt_cosh_envelope_dfa
   procedure :: initialize_inverse_sqrt_cosh_envelope_dfa
   procedure :: set_parameters=>set_parameters_inverse_sqrt_cosh_envelope_dfa
end type inverse_sqrt_cosh_envelope_dfa

type,extends(envelope) :: general_dfa
   real(fp) :: eps,delta,tis,dt,g
   class(envelope),pointer :: shift_envelope
   class(envelope),pointer :: target_envelope
 contains
   procedure :: value_at=>value_at_general_dfa
   procedure :: initialize=>initialize_general_dfa
   procedure :: set_parameters=>set_parameters_general_dfa
   procedure :: set_delta => set_delta_general_dfa
end type general_dfa

type,extends(envelope) :: general_dfa_01
   real(fp) :: eps,delta,tis,dt,g
   class(envelope),pointer :: shift_envelope
   class(envelope),pointer :: target_envelope
 contains
   procedure :: value_at=>value_at_general_dfa_01
   procedure :: initialize=>initialize_general_dfa_01
   procedure :: set_parameters=>set_parameters_general_dfa_01
   procedure :: set_delta => set_delta_general_dfa_01
end type general_dfa_01

type,extends(envelope) :: general_dfa_12
   real(fp) :: eps,delta,tis,dt,g
   class(envelope),pointer :: shift_envelope
   class(envelope),pointer :: target_envelope
 contains
   procedure :: value_at=>value_at_general_dfa_12
   procedure :: initialize=>initialize_general_dfa_12
   procedure :: set_parameters=>set_parameters_general_dfa_12
   procedure :: set_delta => set_delta_general_dfa_12
end type general_dfa_12

type,extends(envelope) :: gaussian_linear_fractional_stirap
   complex(fp) :: A
   real(fp) :: sigma,t0,t2
   real(fp) :: cutoff
   real(fp) :: trnc_point
 contains
   procedure :: value_at=>value_at_gaussian_linear_fractional_stirap
   procedure :: initialize=>initialize_gaussian_linear_fractional_stirap
   procedure :: set_t0=>set_t0_gaussian_linear_fractional_stirap
   procedure :: set_A=>set_A_gaussian_linear_fractional_stirap
   procedure :: set_parameters=>set_parameters_gaussian_linear_fractional_stirap
end type gaussian_linear_fractional_stirap


type,extends(envelope) :: linear_zero_area_pulse
   real(fp) :: t0,tf
   complex(fp) :: A
 contains
   procedure :: value_at=>value_at_linear_zero_area_pulse
   procedure :: initialize=>initialize_linear_zero_area_pulse
   procedure :: set_parameters=>set_parameters_linear_zero_area_pulse
end type linear_zero_area_pulse

type,extends(envelope) :: two_photon_linear_zero_area_pulse
   real(fp) :: t0,tf
   complex(fp) :: A
   real(fp) :: phase_shift
 contains
   procedure :: value_at=>value_at_two_photon_linear_zero_area_pulse
   procedure :: initialize=>initialize_two_photon_linear_zero_area_pulse
   procedure :: set_parameters=>set_parameters_two_photon_linear_zero_area_pulse
end type two_photon_linear_zero_area_pulse

type,extends(envelope) :: two_photon_linear_trnc_pulse
   real(fp) :: t0,tf
   complex(fp) :: A,A2
   real(fp) :: phase_shift
 contains
   procedure :: value_at=>value_at_two_photon_linear_trnc_pulse
   procedure :: initialize=>initialize_two_photon_linear_trnc_pulse
   procedure :: set_parameters=>set_parameters_two_photon_linear_trnc_pulse
end type two_photon_linear_trnc_pulse

type,extends(envelope) :: exp_polinom
   real(fp) :: t0,td,trnc_point,cutoff,n
   complex(fp) :: A
 contains
   procedure :: value_at=>value_at_exp_polinom
   procedure :: initialize=>initialize_exp_polinom
   procedure :: set_parameters=>set_parameters_exp_polinom
end type exp_polinom

type,extends(envelope) :: polynomial_envelope
   complex(fp),dimension(:),allocatable :: A
   real(fp) :: trnc_point
 contains
   procedure :: value_at=>value_at_polynomial_envelope
   procedure :: initialize=>initialize_polynomial_envelope
   procedure :: set_parameters=>set_parameters_polynomial_envelope
end type polynomial_envelope

type,extends(envelope) :: sqrt_polynomial_envelope
   complex(fp),dimension(:),allocatable :: A
   real(fp) :: trnc_point
 contains
   procedure :: value_at=>value_at_sqrt_polynomial_envelope
   procedure :: initialize=>initialize_sqrt_polynomial_envelope
   procedure :: set_parameters=>set_parameters_sqrt_polynomial_envelope
end type sqrt_polynomial_envelope


type,extends(envelope) :: lr_envelope_01
   class(envelope),pointer :: alpha,beta,dalpha,dbeta,Omega02
   complex(fp) :: factor
   real(fp) :: trnc_point
 contains
   procedure :: value_at=>value_at_lr_envelope_01
   procedure :: initialize=>initialize_lr_envelope_01
   procedure :: set_parameters=>set_parameters_lr_envelope_01
end type lr_envelope_01

type,extends(envelope) :: lr_envelope_12
   class(envelope),pointer :: alpha,beta,dalpha,dbeta,Omega02
   complex(fp) :: factor
   real(fp) :: trnc_point
 contains
   procedure :: value_at=>value_at_lr_envelope_12
   procedure :: initialize=>initialize_lr_envelope_12
   procedure :: set_parameters=>set_parameters_lr_envelope_12
end type lr_envelope_12

type,extends(envelope) :: sa_stirap_det_02
   real(fp) :: dt,delta,phi,A_factor,csr
   type(envelope_container),dimension(:),allocatable :: d01,d12
   integer(ip) :: nenv01,nenv12,env_cont_size_01,env_cont_size_12

 contains
   procedure :: value_at=>value_at_sa_stirap_det_02
   procedure :: initialize=>initialize_sa_stirap_det_02
   procedure :: set_parameters=>set_parameters_sa_stirap_det_02
   procedure :: add_envelope_01 => add_envelope_01_sa_stirap_det_02
   procedure :: add_envelope_12 => add_envelope_12_sa_stirap_det_02
   !procedure :: set_delta => set_delta_sa_stirap_det_02
end type sa_stirap_det_02

type,extends(envelope) :: sa_stirap_det_01
   real(fp) :: delta,dt
   class(envelope),pointer :: e01,e12
 contains
   procedure :: value_at=>value_at_sa_stirap_det_01
   procedure :: initialize=>initialize_sa_stirap_det_01
   procedure :: set_parameters=>set_parameters_sa_stirap_det_01
   procedure :: set_delta => set_delta_sa_stirap_det_01
end type sa_stirap_det_01

type,extends(envelope) :: sa_stirap_det_12
   real(fp) :: delta,dt
   class(envelope),pointer :: e01,e12
 contains
   procedure :: value_at=>value_at_sa_stirap_det_12
   procedure :: initialize=>initialize_sa_stirap_det_12
   procedure :: set_parameters=>set_parameters_sa_stirap_det_12
   procedure :: set_delta => set_delta_sa_stirap_det_12
end type sa_stirap_det_12

contains
  !#########################################################
  function value_at_constant_envelope(this,t) result(y)
    use utils,only:fp,ip
    class(constant_envelope), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    y = this%A
  end function value_at_constant_envelope
  !---------------------------------------------------------
  subroutine initialize_constant_envelope(this,A_in)
    use utils,only:fp,ip
    class(constant_envelope) :: this
    complex(fp), intent(in) :: A_in
    this%A = A_in
  end subroutine initialize_constant_envelope
  !---------------------------------------------------------
  subroutine set_parameters_constant_envelope(this,parameter_array)
    use utils,only:fp,ip
    class(constant_envelope), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) < 1) then
       write(*,*) 'Number of parameters incorrect in set_parameters_constant_envelope.'
       call exit(1)
    end if
    this%A = parameter_array(1)
  end subroutine set_parameters_constant_envelope
  !#########################################################
  function value_at_gaussian_envelope(this,t) result(y)
    use utils,only:fp,ip
    class(gaussian_envelope), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    real(fp) :: value_at_cutoff
    if(this%cutoff<=0.0_fp) then
       y = this%A*exp(-(t-this%t0)**2/(2*this%sigma**2))
    else
       if(t < (this%t0 - this%cutoff*this%sigma) .or. t > (this%t0 + this%cutoff*this%sigma) .or. &
            t > this%trnc_point .or. t < this%trnc_point_s) then
          y = 0.0_fp
       else
          value_at_cutoff = exp(-this%cutoff**2/2.0_fp)
          y = this%A/(1 - value_at_cutoff)*(exp(-(t-this%t0)**2/(2*this%sigma**2)) - value_at_cutoff)
       end if
    end if
  end function value_at_gaussian_envelope
  !---------------------------------------------------------
  subroutine initialize_gaussian_envelope(this,A,sigma,t0,cutoff,trnc_point,trnc_point_s)
    use utils,only:fp,ip
    implicit none
    class(gaussian_envelope) :: this
    complex(fp), intent(in) :: A
    real(fp), intent(in) :: sigma,t0
    real(fp),optional,intent(in) :: cutoff,trnc_point,trnc_point_s
    this%A = A
    this%sigma = sigma
    this%t0 = t0
    if(present(cutoff)) then
       this%cutoff = cutoff
    else
       this%cutoff = -1.0_fp
    end if
    if(present(trnc_point)) then
       this%trnc_point = trnc_point
    else
       this%trnc_point = huge(trnc_point)
    end if

    if(present(trnc_point_s)) then
       this%trnc_point_s = trnc_point_s
    else
       this%trnc_point_s = -huge(trnc_point)
    end if
  end subroutine initialize_gaussian_envelope
  !---------------------------------------------------------
  subroutine set_t0_gaussian_envelope(this,t0)
    implicit none
    class(gaussian_envelope) :: this
    real(fp),intent(in) :: t0
    this%t0 = t0
  end subroutine set_t0_gaussian_envelope
  !---------------------------------------------------------
  subroutine set_A_gaussian_envelope(this,A)
    implicit none
    class(gaussian_envelope) :: this
    real(fp),intent(in) :: A
    this%A = A
  end subroutine set_A_gaussian_envelope
  !---------------------------------------------------------
  subroutine set_parameters_gaussian_envelope(this,parameter_array)
    use utils,only:fp,ip
    class(gaussian_envelope), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) < 3) then
       write(*,*) 'Number of parameters incorrect in set_parameters_gaussian_envelope.'
       call exit(1)
    end if
    this%A = parameter_array(1)
    this%sigma = real(parameter_array(2))
    this%t0 = real(parameter_array(3))
    if(size(parameter_array) == 4) then
       this%cutoff = real(parameter_array(4))
    end if
    if(size(parameter_array) == 5) then
       this%trnc_point = real(parameter_array(5))
    end if
    if(size(parameter_array) == 6) then
       this%trnc_point_s = real(parameter_array(6))
    end if
  end subroutine set_parameters_gaussian_envelope
    !#########################################################
  function value_at_gaussian_linear_fractional_stirap(this,t) result(y)
    use utils,only:fp,ip
    class(gaussian_linear_fractional_stirap), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    real(fp) :: f
    if(t>this%trnc_point) then
       if(t<this%trnc_point+this%t2) then
          f = this%A*exp(-(this%trnc_point-this%t0)**2/(2*this%sigma**2))
          y = -f/this%t2*t + f + f*this%trnc_point/this%t2
       else
          y = 0.0_fp
       endif
    else
       y = this%A*exp(-(t-this%t0)**2/(2*this%sigma**2))
    endif
    !write(*,*) y
  end function value_at_gaussian_linear_fractional_stirap
  !---------------------------------------------------------
  subroutine initialize_gaussian_linear_fractional_stirap(this,A,sigma,t0,cutoff,trnc_point,t2)
    use utils,only:fp,ip
    implicit none
    class(gaussian_linear_fractional_stirap) :: this
    complex(fp), intent(in) :: A
    real(fp), intent(in) :: sigma,t0
    real(fp),optional,intent(in) :: cutoff,trnc_point,t2
    this%A = A
    this%sigma = sigma
    this%t0 = t0
    if(present(cutoff)) then
       this%cutoff = cutoff
    else
       this%cutoff = -1.0_fp
    end if
    if(present(trnc_point)) then
       this%trnc_point = trnc_point
    else
       this%trnc_point = huge(trnc_point)
    end if
    if(present(t2)) then
       this%t2 = t2
    else
       this%t2 = 0.0_fp
    end if

  end subroutine initialize_gaussian_linear_fractional_stirap
  !---------------------------------------------------------
  subroutine set_t0_gaussian_linear_fractional_stirap(this,t0)
    implicit none
    class(gaussian_linear_fractional_stirap) :: this
    real(fp),intent(in) :: t0
    this%t0 = t0
  end subroutine set_t0_gaussian_linear_fractional_stirap
  !---------------------------------------------------------
  subroutine set_A_gaussian_linear_fractional_stirap(this,A)
    implicit none
    class(gaussian_linear_fractional_stirap) :: this
    real(fp),intent(in) :: A
    this%A = A
  end subroutine set_A_gaussian_linear_fractional_stirap
  !---------------------------------------------------------
  subroutine set_parameters_gaussian_linear_fractional_stirap(this,parameter_array)
    use utils,only:fp,ip
    class(gaussian_linear_fractional_stirap), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) < 3) then
       write(*,*) 'Number of parameters incorrect in set_parameters_gaussian_envelope.'
       call exit(1)
    end if
    this%A = parameter_array(1)
    this%sigma = real(parameter_array(2))
    this%t0 = real(parameter_array(3))
    if(size(parameter_array) == 4) then
       this%cutoff = real(parameter_array(4))
    end if
    if(size(parameter_array) == 5) then
       this%trnc_point = real(parameter_array(5))
    end if
    if(size(parameter_array) == 6) then
       this%t2 = real(parameter_array(6))
    end if
  end subroutine set_parameters_gaussian_linear_fractional_stirap
  !#########################################################
  function value_at_rectangular_envelope(this,t) result(y)
    use utils,only:fp,ip
    implicit none
    class(rectangular_envelope), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    if(t >= this%t0 .and. (t<= this%t0 + this%duration)) then
       y = this%A
    else
       y = (0.0_fp,0.0_fp)
    endif
  end function value_at_rectangular_envelope
  !---------------------------------------------------------
  subroutine initialize_rectangular_envelope(this,A,duration,t0)
    use utils,only:fp,ip
    class(rectangular_envelope) :: this
    complex(fp), intent(in) :: A
    real(fp), intent(in) :: duration,t0
    this%A = A
    this%duration = duration
    this%t0 = t0
  end subroutine initialize_rectangular_envelope
  !---------------------------------------------------------
  subroutine set_parameters_rectangular_envelope(this,parameter_array)
    use utils,only:fp,ip
    class(rectangular_envelope), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) < 3) then
       write(*,*) 'Number of parameters incorrect in set_parameters_gaussian_envelope.'
       call exit(1)
    end if
    this%A = parameter_array(1)
    this%duration = real(parameter_array(2))
    this%t0 = real(parameter_array(3))
  end subroutine set_parameters_rectangular_envelope
  !#########################################################
  function value_at_dummy_envelope(this,t) result(y)
    use utils,only:fp,ip
    implicit none
    class(dummy_envelope), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    y = (0.0_fp,0.0_fp)
  end function value_at_dummy_envelope
  !---------------------------------------------------------
  subroutine initialize_dummy_envelope(this)
    use utils,only:fp,ip
    class(dummy_envelope) :: this
  end subroutine initialize_dummy_envelope
  !---------------------------------------------------------
  subroutine set_parameters_dummy_envelope(this,parameter_array)
    use utils,only:fp,ip
    class(dummy_envelope), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
  end subroutine set_parameters_dummy_envelope
  !#########################################################
  function value_at_lorenzian_envelope(this,t) result(y)
    use utils,only:fp,ip
    use constants,only:pi
    implicit none
    class(lorenzian_envelope), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    !y=this%A*0.5_fp/pi*this%Gamma/((t-this%t0)**2+0.5_fp*this%Gamma**2) ! Definition from Wolfram
    !y=this%A*0.5_fp*this%Gamma**2/((t-this%t0)**2+0.5_fp*this%Gamma**2) ! Max(y) = A
    !y=this%A*this%Gamma**2/((t-this%t0)**2+this%Gamma**2) ! Max(y) = A, Gamma = 0.5*FWHM
    y=this%A*this%Gamma**2/((t-this%t0)**2+this%Gamma**2) - this%A*this%Gamma**2/((0.0_fp-this%t0)**2+this%Gamma**2) ! Max(y) = A, Gamma = 0.5*FWHM
    !write(*,*) y
  end function value_at_lorenzian_envelope
  !---------------------------------------------------------
  subroutine initialize_lorenzian_envelope(this,A,Gamma,t0)
    use utils,only:fp,ip
    class(lorenzian_envelope) :: this
    complex(fp), intent(in) :: A
    real(fp), intent(in) :: Gamma,t0
    this%A = A
    this%Gamma = Gamma
    this%t0 = t0
  end subroutine initialize_lorenzian_envelope
  !---------------------------------------------------------
  subroutine set_parameters_lorenzian_envelope(this,parameter_array)
    use utils,only:fp,ip
    class(lorenzian_envelope), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) < 3) then
       write(*,*) 'Number of parameters incorrect in set_parameters_gaussian_envelope.'
       call exit(1)
    end if
    this%A = parameter_array(1)
    this%Gamma = real(parameter_array(2))
    this%t0 = real(parameter_array(3))
  end subroutine set_parameters_lorenzian_envelope
  !#########################################################
  function value_at_inverse_cosh_envelope(this,t) result(y)
    use utils,only:fp,ip
    use constants,only:pi
    implicit none
    class(inverse_cosh_envelope), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    if(t<this%t0-this%cutoff) then
       y = 0.0_fp
    elseif(t > this%t0+this%cutoff .or. t > this%trnc_point) then
       y = 0.0_fp
    else
       !y = this%A*this%ts/(this%sigma**2*cosh(this%ts*(t-this%t0)/this%sigma**2))
       y = this%A/cosh(this%ts*(t-this%t0)/this%sigma**2)
    end if
  end function value_at_inverse_cosh_envelope
  !---------------------------------------------------------
  subroutine initialize_inverse_cosh_envelope(this,A,ts,t0,sigma,cutoff,trnc_point)
    use utils,only:fp,ip
    class(inverse_cosh_envelope) :: this
    complex(fp), intent(in) :: A
    real(fp), intent(in) :: ts,sigma,t0
    real(fp), optional,intent(in) :: cutoff,trnc_point
    this%A = A
    this%ts = ts
    this%sigma = sigma
    this%t0 = t0
    if(present(cutoff)) then
       this%cutoff = cutoff
    else
       this%cutoff = huge(1.0_fp)
    end if

    if(present(trnc_point)) then
       this%trnc_point = trnc_point
    else
       this%trnc_point = huge(1.0_fp)
    end if
  end subroutine initialize_inverse_cosh_envelope
  !---------------------------------------------------------
  subroutine set_parameters_inverse_cosh_envelope(this,parameter_array)
    use utils,only:fp,ip
    class(inverse_cosh_envelope), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) > 0) then
       this%A = parameter_array(1)
    endif
    if(size(parameter_array) > 1) then
       this%ts = real(parameter_array(2))
    endif
    if(size(parameter_array) > 2) then
       this%sigma = real(parameter_array(3))
    endif
    if(size(parameter_array) > 3) then
       this%t0 = real(parameter_array(4))
    end if
    if(size(parameter_array) > 4) then
       this%cutoff = real(parameter_array(5))
    end if
    if(size(parameter_array) > 5) then
       this%trnc_point = real(parameter_array(6))
    end if
  end subroutine set_parameters_inverse_cosh_envelope
  !#########################################################
  function value_at_inverse_sqrt_cosh_envelope(this,t) result(y)
    use utils,only:fp,ip
    use constants,only:pi
    implicit none
    class(inverse_sqrt_cosh_envelope), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    if(abs(this%A)>0.0_fp) then
       y = sqrt(abs(this%inverse_cosh_envelope%value_at(t)/this%A))*this%A
    else
       y = 0.0_fp
    end if
  end function value_at_inverse_sqrt_cosh_envelope
  !---------------------------------------------------------
  subroutine initialize_inverse_sqrt_cosh_envelope(this,A,ts,t0,sigma,cutoff,trnc_point)
    use utils,only:fp,ip
    class(inverse_sqrt_cosh_envelope) :: this
    complex(fp), intent(in) :: A
    real(fp), intent(in) :: ts,sigma,t0
    real(fp), optional,intent(in) :: cutoff,trnc_point
    real(fp) :: cutoff_temp,trnc_point_temp
    if(present(cutoff)) then
       cutoff_temp = cutoff
    else
       cutoff_temp = huge(1.0_fp)
    end if
    if(present(trnc_point)) then
       trnc_point_temp = trnc_point
    else
       trnc_point_temp = huge(1.0_fp)
    end if

    call this%inverse_cosh_envelope%initialize(A,ts,t0,sigma,cutoff_temp,trnc_point_temp)
  end subroutine initialize_inverse_sqrt_cosh_envelope
  !---------------------------------------------------------
  subroutine set_parameters_inverse_sqrt_cosh_envelope(this,parameter_array)
    use utils,only:fp,ip
    class(inverse_sqrt_cosh_envelope), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    call this%inverse_cosh_envelope%set_parameters(parameter_array)
  end subroutine set_parameters_inverse_sqrt_cosh_envelope
  !#########################################################
  function value_at_inverse_sqrt_cosh_envelope_dfa(this,t) result(y)
    use utils,only:fp,ip
    use constants,only:pi
    use linfuncs,only:trapz
    implicit none
    class(inverse_sqrt_cosh_envelope_dfa), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    real(fp) :: Omega,dt_c,phi
    integer(ip) :: ii,nis
    real(fp),dimension(:),allocatable :: x,r

    nis = ceiling((t - this%tis)/this%dt)+1
    dt_c = (t - this%tis)/(nis-1)
    allocate(x(nis),r(nis))

    do ii=1,nis
       x(ii) = this%tis + (ii-1)*dt_c
       r(ii) = (abs(this%target_envelope%value_at(x(ii)))*2.0_fp)**2*(1-this%g**2)/(4.0_fp*this%delta)
    end do
    phi = trapz(x,r)
    !write(*,*) phi
    !write(*,*) r
    !write(*,*) nis,x(nis)
    !open(11,file="phi_test.dat",position="append")
    !write(11,*) phi
    !close(11)
    y = this%inverse_sqrt_cosh_envelope%value_at(t)*exp((0.0_fp,1.0_fp)*phi)
  end function value_at_inverse_sqrt_cosh_envelope_dfa
  !---------------------------------------------------------
  subroutine initialize_inverse_sqrt_cosh_envelope_dfa(this,A,ts,t0,sigma,cutoff,trnc_point,target_envelope, eps,delta,g,tis,dt)
    use utils,only:fp,ip
    class(inverse_sqrt_cosh_envelope_dfa) :: this
    class(envelope),intent(in),target :: target_envelope
    complex(fp), intent(in) :: A
    real(fp), intent(in) :: ts,sigma,t0
    real(fp), intent(in) :: eps,delta,tis,dt,g
    real(fp), optional,intent(in) :: cutoff,trnc_point
    real(fp) :: cutoff_temp,trnc_point_temp
    if(present(cutoff)) then
       cutoff_temp = cutoff
    else
       cutoff_temp = huge(1.0_fp)
    end if
    if(present(trnc_point)) then
       trnc_point_temp = trnc_point
    else
       trnc_point_temp = huge(1.0_fp)
    end if

    call this%inverse_sqrt_cosh_envelope%initialize(A,ts,t0,sigma,cutoff_temp,trnc_point_temp)
    this%target_envelope => target_envelope
    this%eps = eps
    this%delta = delta
    this%g = g
    this%tis = tis
    this%dt = dt
  end subroutine initialize_inverse_sqrt_cosh_envelope_dfa
  !---------------------------------------------------------
  subroutine set_parameters_inverse_sqrt_cosh_envelope_dfa(this,parameter_array)
    use utils,only:fp,ip
    class(inverse_sqrt_cosh_envelope_dfa), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    call this%inverse_sqrt_cosh_envelope%set_parameters(parameter_array)
  end subroutine set_parameters_inverse_sqrt_cosh_envelope_dfa
  !#########################################################
  function value_at_general_dfa(this,t) result(y)
    use utils,only:fp,ip
    use constants,only:pi
    use linfuncs,only:trapz
    implicit none
    class(general_dfa), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    real(fp) :: Omega,dt_c,phi
    integer(ip) :: ii,nis
    real(fp),dimension(:),allocatable :: x,r
    real(fp) :: int_r

    nis = ceiling((t - this%tis)/this%dt)+1
    dt_c = (t - this%tis)/(nis-1)
    allocate(x(nis),r(nis))

    int_r = 4.0_fp*((-1*this%g**2)/(4.0_fp*(this%eps-this%delta)) + (+1)/(4.0_fp*this%delta))


    do ii=1,nis
       x(ii) = this%tis + (ii-1)*dt_c
       !r(ii) = (abs(this%shift_envelope%value_at(x(ii)))*2.0_fp)**2*(1-this%g**2)/(4.0_fp*this%delta)
       !r(ii) = (abs(this%shift_envelope%value_at(x(ii)))*2.0_fp)**2*(-1*this%g**2)/(4.0_fp*(this%eps-this%delta)) +&
       !     (abs(this%shift_envelope%value_at(x(ii)))*2.0_fp)**2*(+1)/(4.0_fp*this%delta)

       r(ii) = abs(this%shift_envelope%value_at(x(ii)))**2*int_r

    end do
    phi = trapz(x,r)
    phi = phi*0.5_fp ! Additional 1/2 for two-photon pulse
    !write(*,*) phi
    !write(*,*) r
    !write(*,*) nis,x(nis)
    !open(11,file="phi_test.dat",position="append")
    !write(11,*) phi
    !close(11)
    !y = this%target_envelope%value_at(t)*exp(-(0.0_fp,1.0_fp)*phi)
    !fixed
    y = this%target_envelope%value_at(t)*exp((0.0_fp,1.0_fp)*phi)
  end function value_at_general_dfa
  !---------------------------------------------------------
  subroutine set_delta_general_dfa(this,delta)
    use utils,only:fp,ip
    implicit none
    class(general_dfa), intent(inout) :: this
    real(fp), intent(in) :: delta
    this%delta = delta
  end subroutine set_delta_general_dfa
  !---------------------------------------------------------
  subroutine initialize_general_dfa(this,shift_envelope,target_envelope,eps,delta,tis,dt,g)
    use utils,only:fp,ip
    class(general_dfa) :: this
    class(envelope),intent(in),target :: target_envelope,shift_envelope
    real(fp), intent(in) :: eps,delta,tis,dt,g

    this%target_envelope => target_envelope
    this%shift_envelope => shift_envelope
    this%eps = eps
    this%delta = delta
    this%g = g
    this%tis = tis
    this%dt = dt
  end subroutine initialize_general_dfa
  !---------------------------------------------------------
  subroutine set_parameters_general_dfa(this,parameter_array)
    use utils,only:fp,ip
    class(general_dfa), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    call this%target_envelope%set_parameters(parameter_array)
  end subroutine set_parameters_general_dfa
  !#########################################################
  function value_at_general_dfa_01(this,t) result(y)
    use utils,only:fp,ip
    use constants,only:pi
    use linfuncs,only:trapz
    implicit none
    class(general_dfa_01), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    real(fp) :: Omega,dt_c,phi
    integer(ip) :: ii,nis
    real(fp),dimension(:),allocatable :: x,r
    real(fp) :: int_r

    nis = ceiling((t - this%tis)/this%dt)+1
    dt_c = (t - this%tis)/(nis-1)
    allocate(x(nis),r(nis))

    int_r = 4.0_fp*((1*this%g**2)/(4.0_fp*(this%eps-this%delta)) + (+2)/(4.0_fp*this%delta))

    do ii=1,nis
       x(ii) = this%tis + (ii-1)*dt_c
       ! Corrections from the first oder pertubation theory. Assumes \delta = 2 \Delta
       !r(ii) = (abs(this%shift_envelope%value_at(x(ii)))*2.0_fp)**2*(this%g**2+2)/(4.0_fp*this%delta)
       !r(ii) = (abs(this%shift_envelope%value_at(x(ii)))*2.0_fp)**2*(+1*this%g**2)/(4.0_fp*(this%eps-this%delta)) +&
       !     (abs(this%shift_envelope%value_at(x(ii)))*2.0_fp)**2*(+2)/(4.0_fp*this%delta)
       r(ii) = abs(this%shift_envelope%value_at(x(ii)))**2*int_r

    end do
    phi = trapz(x,r)
    !phi = -0.5*ph
    !phi = 0.5*phi
    !write(*,*) x
    !write(*,*) r
    !write(*,*) nis,x(nis)
    !open(11,file="phi_test.dat",position="append")
    !write(11,*) phi
    !close(11)
    !y = this%target_envelope%value_at(t)*exp(-(0.0_fp,1.0_fp)*phi)
    !fixed
    y = this%target_envelope%value_at(t)*exp((0.0_fp,1.0_fp)*phi)
  end function value_at_general_dfa_01
  !---------------------------------------------------------
  subroutine set_delta_general_dfa_01(this,delta)
    use utils,only:fp,ip
    implicit none
    class(general_dfa_01), intent(inout) :: this
    real(fp), intent(in) :: delta
    this%delta = delta
  end subroutine set_delta_general_dfa_01
  !---------------------------------------------------------
  subroutine initialize_general_dfa_01(this,shift_envelope,target_envelope,eps,delta,tis,dt,g)
    use utils,only:fp,ip
    class(general_dfa_01) :: this
    class(envelope),intent(in),target :: target_envelope,shift_envelope
    real(fp), intent(in) :: eps,delta,tis,dt,g

    this%target_envelope => target_envelope
    this%shift_envelope => shift_envelope
    this%eps = eps
    this%delta = delta
    this%g = g
    this%tis = tis
    this%dt = dt
  end subroutine initialize_general_dfa_01
  !---------------------------------------------------------
  subroutine set_parameters_general_dfa_01(this,parameter_array)
    use utils,only:fp,ip
    class(general_dfa_01), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    call this%target_envelope%set_parameters(parameter_array)
  end subroutine set_parameters_general_dfa_01
  !#########################################################
  function value_at_general_dfa_12(this,t) result(y)
    use utils,only:fp,ip
    use constants,only:pi
    use linfuncs,only:trapz
    implicit none
    class(general_dfa_12), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    real(fp) :: Omega,dt_c,phi
    integer(ip) :: ii,nis
    real(fp),dimension(:),allocatable :: x,r
    real(fp) :: int_r

    nis = ceiling((t - this%tis)/this%dt)+1
    dt_c = (t - this%tis)/(nis-1)
    allocate(x(nis),r(nis))

    int_r = 4.0_fp*((-2*this%g**2)/(4.0_fp*(this%eps-this%delta)) + (-1)/(4.0_fp*this%delta))


    do ii=1,nis
       x(ii) = this%tis + (ii-1)*dt_c
       !r(ii) = (abs(this%shift_envelope%value_at(x(ii)))*2.0_fp)**2*(-2*this%g**2-1)/(4.0_fp*this%delta)
       !r(ii) = (abs(this%shift_envelope%value_at(x(ii)))*2.0_fp)**2*(-2*this%g**2)/(4.0_fp*(this%eps-this%delta)) +&
       !    (abs(this%shift_envelope%value_at(x(ii)))*2.0_fp)**2*(-1)/(4.0_fp*this%delta)
       r(ii) = abs(this%shift_envelope%value_at(x(ii)))**2*int_r

    end do
    phi = trapz(x,r)
    !phi = -0.5*phi
    !write(*,*) phi
    !write(*,*) r
    !write(*,*) nis,x(nis)
    !open(11,file="phi_test.dat",position="append")
    !write(11,*) phi
    !close(11)
    !y = this%target_envelope%value_at(t)*exp(-(0.0_fp,1.0_fp)*phi)
    !fixed
    y = this%target_envelope%value_at(t)*exp((0.0_fp,1.0_fp)*phi)
  end function value_at_general_dfa_12
  !---------------------------------------------------------
  subroutine set_delta_general_dfa_12(this,delta)
    use utils,only:fp,ip
    implicit none
    class(general_dfa_12), intent(inout) :: this
    real(fp), intent(in) :: delta
    this%delta = delta
  end subroutine set_delta_general_dfa_12
  !---------------------------------------------------------
  subroutine initialize_general_dfa_12(this,shift_envelope,target_envelope,eps,delta,tis,dt,g)
    use utils,only:fp,ip
    class(general_dfa_12) :: this
    class(envelope),intent(in),target :: target_envelope,shift_envelope
    real(fp), intent(in) :: eps,delta,tis,dt,g

    this%target_envelope => target_envelope
    this%shift_envelope => shift_envelope
    this%eps = eps
    this%delta = delta
    this%g = g
    this%tis = tis
    this%dt = dt
  end subroutine initialize_general_dfa_12
  !---------------------------------------------------------
  subroutine set_parameters_general_dfa_12(this,parameter_array)
    use utils,only:fp,ip
    class(general_dfa_12), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    call this%target_envelope%set_parameters(parameter_array)
  end subroutine set_parameters_general_dfa_12
  !#########################################################
  function value_at_flat_gaussian_envelope(this,t) result(y)
    use utils,only:fp,ip
    use constants,only:pi
    implicit none
    class(flat_gaussian_envelope), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    y = saturation_curve(real(this%A*exp(-(t-this%t0)**2/(2*this%sigma**2))),&
         this%saturation_start,this%saturation_end)
    !write(*,*) y
  end function value_at_flat_gaussian_envelope
  !---------------------------------------------------------
  function saturation_curve(current_value,saturation_start,saturation_end) result(y)
    implicit none
    real(fp), intent(in) :: current_value,saturation_start,saturation_end
    real(fp) :: y
    if(current_value>saturation_start) then
       y = saturation_start + (saturation_end-saturation_start)&
            *(1-exp(-(current_value-saturation_start)/(saturation_end-saturation_start)))
    else
       y = current_value
    end if
  end function saturation_curve
  !---------------------------------------------------------
  subroutine initialize_flat_gaussian_envelope(this,A,sigma,t0)
    use utils,only:fp,ip
    use constants, only:pi
    class(flat_gaussian_envelope) :: this
    complex(fp), intent(in) :: A
    real(fp), intent(in) :: sigma,t0
    this%A = A
    this%sigma = sigma
    this%t0 = t0
    this%saturation_start = 0.015e9*(2*pi)
    this%saturation_end = 0.02e9*(2*pi)
  end subroutine initialize_flat_gaussian_envelope
  !---------------------------------------------------------
  subroutine set_parameters_flat_gaussian_envelope(this,parameter_array)
    use utils,only:fp,ip
    class(flat_gaussian_envelope), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) < 3) then
       write(*,*) 'Number of parameters incorrect in set_parameters_gaussian_envelope.'
       call exit(1)
    end if
    this%A = parameter_array(1)
    this%sigma = real(parameter_array(2))
    this%t0 = real(parameter_array(3))
    !this%saturation_start = real(parameter_array(4))
    !this%saturation_end = real(parameter_array(5))
  end subroutine set_parameters_flat_gaussian_envelope
  !#########################################################
  function value_at_linear_zero_area_pulse(this,t) result(y)
    use utils,only:fp,ip
    class(linear_zero_area_pulse), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    if(t>this%t0 .and.  t<this%t0+this%tf) then

       if(t<this%t0+this%tf/4.0_fp) then
          y = this%A*(t-this%t0)/this%tf*4.0_fp
       elseif(t <this%t0 + this%tf*0.75_fp) then
          y = this%A*(1 - (t-this%t0 - this%tf*0.25_fp)/this%tf*4.0_fp)
       elseif(t  < this%t0 + this%tf) then
          y = this%A*(-1 + (t -this%t0 - this%tf*0.75_fp)/this%tf*4.0_fp)
       endif
    else
       y = 0.0_fp
    endif
  end function value_at_linear_zero_area_pulse
  !---------------------------------------------------------
  subroutine initialize_linear_zero_area_pulse(this,A,t0,tf)
    use utils,only:fp,ip
    class(linear_zero_area_pulse) :: this
    complex(fp), intent(in) :: A
    real(fp), intent(in) :: t0,tf
    this%A = A
    this%t0 = t0
    this%tf = tf
  end subroutine initialize_linear_zero_area_pulse
  !---------------------------------------------------------
  subroutine set_parameters_linear_zero_area_pulse(this,parameter_array)
    use utils,only:fp,ip
    class(linear_zero_area_pulse), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) < 1) then
       write(*,*) 'Number of parameters incorrect in set_parameters_constant_envelope.'
       call exit(1)
    else
       this%A = parameter_array(1)
       if(size(parameter_array) > 1) then
          this%t0 = real(parameter_array(2),kind=fp)
       endif
       if(size(parameter_array) > 2) then
          this%tf = real(parameter_array(3),kind=fp)
       endif
    end if
  end subroutine set_parameters_linear_zero_area_pulse
  !#########################################################
  function value_at_two_photon_linear_zero_area_pulse(this,t) result(y)
    use utils,only:fp,ip
    class(two_photon_linear_zero_area_pulse), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    if(t>this%t0 .and.  t<this%t0+this%tf) then

       if(t<this%t0+this%tf/4.0_fp) then
          y = this%A*sqrt((t-this%t0)/this%tf*4.0_fp)
       elseif(t <this%t0 + this%tf*0.5_fp) then
          y = this%A*sqrt((1 - (t-this%t0 - this%tf*0.25_fp)/this%tf*4.0_fp))
       elseif(t  < this%t0 + this%tf*0.75_fp) then
          y = this%A*exp(-cmplx(0.0_fp,this%phase_shift,kind=fp))*sqrt((t -this%t0 - this%tf*0.5_fp)/this%tf*4.0_fp)
       elseif(t  < this%t0 + this%tf) then
          y = this%A*exp(-cmplx(0.0_fp,this%phase_shift,kind=fp))*sqrt((1- (t -this%t0 - this%tf*0.75_fp)/this%tf*4.0_fp))
       endif
    else
       y = 0.0_fp
    endif
  end function value_at_two_photon_linear_zero_area_pulse
  !---------------------------------------------------------
  subroutine initialize_two_photon_linear_zero_area_pulse(this,A,t0,tf,phase_shift)
    use utils,only:fp,ip
    class(two_photon_linear_zero_area_pulse) :: this
    complex(fp), intent(in) :: A
    real(fp), intent(in) :: t0,tf,phase_shift
    this%A = A
    this%t0 = t0
    this%tf = tf
    this%phase_shift = phase_shift
  end subroutine initialize_two_photon_linear_zero_area_pulse
  !---------------------------------------------------------
  subroutine set_parameters_two_photon_linear_zero_area_pulse(this,parameter_array)
    use utils,only:fp,ip
    class(two_photon_linear_zero_area_pulse), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) < 1) then
       write(*,*) 'Number of parameters incorrect in set_parameters_constant_envelope.'
       call exit(1)
    else
       this%A = parameter_array(1)
       if(size(parameter_array) > 1) then
          this%t0 = real(parameter_array(2),kind=fp)
       endif
       if(size(parameter_array) > 2) then
          this%tf = real(parameter_array(3),kind=fp)
       endif
       if(size(parameter_array) > 3) then
          this%phase_shift = real(parameter_array(4),kind=fp)
       endif
    end if
  end subroutine set_parameters_two_photon_linear_zero_area_pulse
  !#########################################################
  function value_at_two_photon_linear_trnc_pulse(this,t) result(y)
    use utils,only:fp,ip
    class(two_photon_linear_trnc_pulse), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    real(fp) :: t1,t2,t0
    t1 = this%tf/2.0_fp + this%t0
    t2 = this%t0 + this%tf
    t0 = this%t0

    if(t>this%t0 .and.  t<this%t0+this%tf) then

       if(t<t1) then
          y = this%A*sqrt(1 - (t - t0)/(t1-t0))
       elseif(t<0.5_fp*(t1+t2)) then
          y = this%A2*exp(-cmplx(0.0_fp,this%phase_shift,kind=fp))*sqrt((t - t1)/(t2 - 0.5_fp*(t1+t2)))
       else
          y = this%A2*exp(-cmplx(0.0_fp,this%phase_shift,kind=fp))*sqrt(1 - (t - (t1+t2)/2.0_fp)/(t2 - 0.5_fp*(t2+t1)))
       endif
    else
       y = 0.0_fp
    endif
  end function value_at_two_photon_linear_trnc_pulse
  !---------------------------------------------------------
  subroutine initialize_two_photon_linear_trnc_pulse(this,A,t0,tf,phase_shift,A2)
    use utils,only:fp,ip
    class(two_photon_linear_trnc_pulse) :: this
    complex(fp), intent(in) :: A
    real(fp), intent(in) :: t0,tf,phase_shift
    complex(fp),intent(in),optional :: A2
    this%A = A
    this%t0 = t0
    this%tf = tf
    this%phase_shift = phase_shift
    if(present(A2)) then
       this%A2 = A2
    else
       this%A2 = A
    endif
  end subroutine initialize_two_photon_linear_trnc_pulse
  !---------------------------------------------------------
  subroutine set_parameters_two_photon_linear_trnc_pulse(this,parameter_array)
    use utils,only:fp,ip
    class(two_photon_linear_trnc_pulse), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) < 1) then
       write(*,*) 'Number of parameters incorrect in set_parameters_constant_envelope.'
       call exit(1)
    else
       this%A = parameter_array(1)
       if(size(parameter_array) > 1) then
          this%t0 = real(parameter_array(2),kind=fp)
       endif
       if(size(parameter_array) > 2) then
          this%tf = real(parameter_array(3),kind=fp)
       endif
       if(size(parameter_array) > 3) then
          this%phase_shift = real(parameter_array(4),kind=fp)
       endif
       if(size(parameter_array) > 4) then
          this%A2 = real(parameter_array(5),kind=fp)
       endif
    end if
  end subroutine set_parameters_two_photon_linear_trnc_pulse
  !#########################################################
  function value_at_exp_polinom(this,t) result(y)
    use utils,only:fp,ip
    class(exp_polinom), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    !write(*,*) this%t0,this%n,this%td,this%cutoff,this%trnc_point,this%A
    if(t <= this%trnc_point .and. t >= this%t0-this%cutoff*this%td &
         .and. t <= this%t0+this%cutoff*this%td) then
       y = this%A*exp(-(t - this%t0)**this%n/(2*this%td**this%n))
    else
       y = 0.0_fp
    end if
    !write(*,*) y
  end function value_at_exp_polinom
  !---------------------------------------------------------
  subroutine initialize_exp_polinom(this,A,td,t0,n,cutoff,trnc_point)
    use utils,only:fp,ip
    class(exp_polinom) :: this
    complex(fp), intent(in) :: A
    real(fp), intent(in) :: t0,td,n
    real(fp), intent(in), optional :: cutoff,trnc_point
    this%A = A
    this%t0 = t0
    this%td = td
    this%n = n
    if(present(cutoff)) then
       this%cutoff = cutoff
    else
       this%cutoff = huge(1.0_fp)
    end if
    if(present(trnc_point)) then
       this%trnc_point = trnc_point
    else
       this%trnc_point = huge(1.0_fp)
    endif
  end subroutine initialize_exp_polinom
  !---------------------------------------------------------
  subroutine set_parameters_exp_polinom(this,parameter_array)
    use utils,only:fp,ip
    class(exp_polinom), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) < 1) then
       write(*,*) 'Number of parameters incorrect.'
       call exit(1)
    else
       this%A = parameter_array(1)
       if(size(parameter_array) > 1) then
          this%td = real(parameter_array(2),kind=fp)
       endif
       if(size(parameter_array) > 2) then
          this%t0 = real(parameter_array(3),kind=fp)
       endif
       ! if(size(parameter_array) > 3) then
       !    this%n = real(parameter_array(4),kind=fp)
       ! endif
       if(size(parameter_array) > 3) then
          this%cutoff = real(parameter_array(4),kind=fp)
       endif
       if(size(parameter_array) > 4) then
          this%trnc_point = real(parameter_array(5),kind=fp)
       endif
    end if
  end subroutine set_parameters_exp_polinom
  !#########################################################
  function value_at_polynomial_envelope(this,t) result(y)
    use utils,only:fp,ip
    class(polynomial_envelope), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    integer(ip) :: ii,ll

    ll = size(this%A)
    y = 0.0_fp
    if(t>this%trnc_point) then
       y = 0.0_fp
    else
       do ii=1,ll
          if(t>0.0_fp) then
             y = y + this%A(ii)*t**(ii-1)
          end if
       end do
    end if
  end function value_at_polynomial_envelope
  !---------------------------------------------------------
  subroutine initialize_polynomial_envelope(this,A,trnc_point)
    use utils,only:fp,ip
    class(polynomial_envelope) :: this
    complex(fp),dimension(:), intent(in) :: A
    real(fp), intent(in),optional ::trnc_point
    allocate(this%A(size(A)))
    this%A = A
    if(present(trnc_point)) then
       this%trnc_point = trnc_point
    else
       this%trnc_point = huge(1.0_fp)
    endif
  end subroutine initialize_polynomial_envelope
  !---------------------------------------------------------
  subroutine set_parameters_polynomial_envelope(this,parameter_array)
    use utils,only:fp,ip
    class(polynomial_envelope), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) < 1) then
       write(*,*) 'Number of parameters incorrect.'
       call exit(1)
    else
       if(size(parameter_array) > size(this%A)) then
          this%A = real(parameter_array(1:size(this%A)))
       else
          this%A(1:size(parameter_array)) = parameter_array
       end if
    end if
  end subroutine set_parameters_polynomial_envelope
  !#########################################################
  function value_at_sqrt_polynomial_envelope(this,t) result(y)
    use utils,only:fp,ip
    class(sqrt_polynomial_envelope), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    integer(ip) :: ii,ll

    ll = size(this%A)
    y = 0.0_fp
    if(t>this%trnc_point) then
       y = 0.0_fp
    else
       do ii=1,ll
          if(t>0.0_fp) then
             y = y + this%A(ii)*t**(ii-1)
          end if
       end do
    end if
    y = sqrt(y)
  end function value_at_sqrt_polynomial_envelope
  !---------------------------------------------------------
  subroutine initialize_sqrt_polynomial_envelope(this,A,trnc_point)
    use utils,only:fp,ip
    class(sqrt_polynomial_envelope) :: this
    complex(fp),dimension(:), intent(in) :: A
    real(fp), intent(in),optional ::trnc_point
    allocate(this%A(size(A)))
    this%A = A
    if(present(trnc_point)) then
       this%trnc_point = trnc_point
    else
       this%trnc_point = huge(1.0_fp)
    endif
  end subroutine initialize_sqrt_polynomial_envelope
  !---------------------------------------------------------
  subroutine set_parameters_sqrt_polynomial_envelope(this,parameter_array)
    use utils,only:fp,ip
    class(sqrt_polynomial_envelope), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) < 1) then
       write(*,*) 'Number of parameters incorrect.'
       call exit(1)
    else
       if(size(parameter_array) > size(this%A)) then
          this%A = real(parameter_array(1:size(this%A)))
       else
          this%A(1:size(parameter_array)) = parameter_array
       end if
    end if
  end subroutine set_parameters_sqrt_polynomial_envelope
  !#########################################################
  function value_at_lr_envelope_01(this,t) result(y)
    use utils,only:fp,ip
    class(lr_envelope_01), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    integer(ip) :: ii,ll

    if(t>this%trnc_point) then
       y = 0.0_fp
    else
       y = this%factor*(this%dalpha%value_at(t)*sin(this%alpha%value_at(t))*&
            tan(this%beta%value_at(t)) - this%dbeta%value_at(t)*&
            cos(this%alpha%value_at(t)) - exp((0.0,1.0)*PI/2)*this%Omega02%value_at(t)*&
            sin(this%alpha%value_at(t))*tan(this%beta%value_at(t)))
    end if
  end function value_at_lr_envelope_01
  !---------------------------------------------------------
  subroutine initialize_lr_envelope_01(this,alpha,beta,dalpha,dbeta,Omega02,factor,trnc_point)
    use utils,only:fp,ip
    class(lr_envelope_01) :: this
    class(envelope),target :: alpha,beta,dalpha,dbeta,Omega02
    complex(fp),intent(in) :: factor
    real(fp),intent(in),optional :: trnc_point
    this%alpha => alpha
    this%beta => beta
    this%dbeta => dbeta
    this%dalpha => dalpha
    this%Omega02 => Omega02
    this%factor = factor
    if(present(trnc_point)) then
       this%trnc_point = trnc_point
    else
       this%trnc_point = huge(1.0_fp)
    end if
  end subroutine initialize_lr_envelope_01
  !---------------------------------------------------------
  subroutine set_parameters_lr_envelope_01(this,parameter_array)
    use utils,only:fp,ip
    class(lr_envelope_01), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) > 0) then
       this%factor = parameter_array(1)
    end if
    if(size(parameter_array) > 1) then
       this%trnc_point = real(parameter_array(2))
    end if
  end subroutine set_parameters_lr_envelope_01
  !#########################################################
  function value_at_lr_envelope_12(this,t) result(y)
    use utils,only:fp,ip
    class(lr_envelope_12), intent(in) :: this
    real(fp), intent(in) :: t
    complex(fp) :: y
    integer(ip) :: ii,ll
    if(t > this%trnc_point) then
       y = 0.0_fp
    else
       y = this%factor*(this%dalpha%value_at(t)*cos(this%alpha%value_at(t))*&
            tan(this%beta%value_at(t)) + this%dbeta%value_at(t)*&
            sin(this%alpha%value_at(t)) - exp((0.0,1.0)*PI/2)*this%Omega02%value_at(t)*&
            cos(this%alpha%value_at(t))*tan(this%beta%value_at(t)))
    end if
  end function value_at_lr_envelope_12
  !---------------------------------------------------------
  subroutine initialize_lr_envelope_12(this,alpha,beta,dalpha,dbeta,Omega02,factor,trnc_point)
    use utils,only:fp,ip
    class(lr_envelope_12) :: this
    class(envelope),target :: alpha,beta,dalpha,dbeta,Omega02
    complex(fp),intent(in) :: factor
    real(fp),intent(in),optional :: trnc_point
    this%alpha => alpha
    this%beta => beta
    this%dbeta => dbeta
    this%dalpha => dalpha
    this%Omega02 => Omega02
    this%factor = factor
    if(present(trnc_point)) then
       this%trnc_point = trnc_point
    else
       this%trnc_point = huge(1.0_fp)
    end if
  end subroutine initialize_lr_envelope_12
  !---------------------------------------------------------
  subroutine set_parameters_lr_envelope_12(this,parameter_array)
    use utils,only:fp,ip
    class(lr_envelope_12), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    if(size(parameter_array) > 0) then
       this%factor = parameter_array(1)
    end if
    if(size(parameter_array) > 1) then
       this%trnc_point = real(parameter_array(2))
    end if
  end subroutine set_parameters_lr_envelope_12
  !#########################################################
  function value_at_sa_stirap_det_01(this,t) result(y)
    use utils,only:fp,ip
    use constants,only:pi
    use linfuncs,only:trapz
    implicit none
    class(sa_stirap_det_01), intent(in) :: this
    real(fp), intent(in) :: t
    real(fp) :: theta,dphi,phi0,phi1,t0,t1,Omega010,Omega011
    real(fp) :: Omega120,Omega121,value
    complex(fp) :: y
    real(fp) :: phase

    t0 = t-this%dt/2.0_fp
    t1 = t+this%dt/2.0_fp

    Omega010 = abs(this%e01%value_at(t0)*2)
    Omega011 = abs(this%e01%value_at(t1)*2)
    
    Omega120 = abs(this%e12%value_at(t0)*2)
    Omega121 = abs(this%e12%value_at(t1)*2)
    
    theta = atan(this%e01%value_at(t)/this%e12%value_at(t))
    phi0 = atan(sqrt(Omega010**2 + Omega120**2)/(sqrt(Omega010**2 + Omega120**2 + this%delta**2) + this%delta))
    phi1 = atan(sqrt(Omega011**2 + Omega121**2)/(sqrt(Omega011**2 + Omega121**2 + this%delta**2) + this%delta))
    dphi = (phi1-phi0)/this%dt

    y = this%e01%value_at(t)
    phase = atan2(aimag(y),real(y))
    if(this%e12%value_at(t) /= 0.0_fp) then
       y = y + exp((0.0_fp,1.0_fp)+(0.0_fp,1.0_fp)*phase)*sin(theta)*dphi
    endif
  end function value_at_sa_stirap_det_01
  !---------------------------------------------------------
  subroutine set_delta_sa_stirap_det_01(this,delta)
    use utils,only:fp,ip
    implicit none
    class(sa_stirap_det_01), intent(inout) :: this
    real(fp), intent(in) :: delta
    this%delta = delta
  end subroutine set_delta_sa_stirap_det_01
  !---------------------------------------------------------
  subroutine initialize_sa_stirap_det_01(this,e01,e12,delta,dt)
    use utils,only:fp,ip
    class(sa_stirap_det_01) :: this
    class(envelope),intent(in),target :: e01,e12
    real(fp), intent(in) :: delta,dt
    
    this%e01 => e01
    this%e12 => e12

    this%delta = delta
    this%dt = dt

  end subroutine initialize_sa_stirap_det_01
  !---------------------------------------------------------
  subroutine set_parameters_sa_stirap_det_01(this,parameter_array)
    use utils,only:fp,ip
    class(sa_stirap_det_01), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    call this%e01%set_parameters(parameter_array)
    !write(*,*) parameter_array
  end subroutine set_parameters_sa_stirap_det_01
  !#########################################################
  function value_at_sa_stirap_det_02(this,t) result(y)
    use utils,only:fp,ip
    use constants,only:pi
    use linfuncs,only:trapz
    implicit none
    class(sa_stirap_det_02), intent(in) :: this
    real(fp), intent(in) :: t
    real(fp) :: dphi,phi0,phi1,t0,t1
    complex(fp) :: Omega120,Omega121,value,Omega010,Omega011
    complex(fp) :: y
    real(fp) :: phase01,phase12
    real(fp) :: theta0,theta1,dtheta
    integer(ip) :: ii
    class(envelope),pointer :: env

    t0 = t-this%dt/2.0_fp
    t1 = t+this%dt/2.0_fp

    Omega011 = (0.0_fp,0.0_fp)
    Omega121 = (0.0_fp,0.0_fp)
    Omega010 = (0.0_fp,0.0_fp)
    Omega120 = (0.0_fp,0.0_fp)



    ! In this code Omega means the raw coupling, not Rabi rate
    do ii=1,this%nenv01
       env => this%d01(ii)%obj
       Omega010 = Omega010 + env%value_at(t0)
       Omega011 = Omega011 + env%value_at(t1)
    end do
    do ii=1,this%nenv12
       env => this%d12(ii)%obj
       Omega120 = Omega120 + env%value_at(t0)
       Omega121 = Omega121 + env%value_at(t1)
    end do

    ! if(abs(Omega120) .eq. 0.0_fp .and. abs(Omega010) .eq. 0.0_fp) then
    !    theta0 = 0.0_fp
    ! else
    !    theta0 = atan2(abs(Omega120),abs(Omega010))
    ! end if

    !     if(abs(Omega121) .eq. 0.0_fp .and. abs(Omega011) .eq. 0.0_fp) then
    !    theta1 = 0.0_fp
    ! else
    !    theta1 = atan2(abs(Omega121),abs(Omega011))
    ! end if

    if((abs(Omega120) > tiny(0.0_fp) .and. abs(Omega010) > tiny(0.0_fp)) .and. &
         (abs(Omega121) > tiny(0.0_fp) .and. abs(Omega011) > tiny(0.0_fp))) then
       theta0 = atan2(abs(Omega120)*this%csr,abs(Omega010))
       theta1 = atan2(abs(Omega121)*this%csr,abs(Omega011))
       dtheta = (theta1-theta0)/this%dt
    else
       dtheta = 0.0_fp
    endif

    if(abs(Omega010) .eq. 0.0_fp) then
       phase01 = 0.0_fp
    else
       phase01 = atan2(aimag(Omega010),real(Omega010))
    endif

    if(abs(Omega120) .eq. 0.0_fp) then
       phase12 = 0.0_fp
    else
       phase12 = atan2(aimag(Omega120),real(Omega120))
    endif
    
    dtheta = (theta1-theta0)/this%dt
    y = this%A_factor*sqrt(this%delta*abs(dtheta))*exp((0.0_fp,1.0_fp)*(this%phi+0.5_fp*(phase01+phase12)))
    !y = sqrt(this%delta*abs(dtheta))*exp((0.0_fp,1.0_fp)*this%phi)
     ! if(abs(y) > 0.0_fp) then
     !    write(*,*) dtheta,abs(y)
     ! end if
  end function value_at_sa_stirap_det_02
  !---------------------------------------------------------
  subroutine add_envelope_01_sa_stirap_det_02(this,envelope_in)
    use utils,only:fp,ip
    class(sa_stirap_det_02),intent(inout) :: this
    class(envelope),intent(in),target :: envelope_in
    type(envelope_container),allocatable,dimension(:) :: temp_cont
    if(this%nenv01 >= this%env_cont_size_01) then
       allocate(temp_cont(this%nenv01))
       temp_cont = this%d01(1:this%nenv01)
       deallocate(this%d01)
       allocate(this%d01(this%nenv01*2))
       this%d01(1:this%nenv01) = temp_cont
       deallocate(temp_cont)
       this%env_cont_size_01 = size(this%d01)
       !write(*,*) 'Extending space for envelopes. New size is', this%env_cont_size,'.'
    end if
    this%nenv01 = this%nenv01 + 1
    this%d01(this%nenv01)%obj => envelope_in
  end subroutine add_envelope_01_sa_stirap_det_02
  !---------------------------------------------------------
  subroutine add_envelope_12_sa_stirap_det_02(this,envelope_in)
    use utils,only:fp,ip
    class(sa_stirap_det_02),intent(inout) :: this
    class(envelope),intent(in),target :: envelope_in
    type(envelope_container),allocatable,dimension(:) :: temp_cont
    if(this%nenv12 >= this%env_cont_size_12) then
       allocate(temp_cont(this%nenv12))
       temp_cont = this%d12(1:this%nenv12)
       deallocate(this%d12)
       allocate(this%d12(this%nenv12*2))
       this%d12(1:this%nenv12) = temp_cont
       deallocate(temp_cont)
       this%env_cont_size_12 = size(this%d12)
       !write(*,*) 'Extending space for envelopes. New size is', this%env_cont_size,'.'
    end if
    this%nenv12 = this%nenv12 + 1
    this%d12(this%nenv12)%obj => envelope_in
  end subroutine add_envelope_12_sa_stirap_det_02
  !---------------------------------------------------------
  subroutine initialize_sa_stirap_det_02(this,delta,dt,phi,A_factor,csr)
    use utils,only:fp,ip
    class(sa_stirap_det_02) :: this
    !class(signal),intent(in),target :: d01,d12
    !type(envelope_container),dimension(:),pointer :: d01,d12
    real(fp), intent(in) :: dt,delta,phi,A_factor,csr
    integer(ip) :: initial_env_cont_size
    !this%d01 => d01
    !this%d12 => d12
    this%dt = dt
    this%delta = delta
    this%phi = phi
    this%A_factor = A_factor
    this%csr = csr
    !this%nenv01 => nenv01
    !this%nenv12 => nenv12

    initial_env_cont_size = 2
    this%nenv01 = 0
    allocate(this%d01(initial_env_cont_size))
    this%env_cont_size_01 = initial_env_cont_size
    this%nenv12 = 0
    allocate(this%d12(initial_env_cont_size))
    this%env_cont_size_12 = initial_env_cont_size
  end subroutine initialize_sa_stirap_det_02
  !---------------------------------------------------------
  subroutine set_parameters_sa_stirap_det_02(this,parameter_array)
    use utils,only:fp,ip
    class(sa_stirap_det_02), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    
    if(size(parameter_array) < 4) then
       write(*,*) 'Number of parameters incorrect in set_parameters_constant_sa_stirap_det_02.'
       call exit(1)
    end if
    this%delta = real(parameter_array(1))
    this%dt = real(parameter_array(2))
    this%phi = real(parameter_array(3))
    this%A_factor = real(parameter_array(4))
  end subroutine set_parameters_sa_stirap_det_02
  !#########################################################
  function value_at_sa_stirap_det_12(this,t) result(y)
    use utils,only:fp,ip
    use constants,only:pi
    use linfuncs,only:trapz
    implicit none
    class(sa_stirap_det_12), intent(in) :: this
    real(fp), intent(in) :: t
    real(fp) :: theta,dphi,phi0,phi1,t0,t1,Omega010,Omega011
    real(fp) :: Omega120,Omega121,value
    complex(fp) :: y
    real(fp) :: phase

    t0 = t-this%dt/2.0_fp
    t1 = t+this%dt/2.0_fp

    Omega010 = abs(this%e01%value_at(t0)*2)
    Omega011 = abs(this%e01%value_at(t1)*2)
    
    Omega120 = abs(this%e12%value_at(t0)*2)
    Omega121 = abs(this%e12%value_at(t1)*2)
    
    theta = atan(this%e01%value_at(t)/this%e12%value_at(t))
    phi0 = atan(sqrt(Omega010**2 + Omega120**2)/(sqrt(Omega010**2 + Omega120**2 + this%delta**2) + this%delta))
    phi1 = atan(sqrt(Omega011**2 + Omega121**2)/(sqrt(Omega011**2 + Omega121**2 + this%delta**2) + this%delta))
    dphi = (phi1-phi0)/this%dt
    

    y = this%e12%value_at(t)
    phase = atan2(aimag(y),real(y))

    if(this%e12%value_at(t) /= 0.0_fp) then

       y = y + exp((0.0_fp,1.0_fp)+(0.0_fp,1.0_fp)*phase)*(-cos(theta)*dphi)
    end if
  end function value_at_sa_stirap_det_12
  !---------------------------------------------------------
  subroutine set_delta_sa_stirap_det_12(this,delta)
    use utils,only:fp,ip
    implicit none
    class(sa_stirap_det_12), intent(inout) :: this
    real(fp), intent(in) :: delta
    this%delta = delta
  end subroutine set_delta_sa_stirap_det_12
  !---------------------------------------------------------
  subroutine initialize_sa_stirap_det_12(this,e01,e12,delta,dt)
    use utils,only:fp,ip
    class(sa_stirap_det_12) :: this
    class(envelope),intent(in),target :: e01,e12
    real(fp), intent(in) :: delta,dt
    
    this%e01 => e01
    this%e12 => e12

    this%delta = delta
    this%dt = dt
  end subroutine initialize_sa_stirap_det_12
  !---------------------------------------------------------
  subroutine set_parameters_sa_stirap_det_12(this,parameter_array)
    use utils,only:fp,ip
    class(sa_stirap_det_12), intent(inout) :: this
    complex(fp),dimension(:),intent(in) :: parameter_array
    call this%e12%set_parameters(parameter_array)
    !write(*,*) parameter_array
  end subroutine set_parameters_sa_stirap_det_12
  !#########################################################
end module envelopes
