module linfuncs
use utils, only: fp,ip

interface diagonal
   module procedure diagonal_real,diagonal_complex
end interface diagonal

interface get_diag_matrix
   module procedure get_diag_matrix_real,get_diag_matrix_complex
end interface get_diag_matrix

interface kron
   module procedure kron_real,kron_complex,kron_real_complex,kron_complex_real
end interface kron

interface print_matrix
   module procedure print_matrix_real,print_matrix_complex
end interface print_matrix

interface trace
   module procedure trace_real,trace_complex
end interface trace

interface trace_abs
   module procedure trace_abs_real,trace_abs_complex
end interface trace_abs


interface l2norm
   module procedure l2norm_real,l2norm_complex
end interface l2norm

interface interp1
   module procedure interp1_real,interp1_complex
end interface interp1

contains
subroutine linspace(a, b, n,out)
  implicit none
  real(fp) ::a,b
  integer(ip) :: n,ii
  real(fp), dimension(n) :: out
  real(fp), dimension(n) :: index_vector
  
  do ii = 1,n
     index_vector(ii) = ii-1
  end do
  if(n>1) then
     out = index_vector / (n - 1) * (b - a) + a
     out(n) = b ! Just to make sure that the last index really is correct.
  else
     out(1) = a
  end if

end subroutine linspace
!#############################################
subroutine print_matrix_complex(mat,is_real_in)
implicit none
integer(ip) :: dim
complex(fp),dimension(:,:),intent(in) :: mat
integer(ip) :: ii,jj
logical,optional,intent(in) :: is_real_in
logical :: is_real

if(present(is_real_in)) then
   is_real = is_real_in
else
   is_real = .false.
end if

dim = size(mat,1)
if (dim > 60) then
   write(*,*) 'unable to print larger than 60x60 matrices.'
end if
write(*,*) 'entering. dim is ', dim
do ii=1,dim
   if(is_real) then
      write(*,'(601pe10.3)') (real(mat(ii,jj)),jj=1,dim)
   else
      write(*,'(60(1pe10.3,e10.3))') (real(mat(ii,jj)),aimag(mat(ii,jj)),jj=1,dim)
   end if
end do
end subroutine print_matrix_complex
!-----------------------------------------------------------
subroutine print_matrix_real(mat,is_real_in)
implicit none
integer(ip) :: dim
real(fp),dimension(:,:),intent(in) :: mat
integer(ip) :: ii,jj
logical,intent(in),optional :: is_real_in

dim = size(mat,1)
if (dim > 60) then
   write(*,*) 'unable to print larger than 60x60 matrices.'
end if
write(*,*) 'entering. dim is ', dim
do ii=1,dim
   write(*,'(60e10.3)') (real(mat(ii,jj)),jj=1,dim)
end do
end subroutine print_matrix_real
!#############################################
subroutine print_real_matrix(mat,dim)
implicit none
integer(ip) :: dim
real(fp),dimension(dim,dim) :: mat
integer(ip) :: ii,jj

if (dim > 60) then
   write(*,*) 'unable to print larger than 60x60 matrices.'
end if
write(*,*) 'entering. dim is ', dim
do ii=1,dim
   write(*,'(60e10.3)') (mat(ii,jj),jj=1,dim)
end do
end subroutine print_real_matrix
!#############################################
function kron_complex(a,b) result(k)
    implicit none
    ! k = a\cross b
    complex(fp), dimension(:,:),intent(in) :: a
    complex(fp), dimension(:,:),intent(in) :: b

    integer(ip) :: i, j

    integer(ip) :: ma
    integer(ip) :: na
    integer(ip) :: mb
    integer(ip) :: nb
    complex(fp), dimension(size(a,1)*size(b,1),size(a,2)*size(b,2)) :: k
    ma = size(a,1)
    na = size(a,2)
    mb = size(b,1)
    nb = size(b,2)

    !if (size(k,1) /= ma*mb .or. size(k,2) /= na*nb) then
    !   write(*,*) 'k has invalid size'
       !call stop
    !end if
    k = (0.0_fp,0.0_fp)
    forall(i=1:ma, j=1:na)
        k(mb*(i-1)+1:mb*i,nb*(j-1)+1:nb*j) = a(i,j)*b
    end forall
  end function kron_complex
!-----------------------------------------------------------
function kron_real(a,b) result(k)
    implicit none
    ! k = a\cross b
    real(fp), dimension(:,:),intent(in) :: a
    real(fp), dimension(:,:),intent(in) :: b

    integer(ip) :: i, j

    integer(ip) :: ma
    integer(ip) :: na
    integer(ip) :: mb
    integer(ip) :: nb
    real(fp), dimension(size(a,1)*size(b,1),size(a,2)*size(b,2)) :: k
    ma = size(a,1)
    na = size(a,2)
    mb = size(b,1)
    nb = size(b,2)

    !if (size(k,1) /= ma*mb .or. size(k,2) /= na*nb) then
    !   write(*,*) 'k has invalid size'
       !call stop
    !end if
    k = 0.0_fp
    forall(i=1:ma, j=1:na)
        k(mb*(i-1)+1:mb*i,nb*(j-1)+1:nb*j) = a(i,j)*b
    end forall
  end function kron_real
!-----------------------------------------------------------
function kron_real_complex(a,b) result(k)
    implicit none
    ! k = a\cross b
    real(fp), dimension(:,:),intent(in) :: a
    complex(fp), dimension(:,:),intent(in) :: b

    integer(ip) :: i, j

    integer(ip) :: ma
    integer(ip) :: na
    integer(ip) :: mb
    integer(ip) :: nb
    complex(fp), dimension(size(a,1)*size(b,1),size(a,2)*size(b,2)) :: k
    ma = size(a,1)
    na = size(a,2)
    mb = size(b,1)
    nb = size(b,2)

    !if (size(k,1) /= ma*mb .or. size(k,2) /= na*nb) then
    !   write(*,*) 'k has invalid size'
       !call stop
    !end if
    k = (0.0_fp,0.0_fp)
    forall(i=1:ma, j=1:na)
        k(mb*(i-1)+1:mb*i,nb*(j-1)+1:nb*j) = a(i,j)*b
    end forall
  end function kron_real_complex
!-----------------------------------------------------------
function kron_complex_real(a,b) result(k)
    implicit none
    ! k = a\cross b
    complex(fp), dimension(:,:),intent(in) :: a
    real(fp), dimension(:,:),intent(in) :: b

    integer(ip) :: i, j

    integer(ip) :: ma
    integer(ip) :: na
    integer(ip) :: mb
    integer(ip) :: nb
    complex(fp), dimension(size(a,1)*size(b,1),size(a,2)*size(b,2)) :: k
    ma = size(a,1)
    na = size(a,2)
    mb = size(b,1)
    nb = size(b,2)

    !if (size(k,1) /= ma*mb .or. size(k,2) /= na*nb) then
    !   write(*,*) 'k has invalid size'
       !call stop
    !end if
    k = (0.0_fp,0.0_fp)
    forall(i=1:ma, j=1:na)
        k(mb*(i-1)+1:mb*i,nb*(j-1)+1:nb*j) = a(i,j)*b
    end forall
  end function kron_complex_real
!#############################################
subroutine expm(a,adim,out)
implicit none
integer(ip) :: adim,counter,max_num_terms,ii
complex(fp), dimension(adim,adim) :: a,out,apt
real(fp) :: tol,coef
tol = 1e-2
max_num_terms = 50
out = (0.0,0.0)
do ii=1,adim
   out(ii,ii) = (1.0,0.0)
end do

apt = out

coef = 1.0
counter = 0

do while (maxval(real(matmul(apt,conjg(transpose(apt)))))/coef > tol .and. counter < max_num_terms)
   counter = counter + 1
   apt = matmul(apt,a)
   coef = coef*counter
   out = out + apt/coef
 !  write(*,*) 'counter: ', counter
!   call print_matrix(out,adim,.true.)
end do

end subroutine expm
!#############################################
function diagonal_complex(a) result(diagonal)! COMPLEX
implicit none
complex(fp), dimension(:,:),intent(in) :: a
integer(ip) :: ii,adim
complex(fp), dimension(size(a,1)) :: diagonal
adim = size(a,1)
do ii=1,adim
   diagonal(ii) = a(ii,ii)
end do
end function diagonal_complex
!--------------------------------------------
function diagonal_real(a) result(diagonal)! REAL
implicit none
real(fp), dimension(:,:),intent(in) :: a
integer(ip) :: ii,adim
real(fp), dimension(size(a,1)) :: diagonal
adim = size(a,1)
do ii=1,adim
   diagonal(ii) = a(ii,ii)
end do
end function diagonal_real
!#############################################
function get_eye(dim) result(eye)
    implicit none
    integer(ip), intent(in) :: dim
    real(fp), dimension(dim,dim) :: eye
    integer(ip) :: ii
    eye = 0.0_fp
    do ii=1,dim
       eye(ii,ii) = 1.0_fp
    end do
  end function get_eye
!#############################################
  function get_diag_matrix_complex(diag) result(eye)
    implicit none
    complex(fp), dimension(:),intent(in) :: diag
    complex(fp), dimension(size(diag),size(diag)) :: eye
    integer(ip) :: ii
    eye = (0.0_fp,0.0_fp)
    do ii=1,size(diag)
       eye(ii,ii) = diag(ii)
    end do
  end function get_diag_matrix_complex
!-----------------------------------------------------------
  function get_diag_matrix_real(diag) result(eye)
    implicit none
    real(fp), dimension(:),intent(in) :: diag
    real(fp), dimension(size(diag),size(diag)) :: eye
    integer(ip) :: ii
    eye = 0.0_fp
    do ii=1,size(diag)
       eye(ii,ii) = diag(ii)
    end do
  end function get_diag_matrix_real
  !#############################################
  function get_creation_op(dim) result(cc)
    implicit none
    integer(ip), intent(in) :: dim
    complex(fp), dimension(dim,dim) :: cc
    integer(ip) :: ii
    
    cc = (0.0_fp,0.0_fp)
    do ii=1,dim-1
       cc(ii+1,ii) = cmplx(sqrt(real(ii,kind=fp)),kind=fp)
       !cc(ii+1,ii) = (0.0_fp,1.0_fp)
    end do
  end function get_creation_op
!#############################################
  function get_annihilation_op(dim) result(aa)
    implicit none
    integer(ip), intent(in) :: dim
    complex(fp), dimension(dim,dim) :: aa
    integer(ip) :: ii
    aa = (0.0_fp,0.0_fp)
    do ii=1,dim-1
       aa(ii,ii+1) = cmplx(sqrt(real(ii,kind=fp)),kind=fp)
    end do
  end function get_annihilation_op
!#############################################
  ! Performs cubic Hermite interpolation for complex numbers.
  function interp1_complex(x,y,xquery) result(yquery)
    implicit none
    real(fp),dimension(:),intent(in) :: x
    complex(fp), dimension(:,:),intent(in) :: y
    real(fp),dimension(:),intent(in) :: xquery
    complex(fp), dimension(size(y,1),size(xquery)) :: yquery
    integer(ip) :: dim_in,dim,ii,jj,n
    real(fp),dimension(2*size(y,1),size(x)) :: ya,fd
    real(fp) :: xk,xk1,t,h00,h01,h10,h11
    real(fp) :: yquery_t(2*size(y,1))
    dim_in = size(x)
    dim = size(xquery)
    n = size(y,1)
    ! First calculate the finite differences
    !write(*,*) 'xquery:', xquery
    !write(*,*) 'x', x
    if(x(dim_in) < xquery(dim) .or. xquery(1) < x(1)) then
       write(*,*) 'In interp1, the query point is outside of the range of the data points.'
       stop (1)
    end if

    ya(1:n,:) = real(y)
    ya(n+1:2*n,:) = aimag(y)
    
    fd = finite_difference(x,ya)
    jj = 1
    !write(*,*) 'x:', x
    !write(*,*) 'xq:', xquery
    do ii=1,dim
       do while(x(jj) < xquery(ii) + tiny(x(1)))
          jj = jj+1
       end do
       !write(*,*) 'jj:',jj,', size x:',size(x)
       xk1 = x(jj)
       xk = x(jj-1)
       !write(*,*) 'xk',xk,'xk1: ',xk1
       t = (xquery(ii) - xk)/(xk1-xk)
       h00 = 2*t**3-3*t**2+1
       h10 = t**3 - 2*t**2+t
       h01 = -2*t**3+3*t**2
       h11 = t**3-t**2
       !write(*,*) fd(:,jj-1)
       yquery_t = h00*ya(:,jj-1) + h10*fd(:,jj-1)*(xk1-xk) + h01*ya(:,jj) + h11*fd(:,jj)*(xk1-xk)
       yquery(:,ii) = cmplx(yquery_t(1:n),yquery_t(n+1:2*n),fp)
    end do
  end function interp1_complex
  !------------------------------------------
  ! Performs cubic Hermite interpolation for complex numbers.
  function interp1_real(x,y,xquery) result(yquery)
    implicit none
    real(fp),dimension(:),intent(in) :: x
    real(fp), dimension(:,:),intent(in) :: y
    real(fp),dimension(:),intent(in) :: xquery
    real(fp), dimension(size(y,1),size(xquery)) :: yquery
    integer(ip) :: dim_in,dim,ii,jj,n
    real(fp),dimension(2*size(y,1),size(x)) :: ya,fd
    real(fp) :: xk,xk1,t,h00,h01,h10,h11
    real(fp) :: yquery_t(size(y,1))
    dim_in = size(x)
    dim = size(xquery)
    n = size(y,1)
    ! First calculate the finite differences
    if(x(dim_in) < xquery(dim) .or. xquery(1) < x(1)) then
       write(*,*) 'In interp1, the query point is outside of the range of the data points.'
       stop (1)
    end if

    fd = finite_difference(x,y)
    jj = 1
    !write(*,*) 'x:', x
    !write(*,*) 'xq:', xquery
    do ii=1,dim
       do while(x(jj) < xquery(ii) + tiny(x(1)))
          jj = jj+1
       end do
       !write(*,*) 'jj:',jj,', size x:',size(x)
       xk1 = x(jj)
       xk = x(jj-1)
       !write(*,*) 'xk',xk,'xk1: ',xk1
       t = (xquery(ii) - xk)/(xk1-xk)
       h00 = 2*t**3-3*t**2+1
       h10 = t**3 - 2*t**2+t
       h01 = -2*t**3+3*t**2
       h11 = t**3-t**2
       !write(*,*) fd(:,jj-1)
       yquery_t = h00*y(:,jj-1) + h10*fd(:,jj-1)*(xk1-xk) + h01*y(:,jj) + h11*fd(:,jj)*(xk1-xk)
       yquery(:,ii) = yquery_t
    end do
  end function interp1_real
  !------------------------------------------
  function finite_difference(x,y) result(fd)
    implicit none
    real(fp),dimension(:),intent(in) :: x
    real(fp),dimension(:,:),intent(in) :: y
    real(fp),dimension(size(y,1),size(y,2)) :: fd
    integer(ip) :: ii,dim,n
    n = size(y,1)
    dim = size(y,2)
    ! First and last points
    fd(:,1) = (y(:,2) - y(:,1))/(x(2) - x(1))
    fd(:,dim) = (y(:,dim) - y(:,dim-1))/(x(dim) - x(dim-1))

    ! This is calculated in a stupid way. The below formula simplifies to central difference with delta = x(ii+1)-x(ii-1).
    do ii=2,dim-1
       fd(:,ii) = (y(:,ii+1) - y(:,ii))/(2.0_fp*(x(ii+1) - x(ii))) &
            + (y(:,ii) - y(:,ii-1))/(2.0_fp*(x(ii) - x(ii-1)))
    end do
  end function finite_difference
!#############################################
  function trace_real(A) result(tr)
    implicit none
    real(fp), dimension(:,:),intent(in) :: A
    real(fp) :: tr
    integer(ip) :: ii,nn
    nn = size(A,1)
    tr = 0.0_fp
    do ii=1,nn
       tr = tr + A(ii,ii)
    end do
  end function trace_real
  function trace_complex(A) result(tr)
    implicit none
    complex(fp), dimension(:,:),intent(in) :: A
    complex(fp) :: tr
    integer(ip) :: ii,nn

    nn = size(A,1)
    tr = 0.0_fp
    do ii=1,nn
       tr = tr + A(ii,ii)
    end do
  end function trace_complex
!#############################################
  function trace_abs_real(A) result(tr)
    implicit none
    real(fp), dimension(:,:),intent(in) :: A
    real(fp) :: tr
    integer(ip) :: ii,nn
    nn = size(A,1)
    tr = 0.0_fp
    do ii=1,nn
       tr = tr + abs(A(ii,ii))
    end do
  end function trace_abs_real
  function trace_abs_complex(A) result(tr)
    implicit none
    complex(fp), dimension(:,:),intent(in) :: A
    real(fp) :: tr
    integer(ip) :: ii,nn

    nn = size(A,1)
    tr = 0.0_fp
    do ii=1,nn
       tr = tr + sqrt(conjg(A(ii,ii))*A(ii,ii))
    end do
  end function trace_abs_complex
!###########################################################
  function l2norm_real(v) result(norm)
    implicit none
    real(fp), dimension(:),intent(in) :: v
    real(fp) :: norm
    norm = sqrt(dot_product(v,v))
  end function l2norm_real
!-----------------------------------------------------------
    function l2norm_complex(v) result(norm)
    implicit none
    complex(fp), dimension(:),intent(in) :: v
    real(fp) :: norm
    norm = sqrt(dot_product(v,v))
  end function l2norm_complex
!############################################################
  function dephasing_operator(dim,l1,l2) result(sl1l2)
    implicit none
    integer(ip), intent(in) :: dim,l1,l2
    complex(fp), dimension(dim,dim) :: sl1l2
    if(l1>dim .or. l2>dim) then
       write(*,*) 'Improper values for parameters in dephasing_operator.'
       stop (1)
    endif
    
    sl1l2 = (0.0_fp,0.0_fp)
    !sl1l2(l1,l1) = -1.0_fp/sqrt(2.0_fp)
    !sl1l2(l2,l2) = 1.0_fp/sqrt(2.0_fp)
    !write(*,*) 'Warning, modified stuff in dephasing_operator!.'
    sl1l2(l1,l1) = 0.0_fp
    sl1l2(l2,l2) = sqrt(2.0_fp)
  end function dephasing_operator
!############################################################
  function pauli_x(dim,l1,l2) result(op)
    implicit none
    integer(ip),intent(in) :: dim,l1,l2
    complex(fp),dimension(dim,dim) :: op
    if(l1 > dim .or. l2 > dim .or. l1 < 1 .or. l2 < 1) then
       write(*,*) 'Incorrect parameters in pauli_x.'
       stop
    endif
    op = (0.0_fp,0.0_fp)
    op(l1,l2) = (1.0_fp,0.0_fp)
    op(l2,l1) = (1.0_fp,0.0_fp)
  end function pauli_x

  function pauli_y(dim,l1,l2) result(op)
    implicit none
    integer(ip),intent(in) :: dim,l1,l2
    complex(fp),dimension(dim,dim) :: op
    if(l1 > dim .or. l2 > dim .or. l1 < 1 .or. l2 < 1) then
       write(*,*) 'Incorrect parameters in pauli_y.'
       stop
    endif
    op = (0.0_fp,0.0_fp)
    op(l1,l2) = (0.0_fp,-1.0_fp)
    op(l2,l1) = (0.0_fp,1.0_fp)
  end function pauli_y
!#############################################
  ! Generates the generalized Gell-Mann matrices for a system with dim dimensions. ggm is the output matrix with dimensions (dim,dim,dim**2-1). Uses convention, where sigma_z = [-1,0;0,1]. Equations form "Bloch vectors for qudits", Reinhold A. Bertlmann and Philipp Krammer.
  subroutine generate_ggm(dim,ggm,ggmnorm)
    implicit none
    integer(ip),intent(in) :: dim
    complex(fp),intent(out),dimension(dim,dim,dim**2) :: ggm
    real(fp),intent(out),optional,dimension(dim**2) :: ggmnorm
    integer(ip) :: jj,kk,ll,nn
    real(fp) :: ll_factor

    if(dim<2) then
       write(*,*) 'Dimension < 2 in generate_ggm.'
       stop
    end if

    ggm = (0.0_fp,0.0_fp)
    nn = 1
    
    ! Identity
    do jj=1,dim
       ggm(jj,jj,nn) = 1.0_fp
    end do
    nn = nn+1

    do kk=2,dim
       do jj=1,kk-1
          ggm(jj,kk,nn) = 1.0_fp
          ggm(kk,jj,nn) = 1.0_fp
          nn = nn+1
       end do
    end do
    do kk=2,dim
       do jj=1,kk-1
          ggm(jj,kk,nn) = -(0.0_fp,1.0_fp)
          ggm(kk,jj,nn) = (0.0_fp,1.0_fp)
          nn = nn+1
       end do
    end do
    do ll=1,dim-1
       ll_factor = sqrt(2.0_fp/(ll*(ll+1)))
       do jj=1,ll
          ggm(jj,jj,nn) = -ll_factor
       end do
       ggm(ll+1,ll+1,nn) = ll*ll_factor
       nn = nn+1
    end do

    if(present(ggmnorm)) then
       do kk=1,dim**2
          ggmnorm(kk) = trace(matmul(ggm(:,:,kk),ggm(:,:,kk)))
       end do
    endif
  end subroutine generate_ggm
!#############################################
  subroutine density_matrix_to_ggm(rho,ggm,b)
    implicit none
    complex(fp),intent(in),dimension(:,:) :: rho
    complex(fp),intent(in),dimension(size(rho,1),size(rho,2),size(rho,1)**2) :: ggm
    real(fp),intent(out),dimension(size(rho,1)**2) :: b
    integer(ip) :: ii,nn
    integer(ip) :: c1,c2
    complex(fp),dimension(size(rho,1),size(rho,2)) :: ggmrho
    nn = size(rho,1)
    !call system_clock(c1)
    do ii=1,(nn**2)
       !b(ii) = real(trace(matmul(ggm(:,:,ii),rho)))
       call zgemm('n','n',nn,nn,nn,(1.0_fp,0.0_fp),ggm(:,:,ii),nn,rho,nn,&
            (0.0_fp,0.0_fp),ggmrho,nn)
       b(ii) = real(trace(ggmrho))
    end do
    !call system_clock(c2)
    !write(*,*) 'Density_matrix_to_ggm. Duration: ', c2-c1,'ms.'
  end subroutine density_matrix_to_ggm
!----------------------------------------------!
  subroutine ggm_to_density_matrix(rho,ggm,ggmnorm,b,dim)
    implicit none
    integer(ip),intent(in) :: dim
    real(fp),intent(in),dimension(dim**2) :: b
    complex(fp),intent(out),dimension(dim,dim) :: rho
    complex(fp),intent(out),dimension(dim,dim,dim**2) :: ggm
    real(fp),intent(in),dimension(dim**2) :: ggmnorm
    integer(ip) :: ii
    integer(ip) :: c1,c2
    !rho = get_eye(dim)/dim
    rho = 0.0_fp
    !write(*,*) 'After eye:'
    !call print_matrix(rho)
    ! rho = I/d + \sum_i=1^{d^2-1} b_i A_i/N_i, where N_i = Tr A^2
    !call system_clock(c1)
    do ii=1,dim**2
       !write(*,*) 'b(',ii,')=', b(ii)
       !rho = rho + b(ii)*ggm(:,:,ii)/trace(matmul(ggm(:,:,ii),ggm(:,:,ii))) ! included the normalization
       rho = rho + b(ii)/ggmnorm(ii)*ggm(:,:,ii)! included the normalization

       !call print_matrix(rho)
    end do
    !call system_clock(c2)
    !write(*,*) 'ggm_to_density_matrix. Duration: ', c2-c1,'ms.'
    !rho = rho/dim
  end subroutine ggm_to_density_matrix
!#############################################
  function trapz(x,y) result(r)
    implicit none
    real(fp), dimension(:),intent(in) :: x,y
    real(fp) :: r
    integer(ip) :: n,ii
    n = size(y)
    r = 0.0_fp
    do ii=2,n
       r = r + (x(ii) - x(ii-1))*(y(ii) + y(ii-1))
    end do
    r = r/2.0_fp
  end function trapz

end module linfuncs
