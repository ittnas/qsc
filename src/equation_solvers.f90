module equation_solvers
  use utils, only: fp,ip
  !use hamiltonians ! To fix compilation bug with fclass
  use fclass

  abstract interface
     subroutine nonlin_solver(f,x0,x,tol_in,maxiter_in,status,dx_in)
       use fclass
       use utils
       implicit none
       class(fd),intent(in) :: f
       complex(fp), dimension(:),intent(in) :: x0
       complex(fp), dimension(:),intent(out) :: x
       real(fp),optional,intent(in) :: tol_in
       integer(ip),optional,intent(in) :: maxiter_in
       integer(ip),optional,intent(out) :: status
       real(fp),optional,intent(in) :: dx_in
     end subroutine nonlin_solver
  end interface
contains
  ! Solves the equation f(x) - x = 0, i.e. it does not have the standard form.
  subroutine picard_iteration(f,x0,x,tol_in,maxiter_in,status,dx_in)
    implicit none
    class(fd),intent(in) :: f
    complex(fp), dimension(:),intent(in) :: x0
    complex(fp), dimension(size(x0)),intent(out) :: x
    real(fp),optional,intent(in) :: tol_in
    integer(ip),optional,intent(in) :: maxiter_in
    real(fp) :: delta,tol
    complex(fp),dimension(size(x0)) :: xn
    integer(ip) :: maxiter,niter
    integer(ip),optional,intent(out) :: status
    integer(ip) :: status_temp
    real(fp),optional,intent(in) :: dx_in
    
    if(present(tol_in)) then
       tol = tol_in
    else
       tol = 1.0_fp
    end if

    if(present(maxiter_in)) then
       maxiter = maxiter_in
    else
       maxiter = 5000
    end if

    delta = huge(1.0_fp)
    xn = x0
    x = xn
    niter = 0

    do
       if(delta < tol) then
          !write(*,*) 'niter', niter,', delta',delta
          status_temp = 0
          exit ! We are done
       end if
       if(niter > maxiter) then
          !write(*,*) 'Warning, maximum number of iterations reached in picard solver. Required tolerance not reached.'
          status_temp = 1
          exit
       endif
       call f%value_at(0.0_fp,xn,x)
       delta = real(sqrt(sum((x - xn)*conjg(x-xn))))
       xn = x
       niter = niter+1
    end do
    if(present(status)) then
       status = status_temp
    end if
  end subroutine picard_iteration
!------------------------------------------------------------!
  subroutine quasi_newton(f,x0,x,tol_in,maxiter_in,status,dx_in)
    use linfuncs
    implicit none
    class(fd),intent(in) :: f
    complex(fp), dimension(:),intent(in) :: x0
    complex(fp), dimension(size(x0)),intent(out) :: x
    real(fp),optional,intent(in) :: tol_in
    integer(ip),optional,intent(in) :: maxiter_in
    real(fp),optional,intent(in) :: dx_in
    integer(ip),optional,intent(out) :: status
    real(fp) :: delta,tol,dx
    complex(fp),dimension(size(x0)) :: xn
    integer(ip) :: maxiter,niter
    integer(ip) :: status_temp,status_dgesv
    real(fp),dimension(2*size(x0),2*size(x0)) :: J
    real(fp),dimension(2*size(x0)) :: b
    real(fp),dimension(2*size(x0)) :: bt
    complex(fp),dimension(size(x0)) :: fx
    integer(ip),dimension(2*size(x0)) :: ipiv

    if(present(tol_in)) then
       tol = tol_in
    else
       tol = 1.0_fp
    end if

    if(present(maxiter_in)) then
       maxiter = maxiter_in
    else
       maxiter = 100
    end if

    if(present(dx_in)) then
       dx = dx_in
    else
       dx = 0.1_fp
    endif

    xn = x0
    x = x0

    delta = huge(1.0_fp)
    niter = 0
    do
       if(delta < tol) then
          !write(*,*) 'niter', niter,', delta',delta
          status_temp = 0
          exit ! We are done
       end if

       if(niter > maxiter) then
          write(*,*) 'Warning, maximum number of iterations reached in quasi newton. Required tolerance not reached.'
          !write(*,*) 'niter', niter,', delta',delta
          status_temp = 1
          exit
       endif
       call f%value_at(0.0_fp,xn,fx)
       b(1:size(x0)) = -real(fx)
       b((size(x0)+1):(2*size(x0))) = -aimag(fx)
       bt = b
       J = jacobian_approximation_complex(f,x,dx)
       !SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       call dgesv(size(2*x0),1,J,size(2*x0),ipiv,b,size(2*x0),status_dgesv)
       !write(*,*) b
       if(status_dgesv .ne. 0) then
          status_temp = 2
          write(*,*) 'dgesv could not find solution to the system of equations in quasi-newton. Status:', status_dgesv
          !call print_matrix(J)
          !write(*,*) b
          exit
       endif
       x = xn + b(1:size(x0)) + (0.0_fp,1.0_fp)*b((1+size(x0)):(2*size(x0)))
       !delta = real(sqrt(sum((x - xn)*conjg(x-xn))))
       delta = sqrt(sum((matmul(J,b) - bt)**2))
       xn = x
       niter = niter+1
       !call print_matrix(J)
       write(*,*) delta
       write(*,*) b
       !write(*,*) x
    end do
    if(present(status)) then
       status = status_temp
    end if
end subroutine quasi_newton
!------------------------------------------------------------!
subroutine quasi_newton_real(f,x0,x,tol_in,maxiter_in,status,dx_in)
  use linfuncs
  implicit none
  class(fd),intent(in) :: f
  complex(fp), dimension(:),intent(in) :: x0
  complex(fp), dimension(size(x0)),intent(out) :: x
  real(fp),optional,intent(in) :: tol_in
  integer(ip),optional,intent(in) :: maxiter_in
  real(fp),optional,intent(in) :: dx_in
  integer(ip),optional,intent(out) :: status
  real(fp) :: delta,tol,dx
  complex(fp),dimension(size(x0)) :: xn
  integer(ip) :: maxiter,niter
  integer(ip) :: status_temp,status_dgesv
  real(fp),dimension(size(x0),size(x0)) :: J
  real(fp),dimension(size(x0)) :: b
  real(fp),dimension(size(x0)) :: bt
  complex(fp),dimension(size(x0)) :: fx
  integer(ip),dimension(size(x0)) :: ipiv

  if(present(tol_in)) then
       tol = tol_in
    else
       tol = 1.0_fp
    end if

    if(present(maxiter_in)) then
       maxiter = maxiter_in
    else
       maxiter = 100
    end if

    if(present(dx_in)) then
       dx = dx_in
    else
       dx = 0.1_fp
    endif

    xn = x0
    x = x0

    delta = huge(1.0_fp)
    niter = 0
    do
       if(delta < tol) then
          !write(*,*) 'niter', niter,', delta',delta
          status_temp = 0
          exit ! We are done
       end if

       if(niter > maxiter) then
          write(*,*) 'Warning, maximum number of iterations reached in quasi-newton. Required tolerance not reached.'
          !write(*,*) 'niter', niter,', delta',delta
          status_temp = 1
          exit
       endif
       call f%value_at(0.0_fp,xn,fx)
       b = -real(fx)
       !write(*,*) 'b:', b
       !b((size(x0)+1):(2*size(x0))) = -aimag(fx)
       bt = b ! need a copy because dgesv destroys the original b
       call jacobian_approximation_real(f,x,dx,J)
       !call print_matrix(J)
       !SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       call dgesv(size(x0),1,J,size(x0),ipiv,b,size(x0),status_dgesv)
       !write(*,*) 'dgesv result:', b
       if(status_dgesv .ne. 0) then
          status_temp = 2
          write(*,*) 'dgesv could not find solution to the system of equations in quasi-newton. Status:', status_dgesv
          !call print_matrix(J)
          !write(*,*) b
          exit
       endif
       x = xn + b
       !delta = real(sqrt(sum((x - xn)*conjg(x-xn))))
       delta =sum(real(b)**2)
       !delta = sqrt(sum((matmul(J,b) - bt)**2))
       
       xn = x
       niter = niter+1
       !write(*,*) 'niter',niter
       !call print_matrix(J)
       !write(*,*) delta
       !write(*,*) b
       !write(*,*) x
    end do
    if(present(status)) then
       status = status_temp
    end if
  end subroutine quasi_newton_real
!------------------------------------------------------------!
subroutine spectral_real(f,x0,x,tol_in,maxiter_in,status,dx_in)
  use linfuncs
  implicit none
  class(fd),intent(in) :: f
  complex(fp), dimension(:),intent(in) :: x0
  complex(fp), dimension(:),intent(out) :: x
  real(fp),optional,intent(in) :: tol_in
  integer(ip),optional,intent(in) :: maxiter_in
  real(fp),optional,intent(in) :: dx_in
  integer(ip),optional,intent(out) :: status
  real(fp) :: delta,tol,dx
  integer(ip) :: maxiter,niter,status_temp
  complex(fp),dimension(size(x0)) :: xn,sn
  complex(fp),dimension(size(x0)) :: fn,yn,fn1
  real(fp) :: alpha
  integer(ip) :: c1,c2

  if(present(tol_in)) then
     tol = tol_in
  else
     tol = 1.0_fp
  end if

  if(present(maxiter_in)) then
     maxiter = maxiter_in
  else
     maxiter = 100
  end if

  if(present(dx_in)) then
     dx = dx_in
  else
     dx = 0.1_fp
  endif

  xn = x0
  x = x0
  
  !write(*,*) 'x0', x0
  delta = huge(1.0_fp)
  niter = 0
  
  
  do
     !call system_clock(c1)
     if(niter .eq. 0) then
        call f%value_at(0.0_fp,xn,fn1)
        alpha = min(1.0_fp,1.0_fp/sqrt(real(sum(conjg(fn1)*fn1))))
        if(isnan(alpha)) then
           delta = 0.0_fp ! sum(fn* *fn) = 0
        endif
        
        !write(*,*) 'x0',x0
        !write(*,*) 'fn1',fn1
     end if
     
      if(delta < tol) then
         status_temp = 0
         write(*,*) 'niter', niter,', delta',delta
         exit ! we are done
      end if
      
      if(niter > maxiter) then
         write(*,*) 'Warning, maximum number of iterations reached in spectral solver. Required tolerance not reached.'
         status_temp = 1
         exit
      endif
      
      x = xn - alpha*fn1
      delta = sqrt(real(sum(conjg(fn1)*fn1)))

      sn = x - xn

      xn = x
      fn = fn1

      call f%value_at(0.0_fp,x,fn1)
      yn = fn1 - fn
      alpha = real(sum(conjg(sn)*sn))/real(sum(conjg(sn)*yn))
      !alpha = real(sum(conjg(sn)*yn))/real(sum(conjg(yn)*yn))
      niter = niter + 1
      !write(*,*) 'niter', niter,', delta',delta,', alpha',alpha
      ! write(*,*) 'sn^T*sn',real(sum(conjg(sn)*sn)), ', sn^T*yn', real(sum(conjg(sn)*yn))
      ! write(*,*) 'x',x
      ! write(*,*) 'fn1',fn1
      !call system_clock(c2)
      !write(*,*) 'Duration of one iteration in spectral_real',c2-c1,'ms.'
   end do
    if(present(status)) then
       status = status_temp
    end if
end subroutine spectral_real
!------------------------------------------------------------!
! Finds the roots of equation f(x) = 0 by using the quasi-newton method, by approximating the Jacobian using the Broyden formula. The Jacobian is not fully calculated, but it is updated during every iteration. Consumes more memory than the spectral method.
subroutine broyden_real(f,x0,x,tol_in,maxiter_in,status,dx_in)
  use linfuncs
  implicit none
  class(fd),intent(in) :: f
  complex(fp), dimension(:),intent(in) :: x0
  complex(fp), dimension(size(x0)),intent(out) :: x
  real(fp),optional,intent(in) :: tol_in
  integer(ip),optional,intent(in) :: maxiter_in
  real(fp),optional,intent(in) :: dx_in
  integer(ip),optional,intent(out) :: status
  real(fp) :: delta,tol,dx,delta_x
  integer(ip) :: maxiter,niter,status_temp,status_dgesv
  complex(fp),dimension(size(x0)) :: xn
  complex(fp),dimension(size(x0)) :: fn,fn1
  real(fp),dimension(size(x0),size(x0)) :: J,Jt
  real(fp),dimension(size(x0)) :: bn
  integer(ip),dimension(size(x0)) :: ipiv
  integer(ip) :: c1,c2

  if(present(tol_in)) then
     tol = tol_in
  else
     tol = 1.0_fp
  end if

  if(present(maxiter_in)) then
     maxiter = maxiter_in
  else
     maxiter = 100
  end if

  if(present(dx_in)) then
     dx = dx_in
  else
     dx = 0.1_fp
  endif

  xn = x0
  x = x0

  !write(*,*) 'x0', x0
  delta = huge(1.0_fp)
  niter = 0

  !J = get_eye(size(x0))
  call jacobian_diagonal_real(f,xn,dx,J)
  call f%value_at(0.0_fp,xn,fn)
  !call print_matrix(J)
  !bn = -real(fn)
  delta = sqrt(sum(bn*bn))

  do
     if(delta < tol) then
        write(*,*) 'niter', niter,', delta',delta
        status_temp = 0
        exit ! We are done
     end if
     if(niter > maxiter) then
        write(*,*) 'Warning, maximum number of iterations reached in quasi-newton. Required tolerance not reached.'
        !write(*,*) 'niter', niter,', delta',delta
        status_temp = 1
        exit
     endif
     !write(*,*) 'niter', niter,', delta',delta
     Jt = J
     bn = -real(fn)
     call dgesv(size(x0),1,Jt,size(x0),ipiv,bn,size(x0),status_dgesv)
     if(status_dgesv .ne. 0) then
        status_temp = 2
        write(*,*) 'dgesv could not find solution to the system of equations in quasi-newton. Status:', status_dgesv
        !call print_matrix(J)
        !write(*,*) bn
        exit
     endif

     x = xn + bn
     call f%value_at(0.0_fp,x,fn1)
     !write(*,*) 'J before'
     !call print_matrix(J)
     !write(*,*) 'addition'
     !call print_matrix(matmul(reshape(fn1 - fn - matmul(J,bn),[size(x0),1]),reshape(bn,[1,size(x0)]))/sum(bn**2))
     !write(*,*) 'sum bn',sum(bn**2)
     J = J + matmul(reshape(fn1 - fn - matmul(J,bn),[size(x0),1]),reshape(bn,[1,size(x0)]))/sum(bn**2)
     !write(*,*) 'J after'
     !call print_matrix(J)
     !write(*,*) 'bn',bn
     !write(*,*) 'fn1',fn1
     !write(*,*) 'fn',fn
     !write(*,*) 'fn1-fn - J*bn',reshape(fn1 - fn - matmul(J,bn),[size(x0),1])
     !write(*,*) 'bn^T',reshape(bn,[1,size(x0)])
     delta = sqrt(real(sum(conjg(fn1)*fn1)))
     !delta_x = sqrt(sum(bn*bn))
     !write(*,*) bn
     !write(*,*) x
     fn = fn1
     xn = x
     niter = niter+1
     !call print_matrix(J)
     !write(*,*) 'x:',x
     !write(*,*) 'bn',bn
     !write(*,*) 'J:'
     !call print_matrix(J)
  end do
  if(present(status)) then
     status = status_temp
  end if
end subroutine broyden_real

  !#############################################
  ! This function calculates the Jacobian of the given function at x using dx as the step in difference. Works only for time independent functions f.
function jacobian_approximation_complex(f,x,dx) result(J)
  use linfuncs
  use fclass
  implicit none
  class(fd),intent(in) :: f
  complex(fp),dimension(:),intent(in) :: x
  real(fp),intent(in) :: dx
  real(fp),dimension(size(x)*2,size(x)*2) :: J
  integer(ip) :: ii,jj,nn,kk
  complex(fp),dimension(size(x)) :: fxdx,fx
  complex(fp),dimension(size(x)) :: zv
  J = 0.0_fp
  nn = size(x)
  call f%value_at(0.0_fp,x,fx) ! Makes the calculation at t=0
  !write(*,*) 'fx:',fx
  do ii=1,nn
     do jj=1,2
        zv = (0.0_fp,0.0_fp) ! This is not very efficient
        zv(ii) = (0.0_fp,1.0_fp)**(jj-1)*dx
        write(*,*) 'zv(ii):',zv(ii)
        call f%value_at(0.0_fp,x+zv,fxdx)
        
        !write(*,*) 'fxdx:',fxdx
        ! Calculates the difference, and takes either the real or the imaginary part of it.
        fxdx = (fxdx - fx)/dx
        !J(((nn*(jj-1)+1)):(nn*jj),ii) = real((0.0_fp,1.0_fp)**(jj-1)*(fxdx - fx)/dx)
        J(((nn*(jj-1)+1)):(nn*jj),ii) = real(fxdx)
        J(((nn*(jj-1)+1)):(nn*jj),nn+ii) = aimag(fxdx)
        
        if(jj .eq. 1) then
           write(*,*) 'Changing the ',ii,':th real element. The effect is:'
             do kk=1,nn
                write(*,*) real(fxdx(kk)) + (0.0_fp,1.0_fp)*aimag(fxdx(kk))
             end do
             
          else
             write(*,*) 'Changing the ',ii,':th imag element. The effect is:'           
             do kk=1,nn
                write(*,*) real(fxdx(kk)) + (0.0_fp,1.0_fp)*aimag(fxdx(kk))
             end do
          endif
          !call print_matrix(J)
          !write(*,*) (0.0_fp,1.0_fp)**(jj-1)*(fxdx - fx)/dx
          !write(*,*) ((nn*(jj-1)+1)), ':',(nn*jj)
          !call print_matrix(J)
       end do
    end do
  end function jacobian_approximation_complex
  !#############################################
  subroutine jacobian_approximation_real(f,x,dx,J)
    use linfuncs
    use fclass
    implicit none
    class(fd),intent(in) :: f
    complex(fp),dimension(:),intent(in) :: x
    real(fp),intent(in) :: dx
    real(fp),dimension(size(x),size(x)),intent(out) :: J
    integer(ip) :: ii,jj,nn,kk
    complex(fp),dimension(size(x)) :: fxdx,fx
    complex(fp),dimension(size(x)) :: zv
    J = 0.0_fp
    nn = size(x)
    !write(*,*) 'Entering function.'
    call f%value_at(0.0_fp,x,fx) ! Makes the calculation at t=0
    !write(*,*) 'x',x
    !write(*,*) 'fx:',fx
    do ii=1,nn
       zv = (0.0_fp,0.0_fp) ! This is not very efficient
       zv(ii) = dx
       !write(*,*) 'zv(ii):',zv(ii)
       !write(*,*) 'x+zv',x+zv
       call f%value_at(0.0_fp,x+zv,fxdx)
       
       !write(*,*) 'fxdx:',fxdx
       ! Calculates the difference, and takes either the real or the imaginary part of it.
       fxdx = (fxdx - fx)/dx
       !write(*,*) 'fxdx',real(fxdx)
       !J(((nn*(jj-1)+1)):(nn*jj),ii) = real((0.0_fp,1.0_fp)**(jj-1)*(fxdx - fx)/dx)
       J(:,ii) = real(fxdx)
       !write(*,*) 'J(:,ii)',J(:,ii)
       ! write(*,*) 'Changing the ',ii,':th real element. The effect is:'
       ! do kk=1,nn
       !    write(*,*) real(fxdx(kk))
       ! end do

       !call print_matrix(J)
       !write(*,*) (0.0_fp,1.0_fp)**(jj-1)*(fxdx - fx)/dx
       !write(*,*) ((nn*(jj-1)+1)), ':',(nn*jj)
    end do
  end subroutine jacobian_approximation_real

  !#############################################
  subroutine jacobian_diagonal_real(f,x,dx,J)
    use linfuncs
    use fclass
    implicit none
    class(fd),intent(in) :: f
    complex(fp),dimension(:),intent(in) :: x
    real(fp),intent(in) :: dx
    real(fp),dimension(size(x),size(x)),intent(out) :: J
    integer(ip) :: nn,ii
    complex(fp),dimension(size(x)) :: fx,fxdx
    complex(fp),dimension(size(x)) :: zv

    !zv = (0.0_fp,0.0_fp)
    nn = size(x)
    J = 0.0_fp
    call f%value_at(0.0_fp,x,fx) ! Makes the calculation at t=0
    zv = dx
    call f%value_at(0.0_fp,x+zv,fxdx)

    !write(*,*) 'fxdx:',fxdx
    ! Calculates the difference, and takes either the real or the imaginary part of it.
    fxdx = (fxdx - fx)/dx
    !write(*,*) 'fxdx',real(fxdx)
    do ii=1,nn
       J(ii,ii) = real(fxdx(ii))
    end do
  end subroutine jacobian_diagonal_real

end module equation_solvers
