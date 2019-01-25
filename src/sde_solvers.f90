module sde_solvers
  use utils,only: fp,ip
  
  abstract interface
     subroutine sde_solver(x,y0,f,y,j,rgen,dx_in,step_split_factor_in,tol_in)
       use utils,only:fp,ip
       use hamiltonians
       use sclass
       use random
       implicit none
       real(fp), dimension(:),intent(in) :: x
       complex(fp), dimension(:),intent(in) :: y0
       complex(fp), dimension(:,:),intent(out) :: y
       real(fp),dimension(:,:),intent(out) :: j
       real(fp),intent(in),optional :: dx_in,tol_in
       integer(ip),intent(in),optional :: step_split_factor_in
       class(sd),intent(in) :: f
       class(rng),intent(inout) :: rgen
     end subroutine sde_solver
  end interface
contains
  subroutine test_sde_solver(x,y0,f,y,j,rgen,dx_in,step_split_factor_in,tol_in)
    use utils,only:fp,ip
    use hamiltonians
    use sclass
    use random
    implicit none
    real(fp), dimension(:),intent(in) :: x
    complex(fp), dimension(:),intent(in) :: y0
    complex(fp), dimension(:,:),intent(out) :: y
    real(fp),dimension(:,:),intent(out) :: j
    real(fp),intent(in),optional :: dx_in,tol_in
    integer(ip),intent(in),optional :: step_split_factor_in
    class(sd),intent(in) :: f
    class(rng),intent(inout) :: rgen
    !y = reshape(rgen%norm(size(y)),[size(y,1),size(y,2)])
    !call rgen%norm(y,size(y))
    j = 0.0_fp
  end subroutine test_sde_solver

  subroutine euler_maryama(x,y0,f,y,j,rgen,dx_in,step_split_factor_in,tol_in)
    use utils,only:fp,ip
    use sclass
    use random
    use linfuncs
    implicit none
    real(fp), dimension(:),intent(in) :: x
    complex(fp), dimension(:),intent(in) :: y0
    complex(fp), dimension(:,:),intent(out) :: y
    real(fp),dimension(:,:),intent(out) :: j
    real(fp),intent(in),optional :: dx_in,tol_in
    integer(ip),intent(in),optional :: step_split_factor_in
    class(sd),intent(in) :: f
    class(rng),intent(inout) :: rgen
    integer(ip) :: n,dim,ii,jj,step_split_factor
    real(fp) :: dx,xc,xn,dxc
    complex(fp),dimension(size(y0)) :: fn,yn
    complex(fp),dimension(size(y0),f%nobs) :: s
    logical :: write_next 
    real(fp),dimension(f%nobs) :: dw
    real(fp) :: dn
    n = size(x)
    dim = size(y0)

    if(n < 2) then
       write(*,*) 'At least 2 query points must be defined.'
       stop
    end if
    
    if(size(y,2) .NE. n) then
       write(*,*) 'Size of output array is not equal to query points.'
       stop
    end if

    if(present(dx_in)) then
       dx = dx_in
       if(dx > (x(2) - x(1))/2.0_fp) then
          dx = x(2) - x(1)
       end if
    else
       dx=x(2) - x(1)
    endif

    if(present(step_split_factor_in)) then
       step_split_factor = step_split_factor_in
    else
       step_split_factor = 0
    end if

    if(step_split_factor > 0) then
       dx = dx/real(2**step_split_factor,fp)
       write(*,*) 'Warning: Step splitting is not properly implemented in E-M algorithm.'
    end if

    y(:,1) = y0
    yn = y0
    ii = 2
    write_next = .false.
    dn = 0.0_fp
    do while(ii <= n)
       write_next = .false.
       xc = x(ii-1)
       xn = x(ii)
       do while(.true.)
          if(xc + dx + dx/1000_fp >= xn) then
             dxc = xn - xc 
             write_next = .true.
          else
             dxc = dx
          end if
          !write(*,*) dx,dxc,xc,xn
          call f%value_at(xc,yn,fn,s)
          call rgen%norm(dw,f%nobs,0.0_fp,1.0_fp)
          dw = sqrt(dxc)*dw
          !dw = sqrt(dxc)*rgen%norm(f%nobs,0.0_fp,1.0_fp)
          yn = yn + dxc*fn
          do jj=1,f%nobs
             yn = yn + s(:,jj)*dw(jj)
          end do
          ! Normalization. Makes the solver more stable, but only works for schrödinger equation type
          !yn = yn/l2norm(yn)
          call f%normalize(yn)
          !write(*,*) l2norm(yn)
          !jj = jj+1
          xc = xc + dxc
          if(write_next) then
             !write(*,*) jj
             y(:,ii) = yn
             exit
          end if
       end do
       ii = ii+1
    end do
    !write(*,*) ii-1,jj

  end subroutine euler_maryama

  ! This is an explicit strong order 1.0 stochastic Runge-Kutta. Can be found in Kloeden&Platen pp. XXIX
  subroutine platen(x,y0,f,y,j,rgen,dx_in,step_split_factor_in,tol_in)
    use utils,only:fp,ip
    use sclass
    use random
    use linfuncs
    implicit none
    real(fp), dimension(:),intent(in) :: x
    complex(fp), dimension(:),intent(in) :: y0
    complex(fp), dimension(:,:),intent(out) :: y
    real(fp),dimension(:,:),intent(out) :: j
    real(fp),intent(in),optional :: dx_in,tol_in
    class(sd),intent(in) :: f
    class(rng),intent(inout) :: rgen
    integer(ip),intent(in),optional :: step_split_factor_in
    integer(ip) :: n,dim,ii,jj,kk,nrand,mm,step_split_factor
    real(fp) :: dx,xc,xn,dxc
    real(fp),allocatable,dimension(:,:) :: rnd
    
    complex(fp),dimension(size(y0),f%nobs) :: sn,sp,sm,stilde
    complex(fp),dimension(size(y0)) :: dn,dtilde,dp,dm
    real(fp),dimension(f%nobs) :: dw,j_avg,dw_sum,j_avg_sum
    logical :: write_next
    complex(fp),dimension(size(y0)) :: yn,ym,yp,ytilde
    integer(ip) :: ll,sint
    n = size(x)

    dim = size(y0)
    if(n < 2) then
       write(*,*) 'At least 2 query points must be defined.'
       stop
    end if
    
    if(size(y,2) .NE. n) then
       write(*,*) 'Size of output array is not equal to query points.'
       stop
    end if
    if(present(dx_in)) then
       dx = dx_in
       if(dx > (x(2) - x(1))/2.0_fp) then
          dx = x(2) - x(1)
       end if
    else
       dx=x(2) - x(1)
    endif

    if(present(step_split_factor_in)) then
       step_split_factor = step_split_factor_in
    else
       step_split_factor = 0
    end if

    if(step_split_factor > 0) then
       dx = dx/real(2**step_split_factor,fp)
    end if
    j = 0.0_fp
    y(:,1) = y0
    yn = y0
    ii = 2
    kk = 1
    ll = 1
    !niter = 0

    nrand = (step_split_factor+1)*(ceiling((x(size(x)) - x(1))/dx)+size(x))
    allocate(rnd(f%nobs,nrand))
    !rnd = reshape(rgen%norm(nrand*f%nobs,0.0_fp,1.0_fp),[f%nobs,nrand])
    call rgen%norm(rnd,nrand*f%nobs,0.0_fp,1.0_fp)
    j_avg_sum = 0.0_fp
    dw_sum = 0.0_fp
    write_next = .false.
    do while(ii <= n)
       write_next = .false.
       xc = x(ii-1)
       xn = x(ii)
       do while(.true.)
          if(xc + dx + dx/1000_fp >= xn) then
             dxc = xn - xc 
             write_next = .true.
          else
             dxc = dx
          end if
          call f%value_at(xc,yn,dn,sn)

          ytilde = yn + dn*dxc
          yp = ytilde
          ym = yp
          if(step_split_factor>0) then
             ! Idea for splitting can be found in Matti Silveri's notes.
             ! Split once:
             ! dw_tilde_2n = 1/2*dw_n + 1/2*dz_n
             ! dw_n*dw_n = dt, dQ_n*dQ_n = dt
             ! dw_tilde_n*dw_tilde_n = 1/2*dt
             ! Split twice:
             ! dw_hat_4n = 1/4*dw_n+1/4*dz_n+1/2*dQ_n
             ! dQ_n*dQ_n = 1/2 dt
             
             !dw = 0.5_fp*sqrt(dxc)*(rnd(:,kk) + rnd(:,kk))
             dw = sqrt(2**(step_split_factor)*dxc)*rnd(:,kk)*0.5_fp**step_split_factor
             !dw = 0.0_fp
             ! The contributions are added from last to first (...,dQ,dZ,dW)
             do mm=0,step_split_factor-1
                sint = 1
                !dw = dw + sqrt(2**(step_split_factor - mm + 1)*dxc)*&
                !rnd(:,kk*(step_split_factor - mm+1))*0.5_fp**(step_split_factor-mm+1)*&
                !(-1)**(rshift(and(ll-1,lshift(sint,mm)),mm))
                dw = dw + sqrt(2**(mm + 1)*dxc)*&
                     rnd(:,kk*(step_split_factor - mm+1))*0.5_fp**(mm+1)*&
                     (-1)**(rshift(and(ll-1,lshift(sint,mm)),mm))
                !write(*,*) 'll:',ll,', mm:',mm,', ',(-1)**(rshift(and(ll-1,lshift(sint,mm)),mm))
             end do
          else
             dw = sqrt(dxc)*rnd(:,kk)
             !dw = rnd(:,kk)/sqrt(dxc)
          endif
          do jj=1,f%nobs
             ytilde = ytilde + dw(jj)*sn(:,jj)
             yp = yp + sn(:,jj)*sqrt(dxc)
             ym = yp - sn(:,jj)*sqrt(dxc)
          end do
          call f%value_at(xc+dxc,ytilde,dtilde,stilde,j_avg)
          call f%value_at(xc+dxc,yp,dp,sp)
          call f%value_at(xc+dxc,ym,dm,sm)
          yn = yn + 0.5_fp*(dtilde+dn)*dxc
          do jj=1,f%nobs
             yn = yn + 0.25_fp*(sp(:,jj)+sm(:,jj)+2.0_fp*sn(:,jj))*dw(jj) &
                  +0.25_fp*(sp(:,jj)-sm(:,jj))*(dw(jj)**2 - dxc)/sqrt(dxc)
          end do
          ! Normalization. Makes the solver more stable, but only works for schrödinger equation type
          !yn = yn/l2norm(yn). Should work now, because the normalization is in f.
          call f%normalize(yn)
          j_avg_sum = j_avg_sum + j_avg*dxc
          dw_sum = dw_sum + dw
          !write(*,*) ii,dw_sum,dw,j_avg_sum
          !niter = niter + 1
          xc = xc + dxc
          if(step_split_factor>0) then
             if(ll == 2**step_split_factor) then
                kk = kk+1
                ll = 1
             else
                ll=ll+1
             end if
          else
             kk = kk+1
          end if
          if(write_next) then
             !write(*,*) jj
             y(:,ii) = yn
             j(:,ii) = j_avg_sum/(x(ii)-x(ii-1)) + dw_sum/(x(ii)-x(ii-1))
             !write(*,*) ii,x(ii) - x(ii-1)
             j_avg_sum = 0.0_fp
             dw_sum = 0.0_fp
             exit
          end if
       end do
       ii = ii+1
    end do
    !write(*,*) 'niter: ',niter
  end subroutine platen

  ! Strong order 1.5 implicit Runge-Kutta (Kloeden & Platen pp. 401)
  subroutine isrk15(x,y0,f,y,j,rgen,dx_in,step_split_factor_in,tol_in)
    use utils,only:fp,ip
    use sclass
    use random
    use linfuncs
    use fclass
    use equation_solvers
    implicit none
    
    real(fp), dimension(:),intent(in) :: x
    complex(fp), dimension(:),intent(in) :: y0
    complex(fp), dimension(:,:),intent(out) :: y
    real(fp),dimension(:,:),intent(out) :: j
    real(fp),intent(in),optional :: dx_in,tol_in
    class(sd),intent(in) :: f
    class(rng),intent(inout) :: rgen
    integer(ip),intent(in),optional :: step_split_factor_in

    integer(ip) :: n,dim,ii,jj,kk,ll,nrand,mm,step_split_factor
    real(fp) :: dx,xc,xn,dxc,sqrtdxc,tol
    complex(fp),dimension(size(y0)) :: yn,gp,gm,pp,pm,a,agp,agm,adummy
    complex(fp),dimension(size(y0),f%nobs) :: b,bgp,bgm,bpp,bpm
    real(fp),allocatable,dimension(:,:,:) :: rnd
    real(fp),allocatable,dimension(:,:) :: ng
    real(fp),dimension(f%nobs) :: dw
    real(fp),dimension(f%nobs) :: j_avg,dw_sum,j_avg_sum
    logical :: write_next
    integer(ip) :: sint
    type(fisrk15) :: f15
    integer(ip) :: c1,c2,niter,eq_status
    
    !write(*,*) 'Welcome to isrk method. This is not ready yet.'

    n = size(x)
    dim = size(y0)
    
    if(n < 2) then
       write(*,*) 'At least 2 query points must be defined.'
       stop
    end if
    
    if(size(y,2) .NE. n) then
       write(*,*) 'Size of output array is not equal to query points.'
       stop
    end if

    if(present(dx_in)) then
       dx = dx_in
       if(dx > (x(2) - x(1))/2.0_fp) then
          dx = x(2) - x(1)
       end if
    else
       dx=x(2) - x(1)
    endif

    if(present(tol_in)) then
       tol = tol_in
    else
       tol=1e-5_fp
    endif

    if(present(step_split_factor_in)) then
       step_split_factor = step_split_factor_in
    else
       step_split_factor = 0
    end if

    if(step_split_factor > 0) then
       dx = dx/real(2**step_split_factor,fp)
    end if
    y(:,1) = y0
    yn = y0
    ii = 2
    kk = 1
    ll = 1
    niter = 0
    
    ! second order method requires two random numbers per step    
    nrand = 2*(step_split_factor+1)*(ceiling((x(size(x)) - x(1))/dx)+size(x))
    allocate(rnd(f%nobs,nrand,2))
    rnd = 0.0_fp
    allocate(ng(f%nobs,2))
    call rgen%norm(rnd,nrand*f%nobs,0.0_fp,1.0_fp)
    j_avg_sum = 0.0_fp
    j_avg = 0.0_fp
    j = 0.0_fp
    dw_sum = 0.0_fp
    write_next = .false.
    do while(ii <= n)
       write_next = .false.
       xc = x(ii-1)
       xn = x(ii)
       do while(.true.)
          if(xc + dx + dx/1000.0_fp >= xn) then
             dxc = xn - xc 
             write_next = .true.
          else
             dxc = dx
          end if
          sqrtdxc = sqrt(dxc)
          if(step_split_factor>0) then
             ! Idea for splitting can be found in Matti Silveri's notes.
             ! Splitting is a little bit more difficult for processes that require two random numbers, but lets assume that its enough to generate two numbers which independently fill the scaling conditions 
             ! Split once:
             ! dw_tilde_2n = 1/2*dw_n + 1/2*dz_n
             ! dw_n*dw_n = dt, dQ_n*dQ_n = dt
             ! dw_tilde_n*dw_tilde_n = 1/2*dt
             ! Split twice:
             ! dw_hat_4n = 1/4*dw_n+1/4*dz_n+1/2*dQ_n
             ! dQ_n*dQ_n = 1/2 dt
             
             !dw = 0.5_fp*sqrt(dxc)*(rnd(:,kk) + rnd(:,kk))
             !THIS!dw = sqrt(2**(step_split_factor)*dxc)*rnd(:,kk)*0.5_fp**step_split_factor
             ng = sqrt(2.0_fp**(step_split_factor))*rnd(:,kk,:)*0.5_fp**step_split_factor
             !dw = 0.0_fp
             ! The contributions are added from last to first (...,dQ,dZ,dW)
             do mm=0,step_split_factor-1
                sint = 1
                !dw = dw + sqrt(2**(step_split_factor - mm + 1)*dxc)*&
                !rnd(:,kk*(step_split_factor - mm+1))*0.5_fp**(step_split_factor-mm+1)*&
                !(-1)**(rshift(and(ll-1,lshift(sint,mm)),mm))
                !THIS!dw = dw + sqrt(2**(mm + 1)*dxc)*&
                     !rnd(:,kk*(step_split_factor - mm+1))*0.5_fp**(mm+1)*&
                     !(-1)**(rshift(and(ll-1,lshift(sint,mm)),mm))
                ng = ng + sqrt(2.0_fp**(mm+1))*&
                     rnd(:,kk*(step_split_factor -mm+1),:)*0.5_fp**(mm+1)*&
                     (-1)**(rshift(and(ll-1,lshift(sint,mm)),mm))
                !write(*,*) 'll:',ll,', mm:',mm,', ',(-1)**(rshift(and(ll-1,lshift(sint,mm)),mm))
             end do
          else
             !dw = sqrtdxc*rnd(:,kk)
             !dw = rnd(:,kk) ! This is not really dw, but just the random number
             ng = rnd(:,kk,:) ! Inefficient indexing, but necessary
          endif
          
          ! dw is determined, next calculate all the multipliers
          !call system_clock(c1)
          call f15%initialize(xc,yn,ng,f,dxc,f%nobs,j_avg,dw) ! WARNING step splitting not taken into account
          !call system_clock(c2)
          !write(*,*) 'Init: ', c2-c1
          !call system_clock(c1)
          call picard_iteration(f15,yn,yn,tol,20,eq_status) ! input and output are the same, but this should be safe
          if(eq_status .eq. 1) then
             write(*,*) 'Picard iteration failed to converge: ', eq_status
          end if
          !call system_clock(c2)
          !write(*,*) 'Picard: ', c2-c1

          call f%normalize(yn)
          j_avg_sum = j_avg_sum + j_avg*dxc
          dw_sum = dw_sum + dw
          !niter = niter + 1

          xc = xc + dxc
          if(step_split_factor>0) then
             if(ll == 2**step_split_factor) then
                kk = kk+1
                ll = 1
             else
                ll=ll+1
             end if
          else
             kk = kk+1
          end if
          if(write_next) then
             !write(*,*) jj
             y(:,ii) = yn
             j(:,ii) = j_avg_sum/(x(ii)-x(ii-1)) + dw_sum/(x(ii)-x(ii-1))
             exit
          end if
       end do
       ii = ii+1
    end do
    call f15%delete()
    !write(*,*) 'niter:', niter
  end subroutine isrk15

  subroutine isrk15_adaptive(x,y0,f,y,j,rgen,dx_in,step_split_factor_in,tol_in)
    use utils,only:fp,ip
    use sclass
    use random
    use linfuncs
    use fclass
    use equation_solvers
    implicit none
    
    real(fp), dimension(:),intent(in) :: x
    complex(fp), dimension(:),intent(in) :: y0
    complex(fp), dimension(:,:),intent(out) :: y
    real(fp),dimension(:,:),intent(out) :: j
    real(fp),intent(in),optional :: dx_in,tol_in
    class(sd),intent(in) :: f
    class(rng),intent(inout) :: rgen
    integer(ip),intent(in),optional :: step_split_factor_in

    integer(ip) :: n,dim,ii,jj,kk,ll,nrand,mm,step_split_factor,ss
    real(fp) :: dx,xc,xn,dxc,sqrtdxc,tol
    complex(fp),dimension(size(y0)) :: yn,gp,gm,pp,pm,a,agp,agm,adummy
    complex(fp),dimension(size(y0)) :: yn_t
    complex(fp),dimension(size(y0),f%nobs) :: b,bgp,bgm,bpp,bpm
    real(fp),allocatable,dimension(:,:,:) :: rnd
    real(fp),allocatable,dimension(:,:) :: ng
    real(fp),dimension(f%nobs) :: dw
    logical :: write_next
    integer(ip) :: sint
    type(fisrk15) :: f15
    integer(ip) :: c1,c2,niter,eq_status
    real(fp),dimension(f%nobs) :: j_avg,dw_sum,j_avg_sum
    
    !write(*,*) 'Welcome to isrk method. This is not ready yet.'

    n = size(x)
    dim = size(y0)
    
    if(n < 2) then
       write(*,*) 'At least 2 query points must be defined.'
       stop
    end if
    
    if(size(y,2) .NE. n) then
       write(*,*) 'Size of output array is not equal to query points.'
       stop
    end if

    if(present(dx_in)) then
       dx = dx_in
       if(dx > (x(2) - x(1))/2.0_fp) then
          dx = x(2) - x(1)
       end if
    else
       dx=x(2) - x(1)
    endif

    if(present(tol_in)) then
       tol = tol_in
    else
       tol=1e-5_fp
    endif

    if(present(step_split_factor_in)) then
       step_split_factor = step_split_factor_in
    else
       step_split_factor = 0
    end if
    if(step_split_factor > 0) then
       write(*,*) 'Warning, isrk15_adaptive does not use step splitting.'
    endif

    y(:,1) = y0
    yn = y0
    ii = 2
    kk = 1
    ll = 1
    niter = 0
    ss = 0
    ! second order method requires two random numbers per step
    
    nrand = 2*(ceiling((x(size(x)) - x(1))/dx)+size(x))

    allocate(rnd(f%nobs,nrand,2))
    rnd = 0.0_fp
    allocate(ng(f%nobs,2))
    call rgen%norm(rnd,nrand*f%nobs,0.0_fp,1.0_fp)
    j_avg_sum = 0.0_fp
    j_avg = 0.0_fp
    j = 0.0_fp
    dw_sum = 0.0_fp
    write_next = .false.
    do while(ii <= n)
       write_next = .false.
       xc = x(ii-1)
       xn = x(ii)
       do while(.true.)
          if(xc + dx + dx/1000.0_fp >= xn) then
             dxc = xn - xc 
             write_next = .true.
          else
             dxc = dx
          end if
          sqrtdxc = sqrt(dxc)
             !dw = sqrtdxc*rnd(:,kk)
             !dw = rnd(:,kk) ! This is not really dw, but just the random number
          ! Create more random numbers if needed due to the reduced number of steps.
          if(2*kk > nrand) then
             call rgen%norm(rnd,nrand*f%nobs,0.0_fp,1.0_fp)
             kk = 1
          endif
          ng = rnd(:,kk,:) ! Inefficient indexing, but necessary
          
          ! dw is determined, next calculate all the multipliers
          !call system_clock(c1)
          call f15%initialize(xc,yn,ng,f,dxc,f%nobs,j_avg,dw) ! WARNING step splitting not taken into account
          !call system_clock(c2)
          !write(*,*) 'Init: ', c2-c1
          !call system_clock(c1)
          call picard_iteration(f15,yn,yn_t,tol,20,eq_status) ! We need yn_t to hold the iterated yn, because the result might be rejected.
          
          if(eq_status .eq. 1 .and. ss < 20) then ! reject
             dx = 0.75*dx
             write(*,*) 'Adjusting the step size. The new size is ', dx,'.'
             write_next = .false.
             xc = x(ii-1) ! Return to the point where last value was written
             yn = y(:,ii-1) ! For some reason yn is written during the iteration, therefore make sure its value is reset. THIS IS A BUG.
             ! write(*,*) 'x0=',x(1),'xc = ',xc,'ii =',ii
             ! write(*,*) 'y0 = ',y0
             ! write(*,*) 'yn = ',yn
             ! write(*,*) 'yn_t = ',yn_t
             ! write(*,*) 'y(:,ii-1)',y(:,ii-1)
             ! call sleep(1)
             ss = ss+1
          else ! accept
             yn = yn_t
             call f%normalize(yn)
             j_avg_sum = j_avg_sum + j_avg*dxc
             dw_sum = dw_sum + dw
             xc = xc + dxc
             kk = kk + 1
             ss = 0
          end if
          !write(*,*) 'status: ', eq_status
          !call system_clock(c2)
          !write(*,*) 'Picard: ', c2-c1
          if(write_next) then
             !write(*,*) jj
             y(:,ii) = yn
             j(:,ii) = j_avg_sum/(x(ii)-x(ii-1)) + dw_sum/(x(ii)-x(ii-1))
             exit
          end if
       end do
       ii = ii+1
    end do
    call f15%delete()
    !write(*,*) 'niter:', niter
  end subroutine isrk15_adaptive

end module sde_solvers
   
! program test_sde_solvers
!   use utils,only: fp,ip
!   use sclass
!   use hamiltonians
!   use random
!   use linfuncs
!   use ode_s_solvers
!   use sde_solvers
!   use envelopes
!   use signals
!   implicit none
!   integer(ip) :: nql,ncl
!   type(djch) :: Hdjc
!   type(lcg) :: lcg_t
!   type(diff_schr) :: sfun
!   real(fp),allocatable,dimension(:) :: wq
!   complex(fp),allocatable,dimension(:) :: y0s
!   integer(ip) :: ii,jj,kk,NN
!   complex(fp),allocatable,dimension(:) :: d
!   complex(fp),allocatable,dimension(:,:) :: s
!   real(fp),allocatable,dimension(:) :: x,xquery
!   complex(fp),allocatable,dimension(:,:) :: y,yquery
!   complex(fp),allocatable,dimension(:,:) :: y_step_split
!   real(fp),allocatable,dimension(:,:) :: j
!   type(signal) :: test_signal
!   type(constant_envelope) :: cenv
!   integer(ip) :: c1,c2

!   nql = 2
!   ncl = 1
!   NN = 501

!   allocate(x(NN))
!   call linspace(0.0_fp,300.0_fp,NN,x)

!   allocate(y0s(nql*ncl))
!   y0s = (0.0_fp,0.0_fp)
!   y0s(2) = 1.0_fp

!   allocate(wq(nql))
!   wq = 0.0_fp
!   do ii=1,size(wq)
!      wq(ii) = (ii-1)**0.7
!   end do

!   allocate(d(nql*ncl))
!   allocate(s(nql*ncl,2))

!   allocate(y(nql*ncl,NN))
!   allocate(j(1,NN))

!   call Hdjc%initialize_djch(nql,ncl,wq,0.8_fp,0.01_fp,hbar_in = 1.0_fp)
!   call test_signal%initialize(Hdjc%op%cc(:,:,1),1.0_fp)
!   call cenv%initialize_constant_envelope((0.05_fp,0.0_fp))
!   call test_signal%add_envelope(cenv)
!   call Hdjc%add_signal(test_signal)
!   !write(*,*) Hdjc%value_at(0.0_fp)
!   call sfun%init_diff_schr(Hdjc)
!   call sfun%add_observable(Hdjc%op%aa(:,:,1),(0.0_fp,0.02_fp))
!   !call sfun%value_at(0.0_fp,y0s,d,s)
!   !call lcg_t%init(generate_random_seed())
!   call lcg_t%init(0)
!   !call test_sde_solver(x,y0s,sfun,y,j,lcg_t,0.01_fp)

!   call platen(x,y0s,sfun,y,j,lcg_t,0.2_fp)
!   !call euler_maryama(x,y0s,sfun,y,j,lcg_t,0.01_fp)

!   allocate(xquery(NN/3),yquery(size(y0s),NN/3))
!   call linspace(0.1_fp,290.5_fp,NN/3,xquery)
!   yquery = interp1(x,y,xquery)
  
!   !do ii=1,size(x)
!   !   y(:,ii) = y(:,ii)/l2norm(y(:,ii))
!   !end do
  
!   open(unit = 1, file = 'st_data_s.dat')
!   open(unit = 2, file = 'st_data_s_interp.dat')
!   do ii=1,size(x)
!      write(1,*) x(ii), real(conjg(y(:,ii))*y(:,ii))
!   end do
!   do ii=1,size(xquery)
!      write(2,*) xquery(ii),real(conjg(yquery(:,ii))*yquery(:,ii))
!   end do
!   close(1)
!   close(2)

!   call system_clock(c1)
!   !call basic_ode_s_solver(x,y0s,sfun,y,euler_maryama,100,0.1_fp)
!   !call basic_ode_s_solver(x,y0s,sfun,y,platen,1000,0.2_fp)
!   call system_clock(c2)
!   write(*,*) 'Linear', (c2-c1), 'milliseconds'
!   call system_clock(c1)
!   !call basic_ode_s_solver(x,y0s,sfun,y,euler_maryama,100,0.1_fp)
!   !call basic_ode_s_solver(x,y0s,sfun,y,platen,100,2.0_fp)
!   call parallel_ode_s_solver(x,y0s,sfun,y,platen,40,0.2_fp)
!   call system_clock(c2)
!   write(*,*) 'Parallel', (c2-c1), 'milliseconds'


!   open(unit = 1, file = 'st_data_ode_s.dat')
!   do ii=1,size(x)
!      write(1,*) x(ii), real(y(:,ii))
!   end do
!   close(1)

!   allocate(y_step_split(NN,4))
!   do ii=1,4
!      call lcg_t%init(0)
!      call platen(x,y0s,sfun,y,j,lcg_t,1.0_fp,2*(ii-1))
!      y_step_split(:,ii) = y(1,:)
!   end do
  
!   open(unit = 1, file = 'step_split.dat')
!   do ii=1,size(x)
!      write(1,*) x(ii), real(conjg(y_step_split(ii,:))*y_step_split(ii,:))
!   end do
!   close(1)

! end program test_sde_solvers
