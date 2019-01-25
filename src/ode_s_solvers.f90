module ode_s_solvers
  use utils,only: fp,ip
  use hamiltonians ! Adding this makes everything compile
  use linfuncs
  
  abstract interface
     subroutine ode_s_solver(x,y0,f,y,solver,N,dx_in,seed_in,tol_in)
       use utils,only:fp,ip
       use sclass
       use sde_solvers
       implicit none
       real(fp), dimension(:),intent(in) :: x
       complex(fp), dimension(:),intent(in) :: y0
       complex(fp), dimension(:,:),intent(out) :: y
       real(fp),intent(in),optional :: dx_in,tol_in
       class(sd),intent(in) :: f
       procedure(sde_solver) :: solver
       integer(ip), intent(in) :: N
       integer(ip), intent(in),optional :: seed_in
     end subroutine ode_s_solver
  end interface
contains
  subroutine basic_ode_s_solver(x,y0,f,y,solver,N,dx_in,seed_in,tol_in)
    use utils,only:fp,ip
    use hamiltonians
    use random
    use sclass
    use sde_solvers
    implicit none
    real(fp), dimension(:),intent(in) :: x
    complex(fp), dimension(:),intent(in) :: y0
    complex(fp), dimension(:,:),intent(out) :: y
    real(fp),intent(in),optional :: dx_in,tol_in
    class(sd),intent(in) :: f
    procedure(sde_solver) :: solver
    integer(ip), intent(in) :: N
    integer(ip), intent(in),optional :: seed_in
    integer(ip) :: ii,dim,jj,seed
    real(fp) :: dx,tol
    complex(fp), dimension(size(y,1),size(y,2)),target :: y_temp
    complex(fp), dimension(:,:,:),pointer :: y_reshaped,y_rs_temp
    real(fp),dimension(f%nobs,size(x)) :: j
    type(lcg) :: lcg_t

    if(size(x) < 2) then
       write(*,*) 'At least 2 query points must be defined.'
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
    if(present(seed_in)) then
       seed = seed_in
    else
       seed = generate_random_seed()
    end if

    if(present(tol_in)) then
       tol = tol_in
    else
       tol=1e-5_fp
    endif


    !dim = sqrt(real(size(y,1),fp))
    !if(dim*dim .ne. size(y,1)) then
    !write(*,*) 'y must be a density matrix.'
    !stop
 !end if
    call lcg_t%init(seed)
    y = (0.0_fp,0.0_fp)
    do ii=1,N
       call solver(x,y0,f,y_temp,j,lcg_t,dx)
       !y_rs_temp(1:dim,1:dim,1:size(x)) => y_temp
       !do jj=1,dim
       !y_reshaped(:
    !end do
       !y_reshaped = y_reshaped + 
       !y = y + y_temp
       y = y + conjg(y_temp)*y_temp
    end do
    y = y/N
  end subroutine basic_ode_s_solver
!###################################################
  subroutine parallel_ode_s_solver(x,y0,f,y,solver,N,dx_in,seed_in,tol_in)
    use utils,only:fp,ip
    use hamiltonians
    use random
    use sclass
    use sde_solvers
    use omp_lib
    implicit none
    real(fp), dimension(:),intent(in) :: x
    complex(fp), dimension(:),intent(in) :: y0
    complex(fp), dimension(:,:),intent(out) :: y
    real(fp),intent(in),optional :: dx_in,tol_in
    class(sd),intent(in) :: f
    procedure(sde_solver) :: solver
    integer(ip), intent(in) :: N
    integer(ip), intent(in),optional :: seed_in
    integer(ip) :: ii,dim,jj
    real(fp) :: dx,tol
    complex(fp), dimension(size(y0),size(y,2)),target :: y_temp
    complex(fp), dimension(size(y0)) :: y_vec
    !complex(fp), dimension(:,:,:),pointer :: y_reshaped,y_rs_temp
    real(fp),dimension(f%nobs,size(x)) :: j
    type(lcg) :: lcg_t
    integer(ip) :: seed
    if(size(x) < 2) then
       write(*,*) 'At least 2 query points must be defined.'
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
    
    if(present(seed_in)) then
       seed = seed_in
    else
       seed = generate_random_seed()
    end if

    if(present(tol_in)) then
       tol = tol_in
    else
       tol=1e-5_fp
    endif
    
    y = (0.0_fp,0.0_fp)
    !---------------------------------------------------------------------------------------
    !$OMP PARALLEL DEFAULT(shared),&
    !$OMP PRIVATE(ii,y_temp,j,lcg_t,y_vec),&
    !$OMP SHARED(y,N,x,y0,f,dx,seed,tol)
    call lcg_t%init(lcg_t%get_parallel_seed(seed,omp_get_thread_num()))
    !$OMP DO
    !---------------------------------------------------------------------------------------
    do ii=1,N
       y_temp = (0.0_fp,0.0_fp)
       call solver(x,y0,f,y_temp,j,lcg_t,dx,0,tol)
       !$OMP CRITICAL
       ! That's not how it should work. Y only gives the diagonal elements, but hides the non-diagonal.
       do jj=1,size(y,2)
          y_vec = y_temp(:,jj)
          y(:,jj) = y(:,jj) + reshape(matmul(conjg(reshape(y_vec,[size(y_vec),1])),reshape(y_vec,[1,size(y_vec)])),[size(y,1)])
       end do
       !y = y + conjg(y_temp)*y_temp
       !call print_matrix(matmul(conjg(reshape(y_temp,[size(y_temp),1])),reshape(y_temp,[1,size(y_temp)])))
       !$OMP END CRITICAL
    end do
    !---------------------------------------------------------------------------------------    
    !$OMP END DO
    !$OMP END PARALLEL
    !---------------------------------------------------------------------------------------
    y = y/N
  end subroutine parallel_ode_s_solver
end module ode_s_solvers
