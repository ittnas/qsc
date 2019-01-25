program sparse_blas_test
use utils,only: fp,ip
use blas_sparse
implicit none
integer(ip) :: a,istat,ii
real(fp),dimension(5) :: val
integer(ip),dimension(5) :: idx,jdx
real(fp),dimension(4) :: x,y

val = [1.0_fp,2.0_fp,3.0_fp,4.0_fp,5.0_fp]
idx = [1, 2, 3, 4, 5]
jdx = [2, 3, 4, 5, 6]
x = [1.0_fp, 1.0_fp, 1.0_fp, 1.0_fp]
y = 0.0_fp


call duscr_begin(4,4,a,istat)
do ii=1,5
   call duscr_insert_entry(a,val(ii),idx(ii),jdx(ii),istat)
end do
call uscr_end(a,istat)

call usmv(a,x,y,istat)
call usds(a,istat)
write(*,*) y

end program sparse_blas_test
