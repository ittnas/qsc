program ptrtest
   real, pointer :: a(:)
   real, pointer :: b(:,:)
   integer :: n = 10
   allocate(a(n**2))
   a = 42
   b (1:n, 1:n) => a

   WRITE(*,*) a
end program ptrtest
