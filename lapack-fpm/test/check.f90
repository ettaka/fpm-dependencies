program check
implicit none

print *, "Put some tests in here!"

call test_gemvrsp()

contains


  subroutine test_gemvrsp()
      real :: A(3,3),x(3),y(3),ylap(3),yintr(3),alpha,beta
      call random_number(alpha)
      call random_number(beta)
      call random_number(A)
      call random_number(x)
      call random_number(y)
      ylap = y
      call sgemv('No transpose',size(A,1),size(A,2),alpha,A,size(A,1),x,1,beta,ylap,1)
      yintr = alpha*matmul(A,x)+beta*y

      print *, "A="
      print '(3(f3.1,x))', A

      print '("yintr = ",3(f8.5,x))', yintr
      print '("ylap = ",3(f8.5,x))', ylap

  end subroutine test_gemvrsp

end program check
