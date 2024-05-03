module lapack_fpm
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, lapack-fpm!"
  end subroutine say_hello
end module lapack_fpm
