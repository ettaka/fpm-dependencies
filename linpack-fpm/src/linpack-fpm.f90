module linpack_fpm
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, linpack-fpm!"
  end subroutine say_hello
end module linpack_fpm
