program testzvode
  use iso_fortran_env, only: dp => real64
  implicit none
  external fex, jex
  complex(dp), parameter :: im = (0_dp, 1_dp)
  integer :: neq, itol, itask, istate, iopt, lzw, lrw, liw, mf, iout, ipar
  integer :: iwork(32)
  real(dp) :: t, dtout, tout, rtol, atol, aemax
  complex(dp) :: y(2), zwork(24), rpar, err, wtru
  real(dp), dimension(22) :: rwork

  neq = 2
  y(1) = complex(1.0_dp, 0.0_dp) / 2.1_dp
  y(2) = complex(1.0_dp, 0.0_dp)
  t = 0.0_dp
  dtout = 0.1570796326794896_dp
  tout = dtout
  itol = 1
  rtol = 1.0e-9_dp
  atol = 1.0e-8_dp
  itask = 1
  istate = 1
  iopt = 0
  lzw = 24
  lrw = 22
  liw = 32
  mf = 21
  rpar = complex(0.0_dp,1.0_dp)
  aemax = 0.0_dp
  ipar = 0

  write(*, 10)
  10 format('   t', 11x, 'w', 26x, 'z')

  do iout = 1, 40
    call zvode(fex, neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, &
          zwork, lzw, rwork, lrw, iwork, liw, jex, mf, rpar, ipar)
    wtru = 1_dp/complex(cos(t) + 1.1_dp, sin(t))
    err = y(1) - wtru
    aemax = max(aemax, abs(err))
    write(*, 20) t, real(y(1)), aimag(y(1)), real(y(2)), aimag(y(2))
  20 format(f9.5, 2x, 2f12.7, 3x, 2f12.7)

    if (istate < 0) then
      write(*, 30) istate
   30 format(//'***** Error halt.  istate =', i3)
      stop
    endif

  tout = tout + dtout
  end do
  write(*, 50) iwork(11), iwork(12), iwork(13), iwork(20), &
            iwork(21), iwork(22), iwork(23), aemax
  50 format(/' No. steps =', i4, '   No. f-s =', i5, &
        '   No. J-s =', i4, '   No. LU-s =', i4, /, &
        ' No. nonlinear iterations =', i4, /, &
        ' No. nonlinear convergence failures =', i4, /, &
        ' No. error test failures =', i4, /, &
        ' Max. abs. error in w =', d10.2)


end program testzvode

subroutine fex(neq, t, y, ydot, rpar, ipar)
  use iso_fortran_env, only: dp => real64
  integer, intent(in) :: neq
  real(dp), intent(in) :: t
  complex(dp), dimension(neq), intent(in) :: y
  complex(dp), dimension(neq), intent(out) :: ydot
  complex(dp), intent(in) :: rpar
  integer, intent(in) :: ipar
  ydot(1) = -rpar * y(1) * y(1) * y(2)
  ydot(2) = rpar * y(2)
end subroutine fex

subroutine jex(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
  use iso_fortran_env, only: dp => real64
  integer, intent(in) :: neq, ml, mu, nrpd
  real(dp), intent(in) :: t
  complex(dp), dimension(neq), intent(in) :: y
  complex(dp), dimension(nrpd, neq), intent(out) :: pd
  complex(dp), intent(in) :: rpar
  integer, intent(in) :: ipar
  pd(1, 1) = -2.0_dp * rpar * y(1) * y(2)
  pd(1, 2) = -rpar * y(1) * y(1)
  pd(2, 2) = rpar
end subroutine jex

