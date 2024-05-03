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

module zvode_test_program

  contains

    subroutine testzvode
      use iso_fortran_env, only: dp => real64, stdout => output_unit, stderr => error_unit
      implicit none
      external fex, jex
      integer, parameter :: neq = 2                  ! number of equations
      integer, parameter :: mf = 21                  ! method: 21 for stiff (BDF) method, user-supplied full Jacobian.
      integer, parameter :: lzw = 8*neq + 2*neq**2   ! zwork dimension (different for each method) 
      integer, parameter :: lrw = 20 + neq           ! rwork array length
      integer, parameter :: liw = 30 + neq           ! integer work array length (different for mf = 10)
      integer, parameter :: itol = 1                 ! estimated local error: EWT(i) = RTOL*abs(Y(i)) + ATOL
      integer, parameter :: itask = 1                ! set normal computation of output values of y(t) at t = TOUT
                                                     ! (by overshooting and interpolating)
      integer, parameter :: iopt = 0                 ! no optional input used
      real(kind=dp), parameter :: rtol = 1.0e-9_dp   ! relative tolerance
      real(kind=dp), parameter :: atol = 1.0e-8_dp   ! absolute tolerance
        
      integer :: istate, iout, ipar
      integer :: iwork(liw)
      real(dp) :: t, time_step, tout, aemax
      complex(dp) :: y(neq), zwork(lzw), rpar, err, wtru
      real(dp) :: rwork(lrw)
      character(*), parameter :: headerfmt = "('   t', 11x, 'w', 26x, 'z')"
      character(*), parameter :: resfmt = "(f9.5, 2x, 2f12.7, 3x, 2f12.7)"
      character(*), parameter :: errorfmt = "(//'***** Error halt.  istate =', i3)"
      character(*), parameter :: summaryfmt = &
                "(/' No. steps =', i4, '   No. f-s =', i5, '   No. J-s =', i4, '   No. LU-s =', i4, /, &
                 & ' No. nonlinear iterations =', i4, /,                                               &
                 & ' No. nonlinear convergence failures =', i4, /,                                     &
                 & ' No. error test failures =', i4, /,                                                &
                 & ' Max. abs. error in w =', d10.2)"

      y(1) = complex(1.0_dp, 0.0_dp) / 2.1_dp       ! initial value for y(1)
      y(2) = complex(1.0_dp, 0.0_dp)                ! initial value for y(2)
      t = 0.0_dp                                    ! initial time
      time_step = 0.1570796326794896_dp             ! time step
      tout = time_step                              ! time after one step
      istate = 1                                    ! set to 1 at first (meaning: first call to zvode)
      rpar = complex(0.0_dp,1.0_dp)                 ! user-defined parameter passed to f and jacobian functions
      ipar = 0                                      ! user-defined parameter passed to f and jacobian functions

      aemax = 0.0_dp                                ! initialize maximum error (computed later in this routine)

      write(stdout, headerfmt)
      time_loop: do iout = 1, 40
        !-----------------------------------------------------------------------
        ! `zvode` is a Variable-coefficient Ordinary Differential Equation (ODE)
        ! solver, with fixed-leading-coefficient implementation. It provides
        ! implicit Adams method (for non-stiff problems) and a method based on
        ! backward differentiation formulas (BDF) (for stiff problems).
        ! The ODE is of the form:
        !
        !     dy/dt = f(t,y) ,  or, in component form,
        !     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq),
        !
        ! Parameters to the subroutine are (see zvode.f:99):
        ! f      = Name of subroutine for right-hand side vector f.
        !          This name must be declared external in calling program.
        ! neq    = Number of first order ODEs.
        ! y      = Double complex array of initial values, of length neq.
        ! t      = The initial value of the independent variable.
        ! tout   = First point where output is desired (.ne. t).
        ! itol   = 1 or 2 according as ATOL (below) is a scalar or array.
        ! rtol   = Relative tolerance parameter (scalar).
        ! atol   = Absolute tolerance parameter (scalar or array).
        !          The estimated local error in Y(i) will be controlled so as
        !          to be roughly less (in magnitude) than
        !             EWT(i) = RTOL*abs(Y(i)) + ATOL     if ITOL = 1, or
        !             EWT(i) = RTOL*abs(Y(i)) + ATOL(i)  if ITOL = 2.
        !          Thus the local error test passes if, in each component,
        !          either the absolute error is less than ATOL (or ATOL(i)),
        !          or the relative error is less than RTOL.
        !          Use RTOL = 0.0 for pure absolute error control, and
        !          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
        !          control.  Caution: Actual (global) errors may exceed these
        !          local tolerances, so choose them conservatively.
        ! itask  = 1 for normal computation of output values of Y at t = tout.
        ! istate = Integer flag (input and output).  Set istate = 1.
        ! iopt   = 0 to indicate no optional input used.
        ! zwork  = Double precision complex work array of length at least:
        !             15*neq                      for mf = 10,
        !             8*neq + 2*neq**2            for mf = 21 or 22,
        !             10*neq + (3*ml + 2*mu)*neq  for mf = 24 or 25.
        ! lzw    = Declared length of zwork (in user's dimension statement).
        ! rwork  = Real work array of length at least 20 + neq.
        ! lrw    = Declared length of rwork (in user's dimension statement).
        ! iwork  = Integer work array of length at least:
        !             30        for mf = 10,
        !             30 + neq  for mf = 21, 22, 24, or 25.
        !          If mf = 24 or 25, input in iwork(1),iwork(2) the lower
        !          and upper half-bandwidths ml,mu.
        ! liw    = Declared length of iwork (in user's dimension statement).
        ! jac    = Name of subroutine for Jacobian matrix (mf = 21 or 24).
        !          If used, this name must be declared external in calling
        !          program.  If not used, pass a dummy name.
        ! mf     = Method flag.  Standard values are:
        !          10 for nonstiff (Adams) method, no Jacobian used.
        !          21 for stiff (BDF) method, user-supplied full Jacobian.
        !          22 for stiff method, internally generated full Jacobian.
        !          24 for stiff method, user-supplied banded Jacobian.
        !          25 for stiff method, internally generated banded Jacobian.
        ! rpar   = user-defined real or complex array passed to f and jac.
        ! ipar   = user-defined integer array passed to f and jac.
        ! Note that the main program must declare arrays y, zwork, rwork, iwork,
        ! and possibly atol, rpar, and ipar.  rpar may be declared real, double,
        ! complex, or double complex, depending on the user's needs.
        call zvode(fex, neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, &
                   zwork, lzw, rwork, lrw, iwork, liw, jex, mf, rpar, ipar)
        wtru = 1_dp/complex(cos(t) + 1.1_dp, sin(t))
        err = y(1) - wtru
        aemax = max(aemax, abs(err))

        write(stdout, resfmt) t, y

        if (istate < 0) then
          write(stderr, errorfmt) istate
          exit time_loop
        endif

        tout = tout + time_step
      end do time_loop

      write(stdout, summaryfmt) iwork(11:13), iwork(20:23), aemax
    end subroutine testzvode

end module zvode_test_program
