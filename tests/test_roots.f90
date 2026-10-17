!> Assertion-based test suite for the bracket_root module.
!>
!> Exits with a nonzero status (error stop 1) if any check fails, so it can
!> gate CI via `make test`.
module test_functions
   use bracket_root, only: rk
   implicit none
contains

   pure real(rk) function f_cubic(x)      ! classic x^3 - 2x - 5 = 0
      real(rk), intent(in) :: x
      f_cubic = x**3 - 2.0_rk*x - 5.0_rk
   end function f_cubic

   pure real(rk) function f_cosx(x)       ! cos x = x
      real(rk), intent(in) :: x
      f_cosx = cos(x) - x
   end function f_cosx

   pure real(rk) function f_sinx(x)       ! sin x = 1/2, analytic root pi/6
      real(rk), intent(in) :: x
      f_sinx = sin(x) - 0.5_rk
   end function f_sinx

   pure real(rk) function f_omega(x)      ! x e^x = 1, root = W(1)
      real(rk), intent(in) :: x
      f_omega = x*exp(x) - 1.0_rk
   end function f_omega

   pure real(rk) function f_samesign(x)   ! strictly positive: no sign change
      real(rk), intent(in) :: x
      f_samesign = x*x + 1.0_rk
   end function f_samesign

   pure real(rk) function f_domainerr(x)  ! NaN for x < 4 (sqrt of negative)
      real(rk), intent(in) :: x
      f_domainerr = sqrt(x - 4.0_rk) - 1.0_rk
   end function f_domainerr

   pure real(rk) function f_exp2(x)       ! convex; hard for plain regula falsi
      real(rk), intent(in) :: x
      f_exp2 = exp(x) - 2.0_rk
   end function f_exp2

   pure real(rk) function f_step(x)       ! sign jump at x = 1, no finite zero
      real(rk), intent(in) :: x
      if (x < 1.0_rk) then
         f_step = -1.0_rk
      else
         f_step = 1.0_rk
      end if
   end function f_step

   pure real(rk) function f_line(x)       ! root exactly at 0
      real(rk), intent(in) :: x
      f_line = x
   end function f_line

   pure real(rk) function f_asym(x)       ! root 2.01, near the left end of [2,3]
      real(rk), intent(in) :: x
      f_asym = x - 2.01_rk
   end function f_asym

end module test_functions

program test_roots
   use, intrinsic :: iso_fortran_env, only: output_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use bracket_root
   use test_functions
   implicit none

   integer :: nchecks, nfail, i
   type(root_result) :: rb, ri, rs, rr
   type(step_record) :: hist(25)
   real(rk), parameter :: r_cubic = 2.0945514815423265_rk
   real(rk), parameter :: r_cosx  = 0.7390851332151607_rk
   real(rk), parameter :: r_sinx  = 0.5235987755982989_rk   ! pi/6
   real(rk), parameter :: r_omega = 0.56714329040978387_rk  ! W(1)
   real(rk) :: ratio, nan_val

   nchecks = 0
   nfail = 0
   nan_val = ieee_value(1.0_rk, ieee_quiet_nan)

   write (output_unit, '(a)') '== bracket_root test suite =='

   ! ---- 1. known roots, three solvers, four equations -------------------
   call solve3('cubic', f_cubic, 2.0_rk, 3.0_rk, r_cubic)
   call solve3('cosx ', f_cosx, 0.0_rk, 1.0_rk, r_cosx)
   call solve3('sinx ', f_sinx, 0.0_rk, 1.5_rk, r_sinx)
   call solve3('omega', f_omega, 0.0_rk, 2.0_rk, r_omega)

   ! ---- 2. invalid brackets ----------------------------------------------
   call invalid_case('same-sign endpoints', f_samesign, 0.0_rk, 1.0_rk)
   call invalid_case('reversed bracket    ', f_cubic, 3.0_rk, 2.0_rk)
   call invalid_case('non-finite f(a)     ', f_domainerr, 0.0_rk, 9.0_rk)
   call invalid_case('NaN endpoint        ', f_cubic, nan_val, 3.0_rk)

   ! ---- 3. tiny iteration budget -> MAX_ITERATIONS ------------------------
   call bisection(f_cubic, 2.0_rk, 3.0_rk, 0.0_rk, 0.0_rk, 3, rb)
   call check(rb%flag == root_err_max_iterations, 'tiny budget: MAX_ITERATIONS flag')
   call check(rb%iterations == 3, 'tiny budget: iterations == 3')
   call check(rb%root >= 2.0_rk .and. rb%root <= 3.0_rk, 'tiny budget: iterate within bracket')
   call safeguarded_secant(f_cubic, 2.0_rk, 3.0_rk, 0.0_rk, 0.0_rk, 3, rs)
   call check(rs%flag == root_err_max_iterations .and. rs%iterations == 3, &
              'safeguarded tiny budget: MAX_ITERATIONS flag, iterations == 3')

   ! ---- 4. zero iteration budget ------------------------------------------
   call bisection(f_cubic, 2.0_rk, 3.0_rk, 1.0e-14_rk, 1.0e-14_rk, 0, rb)
   call check(rb%flag == root_err_max_iterations .and. rb%iterations == 0, &
              'zero budget: MAX_ITERATIONS flag with 0 iterations')

   ! ---- 5. machine-precision stagnation -> NOT_CONVERGED ------------------
   call bisection(f_step, 0.5_rk, 1.5_rk, 0.0_rk, 1.0e-12_rk, 100, rb)
   call check(rb%flag == root_err_not_converged, 'stagnating bisection: NOT_CONVERGED flag')
   call safeguarded_secant(f_step, 0.5_rk, 1.5_rk, 0.0_rk, 1.0e-12_rk, 100, rs)
   call check(rs%flag == root_err_not_converged, 'stagnating safeguarded: NOT_CONVERGED flag')

   ! ---- 6. bisection contraction law + history invariants -----------------
   call bisection(f_cubic, 2.0_rk, 3.0_rk, 0.0_rk, 0.0_rk, 25, rb, hist)
   call check(rb%flag == root_err_max_iterations .and. rb%iterations == 25, &
              'contraction run: 25 iterations used, MAX_ITERATIONS flag')
   do i = 1, 25
      call check(hist(i)%a <= hist(i)%x .and. hist(i)%x <= hist(i)%b, &
                 'history: iterate lies inside its bracket')
   end do
   do i = 2, 25
      ratio = hist(i)%halfwidth / hist(i-1)%halfwidth
      call check(abs(ratio - 0.5_rk) <= 1.0e-9_rk, &
                 'history: halfwidth contracts by factor 1/2 per iteration')
   end do

   ! ---- 7. Illinois vs plain regula falsi on a hard convex case -----------
   call regula_falsi(f_exp2, 0.0_rk, 5.0_rk, 1.0e-12_rk, 1.0e-12_rk, 1000, rr)
   call illinois(f_exp2, 0.0_rk, 5.0_rk, 1.0e-12_rk, 1.0e-12_rk, 200, ri)
   call safeguarded_secant(f_exp2, 0.0_rk, 5.0_rk, 1.0e-12_rk, 1.0e-12_rk, 200, rs)
   call check(rr%flag == root_ok, 'regula falsi, hard convex: eventually converges (flag OK)')
   call check(ri%flag == root_ok, 'illinois, hard convex: flag OK')
   call check(rs%flag == root_ok, 'safeguarded, hard convex: flag OK')
   call check(abs(ri%root - log(2.0_rk)) <= 1.0e-12_rk, 'illinois, hard convex: root == ln(2)')
   call check(ri%iterations < rr%iterations, 'illinois needs fewer iterations than regula falsi')
   call check(3*ri%iterations <= rr%iterations, &
              'regula falsi one-sided stagnation: >= 3x illinois iterations')
   write (output_unit, '(a,i0,a,i0)') 'note: illinois iterations = ', ri%iterations, &
                                      '  regula falsi iterations = ', rr%iterations

   ! ---- 8. root exactly at an endpoint -------------------------------------
   call bisection(f_line, -1.0_rk, 0.0_rk, 1.0e-14_rk, 1.0e-14_rk, 100, rb)
   call check(rb%flag == root_ok .and. rb%iterations == 0 .and. abs(rb%root) <= 0.0_rk, &
              'endpoint root: immediate OK, 0 iterations')

   ! ---- 9. flag_message mapping --------------------------------------------
   call check(trim(flag_message(root_ok)) == 'OK', 'flag_message: OK')
   call check(trim(flag_message(root_err_invalid_bracket)) == 'INVALID_BRACKET', &
              'flag_message: INVALID_BRACKET')
   call check(trim(flag_message(root_err_max_iterations)) == 'MAX_ITERATIONS', &
              'flag_message: MAX_ITERATIONS')
   call check(trim(flag_message(root_err_not_converged)) == 'NOT_CONVERGED', &
              'flag_message: NOT_CONVERGED')

   ! ---- 10. xtol contract: width-stop returns the bracket midpoint --------
   ! f(x) = x - 2.01 on [2,3], xtol = 0.1, ftol = 0: bisection stops after
   ! 3 iterations with bracket [2, 2.125].  Returning the probed endpoint
   ! 2.125 would give |root - r| = 0.115 > xtol; the midpoint 2.0625 gives
   ! 0.0525 <= xtol.  (Regression test for the midpoint-return fix.)
   call bisection(f_asym, 2.0_rk, 3.0_rk, 0.1_rk, 0.0_rk, 100, rb)
   call check(rb%flag == root_ok, 'xtol contract: loose-xtol bisection returns OK')
   call check(rb%iterations == 3, 'xtol contract: stops at first half-width <= xtol')
   call check(abs(rb%root - 2.0625_rk) <= 1.0e-12_rk, &
              'xtol contract: width-stop returns final bracket midpoint 2.0625')
   call check(abs(rb%root - 2.01_rk) <= 0.1_rk, &
              'xtol contract: |root - r| <= xtol (endpoint return would give 0.115)')
   call check(abs(rb%residual - (2.0625_rk - 2.01_rk)) <= 1.0e-12_rk, &
              'xtol contract: reported residual is the true midpoint residual')
   call illinois(f_asym, 2.0_rk, 3.0_rk, 0.1_rk, 0.0_rk, 100, ri)
   call check(ri%flag == root_ok .and. abs(ri%root - 2.01_rk) <= 0.1_rk + 1.0e-12_rk, &
              'xtol contract: illinois midpoint return within xtol of root')
   call safeguarded_secant(f_asym, 2.0_rk, 3.0_rk, 0.1_rk, 0.0_rk, 100, rs)
   call check(rs%flag == root_ok .and. abs(rs%root - 2.01_rk) <= 0.1_rk + 1.0e-12_rk, &
              'xtol contract: safeguarded midpoint return within xtol of root')

   write (output_unit, '()')
   write (output_unit, '(a,i0,a,i0,a)') 'passed ', nchecks - nfail, ' / ', nchecks, ' checks'
   if (nfail > 0) then
      write (output_unit, '(a,i0,a)') 'FAILURES: ', nfail, ' -- test suite FAILED'
      error stop 1
   end if
   write (output_unit, '(a)') 'all tests passed'

contains

   subroutine check(cond, label)
      logical, intent(in)      :: cond
      character(*), intent(in) :: label
      nchecks = nchecks + 1
      if (cond) then
         write (output_unit, '(a)') 'PASS  ' // label
      else
         nfail = nfail + 1
         write (output_unit, '(a)') 'FAIL  ' // label
      end if
   end subroutine check

   subroutine solve3(eqname, f, a, b, refroot)
      character(*), intent(in) :: eqname
      procedure(func)          :: f
      real(rk), intent(in)     :: a, b, refroot

      call bisection(f, a, b, 1.0e-14_rk, 1.0e-14_rk, 200, rb)
      call check(rb%flag == root_ok, trim(eqname) // '/bisection: flag OK')
      call check(abs(rb%root - refroot) <= 1.0e-12_rk, &
                 trim(eqname) // '/bisection: root within 1e-12 of reference')

      call illinois(f, a, b, 1.0e-14_rk, 1.0e-14_rk, 200, ri)
      call check(ri%flag == root_ok, trim(eqname) // '/illinois: flag OK')
      call check(abs(ri%root - refroot) <= 1.0e-12_rk, &
                 trim(eqname) // '/illinois: root within 1e-12 of reference')

      call safeguarded_secant(f, a, b, 1.0e-14_rk, 1.0e-14_rk, 200, rs)
      call check(rs%flag == root_ok, trim(eqname) // '/safeguarded: flag OK')
      call check(abs(rs%root - refroot) <= 1.0e-12_rk, &
                 trim(eqname) // '/safeguarded: root within 1e-12 of reference')
   end subroutine solve3

   subroutine invalid_case(label, f, a, b)
      character(*), intent(in) :: label
      procedure(func)          :: f
      real(rk), intent(in)     :: a, b

      call bisection(f, a, b, 1.0e-14_rk, 1.0e-14_rk, 100, rb)
      call illinois(f, a, b, 1.0e-14_rk, 1.0e-14_rk, 100, ri)
      call safeguarded_secant(f, a, b, 1.0e-14_rk, 1.0e-14_rk, 100, rs)
      call check(rb%flag == root_err_invalid_bracket .and. &
                 ri%flag == root_err_invalid_bracket .and. &
                 rs%flag == root_err_invalid_bracket, &
                 trim(label) // ': INVALID_BRACKET for all solvers')
   end subroutine invalid_case

end program test_roots
