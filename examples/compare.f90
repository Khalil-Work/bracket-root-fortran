!> compare.f90 — Bracket Root method-comparison driver.
!>
!> Solves a documented set of scalar test equations with the three bracketed
!> solvers (bisection, Illinois, safeguarded secant) and writes, per equation
!> and method, a CSV iteration history
!>
!>    iteration, a, b, x, fx, halfwidth
!>
!> into build/csv/<equation>_<method>.csv, plus a summary table in
!> build/csv/summary.csv.  Run from the repository root via `make run`.
module equations
   use bracket_root, only: rk, func
   implicit none
contains

   !> Wallis/Newton classic: x^3 - 2x - 5 = 0, root 2.0945514815423265...
   pure real(rk) function f_cubic(x)
      real(rk), intent(in) :: x
      f_cubic = x**3 - 2.0_rk*x - 5.0_rk
   end function f_cubic

   !> Dottie number: cos x = x, root 0.7390851332151607...
   pure real(rk) function f_cosx(x)
      real(rk), intent(in) :: x
      f_cosx = cos(x) - x
   end function f_cosx

   !> Smooth transcendental with an exactly known analytic root:
   !> sin x = 1/2  =>  x* = pi/6 = 0.5235987755982989...
   pure real(rk) function f_sinx(x)
      real(rk), intent(in) :: x
      f_sinx = sin(x) - 0.5_rk
   end function f_sinx

   !> Omega constant: x e^x = 1  =>  x* = W(1) = 0.56714329040978387...
   pure real(rk) function f_omega(x)
      real(rk), intent(in) :: x
      f_omega = x*exp(x) - 1.0_rk
   end function f_omega

end module equations

program compare
   use, intrinsic :: iso_fortran_env, only: output_unit
   use bracket_root
   use equations
   implicit none

   integer, parameter :: maxit = 200
   real(rk), parameter :: xtol = 1.0e-14_rk
   real(rk), parameter :: ftol = 1.0e-14_rk

   integer :: usum
   type(step_record), allocatable :: hist(:)

   allocate (hist(maxit))

   open (newunit=usum, file='build/csv/summary.csv', status='replace', action='write')
   write (usum, '(a)') 'method,equation,iterations,final_residual,flag,flag_name,' // &
                       'root,reference_root,abs_error'

   write (output_unit, '(a)') '=== Bracket Root: method comparison driver ==='
   write (output_unit, '(a,es9.2,a,es9.2,a,i0)') 'tolerances: xtol = ', xtol, &
      '  ftol = ', ftol, '  maxiter = ', maxit

   call run_case('cubic', 'x^3 - 2x - 5 = 0', f_cubic, 2.0_rk, 3.0_rk, 2.0945514815423265_rk)
   call run_case('cosx', 'cos(x) - x = 0', f_cosx, 0.0_rk, 1.0_rk, 0.7390851332151607_rk)
   call run_case('sinx', 'sin(x) - 1/2 = 0 (root pi/6)', f_sinx, 0.0_rk, 1.5_rk, 0.5235987755982989_rk)
   call run_case('omega', 'x*exp(x) - 1 = 0 (root W(1))', f_omega, 0.0_rk, 2.0_rk, 0.56714329040978387_rk)

   close (usum)
   write (output_unit, '(a)') 'wrote build/csv/<equation>_<method>.csv and build/csv/summary.csv'

contains

   subroutine run_case(name, formula, f, a, b, refroot)
      character(*), intent(in) :: name, formula
      procedure(func)          :: f
      real(rk), intent(in)     :: a, b, refroot
      type(root_result)        :: res

      write (output_unit, '()')
      write (output_unit, '(a,a,a,f7.4,a,f7.4,a,a)') 'Equation ', trim(name), &
         ': [', a, ',', b, ']  ', trim(formula)

      call bisection(f, a, b, xtol, ftol, maxit, res, hist)
      call report('bisection', name, res, refroot)

      call illinois(f, a, b, xtol, ftol, maxit, res, hist)
      call report('illinois', name, res, refroot)

      call safeguarded_secant(f, a, b, xtol, ftol, maxit, res, hist)
      call report('safeguarded_secant', name, res, refroot)
   end subroutine run_case

   subroutine report(method, eqname, res, refroot)
      character(*), intent(in)    :: method, eqname
      type(root_result), intent(in) :: res
      real(rk), intent(in)        :: refroot

      character(len=80) :: fname
      integer :: u, k

      fname = 'build/csv/' // trim(eqname) // '_' // trim(method) // '.csv'
      open (newunit=u, file=fname, status='replace', action='write')
      write (u, '(a)') 'iteration,a,b,x,fx,halfwidth'
      do k = 1, res%iterations
         write (u, '(i0,a,es24.16,a,es24.16,a,es24.16,a,es24.16,a,es24.16)') &
            hist(k)%iteration, ',', hist(k)%a, ',', hist(k)%b, ',', &
            hist(k)%x, ',', hist(k)%fx, ',', hist(k)%halfwidth
      end do
      close (u)

      write (usum, '(a,",",a,",",i0,",",es24.16,",",i0,",",a,",",es24.16,",",es24.16,",",es24.16)') &
         trim(method), trim(eqname), res%iterations, res%residual, res%flag, &
         trim(flag_message(res%flag)), res%root, refroot, abs(res%root - refroot)

      write (output_unit, '(a22,i8,es13.3,2x,a,a,es10.2)') method, res%iterations, &
         res%residual, trim(flag_message(res%flag)), ' |x-r| = ', abs(res%root - refroot)
   end subroutine report

end program compare
