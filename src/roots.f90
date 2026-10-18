!> bracket_root — native Fortran bracketed scalar root solvers.
!>
!> Written from scratch for the Bracket Root project.  No third-party solver
!> code (in particular no published Brent implementation) was consulted or
!> copied; the safeguard logic in `safeguarded_secant` is an original design
!> described in docs/model.md.
!>
!> Provided solvers
!>   bisection           — guaranteed half-width contraction per iteration
!>   regula_falsi        — classic false position (baseline, stagnates on
!>                         one-sided convex/concave problems)
!>   illinois            — false position with the Illinois modification
!>                         (retained-endpoint value halving)
!>   safeguarded_secant  — secant steps constrained to the bracket with a
!>                         forced bisection backstop (original logic)
!>
!> All solvers share one verified kernel: validate the bracket, iterate,
!> record history, and report (root, residual, iterations, error flag).
module bracket_root
   use, intrinsic :: ieee_arithmetic, only: &
      ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   private

   !> Working precision: IEEE double precision everywhere.
   integer, parameter, public :: rk = kind(1.0d0)

   !> Error-flag codes (see docs/model.md for full semantics).
   integer, parameter, public :: root_ok                  = 0
   integer, parameter, public :: root_err_invalid_bracket = 1
   integer, parameter, public :: root_err_max_iterations  = 2
   integer, parameter, public :: root_err_not_converged   = 3

   !> Private method selectors for the shared kernel.
   integer, parameter :: method_bisect    = 1
   integer, parameter :: method_regula    = 2
   integer, parameter :: method_illinois  = 3
   integer, parameter :: method_safeguard = 4

   !> One row of the per-iteration history.  The bracket [a, b] is the state
   !> AFTER the iteration's update; x is the point probed during it.
   type, public :: step_record
      integer  :: iteration = 0
      real(rk) :: a         = 0.0_rk
      real(rk) :: b         = 0.0_rk
      real(rk) :: x         = 0.0_rk
      real(rk) :: fx        = 0.0_rk
      real(rk) :: halfwidth = 0.0_rk
   end type step_record

   !> Solver outcome: best iterate, its residual, iterations used, flag.
   type, public :: root_result
      real(rk) :: root       = 0.0_rk
      real(rk) :: residual   = 0.0_rk
      integer  :: iterations = 0
      integer  :: flag       = root_ok
   end type root_result

   !> Interface for user-supplied scalar functions.
   abstract interface
      pure real(rk) function func(x)
         import :: rk
         real(rk), intent(in) :: x
      end function func
   end interface

   public :: bisection, regula_falsi, illinois, safeguarded_secant
   public :: flag_message, func

contains

   !> Bisection: x_{k+1} = a_k + (b_k - a_k)/2, bracket halved every step.
   pure subroutine bisection(f, a, b, xtol, ftol, maxiter, res, hist)
      procedure(func)                    :: f
      real(rk), intent(in)               :: a, b, xtol, ftol
      integer, intent(in)                :: maxiter
      type(root_result), intent(out)     :: res
      type(step_record), intent(out), optional :: hist(:)

      call solve_kernel(f, a, b, xtol, ftol, maxiter, method_bisect, res, hist)
   end subroutine bisection

   !> Classic regula falsi (false position).  Retained as a baseline to
   ! demonstrate one-sided stagnation; use `illinois` in practice.
   pure subroutine regula_falsi(f, a, b, xtol, ftol, maxiter, res, hist)
      procedure(func)                    :: f
      real(rk), intent(in)               :: a, b, xtol, ftol
      integer, intent(in)                :: maxiter
      type(root_result), intent(out)     :: res
      type(step_record), intent(out), optional :: hist(:)

      call solve_kernel(f, a, b, xtol, ftol, maxiter, method_regula, res, hist)
   end subroutine regula_falsi

   !> Regula falsi with the Illinois modification: when the same endpoint is
   ! retained on two (or more) consecutive iterations, its stored function
   ! value is halved; alternating phases keep true values, giving secant-like
   ! superlinear behaviour (order about 1.44).
   pure subroutine illinois(f, a, b, xtol, ftol, maxiter, res, hist)
      procedure(func)                    :: f
      real(rk), intent(in)               :: a, b, xtol, ftol
      integer, intent(in)                :: maxiter
      type(root_result), intent(out)     :: res
      type(step_record), intent(out), optional :: hist(:)

      call solve_kernel(f, a, b, xtol, ftol, maxiter, method_illinois, res, hist)
   end subroutine illinois

   !> Safeguarded secant: secant proposals from the two latest iterates,
   ! accepted only when they land strictly inside the current bracket; a
   ! forced midpoint step follows any two consecutive non-halving iterations,
   ! which caps the worst case at a small multiple of the bisection count.
   pure subroutine safeguarded_secant(f, a, b, xtol, ftol, maxiter, res, hist)
      procedure(func)                    :: f
      real(rk), intent(in)               :: a, b, xtol, ftol
      integer, intent(in)                :: maxiter
      type(root_result), intent(out)     :: res
      type(step_record), intent(out), optional :: hist(:)

      call solve_kernel(f, a, b, xtol, ftol, maxiter, method_safeguard, res, hist)
   end subroutine safeguarded_secant

   !> Human-readable name for an error-flag code.
   pure function flag_message(flag) result(msg)
      integer, intent(in)   :: flag
      character(len=20)     :: msg

      select case (flag)
      case (root_ok)
         msg = 'OK'
      case (root_err_invalid_bracket)
         msg = 'INVALID_BRACKET'
      case (root_err_max_iterations)
         msg = 'MAX_ITERATIONS'
      case (root_err_not_converged)
         msg = 'NOT_CONVERGED'
      case default
         msg = 'UNKNOWN'
      end select
   end function flag_message

   !> Shared iteration kernel for all four methods.
   pure subroutine solve_kernel(f, a, b, xtol, ftol, maxiter, method, res, hist)
      procedure(func)                :: f
      real(rk), intent(in)           :: a, b, xtol, ftol
      integer, intent(in)            :: maxiter, method
      type(root_result), intent(out) :: res
      type(step_record), intent(out), optional :: hist(:)

      real(rk) :: aa, bb, faw, fbw   ! working bracket + stored endpoint values
      real(rk) :: xc, fc, xp, fp     ! secant pair: current and previous iterate
      real(rk) :: xn, fn             ! this iteration's point and residual
      real(rk) :: xm, fm             ! midpoint returned on xtol convergence
      real(rk) :: s, m, denom        ! secant proposal, midpoint, denominator
      real(rk) :: wprev, wnew        ! bracket widths before/after the update
      integer  :: it, nslow          ! iteration index, non-halving streak
      integer  :: last_side          ! 1 = a replaced last, 2 = b replaced last
      logical  :: accept, done

      res%root       = ieee_value(1.0_rk, ieee_quiet_nan)
      res%residual   = ieee_value(1.0_rk, ieee_quiet_nan)
      res%iterations = 0
      res%flag       = root_ok

      ! ---- bracket validation ------------------------------------------
      if (.not. (ieee_is_finite(a) .and. ieee_is_finite(b))) then
         res%flag = root_err_invalid_bracket
         return
      end if
      if (a >= b) then
         res%flag = root_err_invalid_bracket
         return
      end if
      faw = f(a)
      fbw = f(b)
      if (.not. (ieee_is_finite(faw) .and. ieee_is_finite(fbw))) then
         res%flag = root_err_invalid_bracket
         return
      end if
      if (abs(faw) <= 0.0_rk) then        ! endpoint already a root
         res%root = a
         res%residual = faw
         res%flag = root_ok
         return
      end if
      if (abs(fbw) <= 0.0_rk) then
         res%root = b
         res%residual = fbw
         res%flag = root_ok
         return
      end if
      if (sign(1.0_rk, faw)*sign(1.0_rk, fbw) > 0.0_rk) then   ! same sign
         res%flag = root_err_invalid_bracket
         return
      end if

      ! ---- method-specific initialisation -------------------------------
      aa = a
      bb = b
      xc = 0.0_rk
      fc = 0.0_rk
      xp = 0.0_rk
      fp = 0.0_rk
      nslow = 0
      last_side = 0
      if (method == method_safeguard) then
         xc = a                           ! first secant pair spans the bracket,
         fc = faw                         ! so the first proposal is the
         xp = b                           ! false-position point
         fp = fbw
      end if

      ! ---- iteration loop ------------------------------------------------
      done = .false.
      do it = 1, maxiter
         wprev = bb - aa

         ! propose the next point according to the method
         select case (method)
         case (method_bisect)
            xn = aa + 0.5_rk*(bb - aa)

         case (method_regula, method_illinois)
            denom = faw - fbw
            if (abs(denom) > 0.0_rk) then
               xn = (faw*bb - fbw*aa)/denom
            else
               xn = aa + 0.5_rk*(bb - aa)
            end if

         case (method_safeguard)
            m = aa + 0.5_rk*(bb - aa)
            accept = .false.
            s = ieee_value(1.0_rk, ieee_quiet_nan)
            if (abs(xc - xp) > 0.0_rk .and. abs(fc - fp) > 0.0_rk) then
               s = xc - fc*(xc - xp)/(fc - fp)
            end if
            if (ieee_is_finite(s)) then
               ! S1: the secant step must land strictly inside the current
               ! bracket, so every accepted step shrinks the bracket.
               if (s > aa .and. s < bb) then
                  accept = .true.
               end if
            end if
            ! S2: after two consecutive iterations that failed to halve the
            ! bracket, fall back to the midpoint (always halves).
            if (nslow >= 2) accept = .false.
            if (accept) then
               xn = s
            else
               xn = m
            end if
         end select

         fn = f(xn)
         if (.not. ieee_is_finite(fn)) then  ! non-finite function behaviour
            res%root = xn
            res%residual = fn
            res%iterations = it
            res%flag = root_err_invalid_bracket
            done = .true.
            exit
         end if

         ! update the bracket by sign; Illinois halves the stored value of
         ! an endpoint that is retained for a second consecutive iteration
         if (sign(1.0_rk, fn)*sign(1.0_rk, faw) > 0.0_rk) then
            if (method == method_illinois .and. last_side == 1) fbw = 0.5_rk*fbw
            aa = xn
            faw = fn
            last_side = 1
         else
            if (method == method_illinois .and. last_side == 2) faw = 0.5_rk*faw
            bb = xn
            fbw = fn
            last_side = 2
         end if

         ! shift the secant memory
         if (method == method_safeguard) then
            xp = xc
            fp = fc
            xc = xn
            fc = fn
         end if

         wnew = bb - aa
         res%root = xn
         res%residual = fn
         res%iterations = it
         if (present(hist)) then
            if (it <= size(hist)) then
               hist(it) = step_record(it, aa, bb, xn, fn, 0.5_rk*wnew)
            end if
         end if

         ! convergence, residual test (covers the exact-zero case for
         ! ftol >= 0): the probed point itself is returned
         if (abs(fn) <= ftol) then
            res%flag = root_ok
            done = .true.
            exit
         end if

         ! convergence, bracket-width test: return the midpoint of the
         ! updated bracket.  The probed point xn is an ENDPOINT of [aa, bb],
         ! so it carries no half-width error bound; the midpoint m satisfies
         ! |m - r| <= (bb - aa)/2 <= xtol for every r in the bracket, up to
         ! the rounding of the midpoint sum.  This performs one additional
         ! function evaluation that is not counted in res%iterations (which
         ! counts probing iterations), and the history keeps recording the
         ! probed points, so on this exit res%root need not equal the last
         ! history row's x.
         if (0.5_rk*wnew <= xtol) then
            xm = aa + 0.5_rk*(bb - aa)
            fm = f(xm)
            if (.not. ieee_is_finite(fm)) then
               res%root = xm
               res%residual = fm
               res%flag = root_err_invalid_bracket
               done = .true.
               exit
            end if
            res%root = xm
            res%residual = fm
            res%flag = root_ok
            done = .true.
            exit
         end if

         ! update the safeguard streak counter
         if (method == method_safeguard) then
            if (wnew <= 0.5_rk*wprev) then
               nslow = 0
            else
               nslow = nslow + 1
            end if
         end if

         ! stagnation guard: every valid iteration must shrink the bracket
         if (wnew >= wprev) then
            res%flag = root_err_not_converged
            done = .true.
            exit
         end if
      end do

      if (.not. done) then                  ! budget exhausted unconverged
         res%flag = root_err_max_iterations
         res%iterations = max(maxiter, 0)
      end if
   end subroutine solve_kernel

end module bracket_root
