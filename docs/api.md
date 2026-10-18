# Fortran API — module `bracket_root`

Everything lives in `src/roots.f90`, module `bracket_root`, compiled with
`gfortran -std=f2018 -Wall -Wextra -O2 -fall-intrinsics`. All real
computation is double precision (`rk = kind(1.0d0)`), all procedures have
explicit interfaces, and the solvers are `pure` — they perform no I/O and
touch no global state, so they are usable from `pure`/`elemental` contexts
and are trivially thread-safe when called with different arguments.

## Public entities

| entity | kind | description |
|---|---|---|
| `rk` | integer parameter | double-precision kind |
| `root_ok` | integer parameter | flag value 0 |
| `root_err_invalid_bracket` | integer parameter | flag value 1 |
| `root_err_max_iterations` | integer parameter | flag value 2 |
| `root_err_not_converged` | integer parameter | flag value 3 |
| `func` | abstract interface | `pure real(rk) function func(x)` with `real(rk), intent(in) :: x` |
| `step_record` | derived type | one iteration-history row |
| `root_result` | derived type | solver outcome |
| `bisection`, `regula_falsi`, `illinois`, `safeguarded_secant` | pure subroutines | the solvers |
| `flag_message` | pure function | maps a flag code to its documented name |

## Types

```fortran
type :: step_record          ! one row of the iteration history
   integer  :: iteration     ! 1-based iteration index
   real(rk) :: a, b          ! bracket AFTER this iteration's update
   real(rk) :: x             ! point probed this iteration
   real(rk) :: fx            ! f(x) (true residual, always freshly evaluated)
   real(rk) :: halfwidth     ! (b - a)/2 after the update
end type

type :: root_result
   real(rk) :: root          ! best iterate (NaN if the bracket was invalid)
   real(rk) :: residual      ! f(root)
   integer  :: iterations    ! probing iterations used
   integer  :: flag          ! 0 OK, 1 INVALID_BRACKET, 2 MAX_ITERATIONS, 3 NOT_CONVERGED
end type
```

## Solver signatures

All four solvers share one signature; only the proposal logic differs:

```fortran
pure subroutine <solver>(f, a, b, xtol, ftol, maxiter, res, hist)
   procedure(func)              :: f        ! user function
   real(rk), intent(in)         :: a, b     ! initial bracket, a < b
   real(rk), intent(in)         :: xtol     ! stop when (b-a)/2 <= xtol   (>= 0)
   real(rk), intent(in)         :: ftol     ! stop when |f(x)| <= ftol    (>= 0)
   integer,  intent(in)         :: maxiter  ! iteration budget (>= 1)
   type(root_result), intent(out)      :: res
   type(step_record), intent(out), optional :: hist(:)
```

`hist` is optional. If present it should have `size(hist) >= maxiter`; rows
`1..res%iterations` are filled (rows beyond that keep their default-zero
initialisation). Passing `hist` never changes the numerical result.

Solvers: `bisection`, `regula_falsi` (baseline — retained for the stagnation
comparison), `illinois`, `safeguarded_secant`. See `docs/model.md` for the
iteration formula, guarantees, and flag semantics of each.

`flag_message`:

```fortran
pure function flag_message(flag) result(msg)   ! character(len=20)
   integer, intent(in) :: flag                 ! returns OK / INVALID_BRACKET /
                                               ! MAX_ITERATIONS / NOT_CONVERGED
```

## Minimal example

```fortran
program demo
   use bracket_root
   implicit none
   type(root_result) :: res

   call safeguarded_secant(f, 2.0_rk, 3.0_rk, 1.0e-14_rk, 1.0e-14_rk, 100, res)
   print *, res%root, res%residual, res%iterations, trim(flag_message(res%flag))

contains
   pure real(rk) function f(x)
      real(rk), intent(in) :: x
      f = x**3 - 2.0_rk*x - 5.0_rk
   end function f
end program
```

Compile and link:

```sh
gfortran -std=f2018 -Wall -Wextra -O2 -fall-intrinsics -c src/roots.f90
gfortran -std=f2018 -Wall -Wextra -O2 demo.f90 roots.o -o demo
```

(The function must be `pure` because the solvers are `pure`; internal and
module procedures with that shape are accepted as the actual argument.)

## Build targets

| command | effect |
|---|---|
| `make` | compile `src/roots.f90` and the test binary |
| `make test` | build and run `tests/test_roots.f90` (nonzero exit on failure) |
| `make run` | build and run `examples/compare.f90` → `build/csv/*.csv` |
| `make results` | `run` + `tools/make_results.py` → figures, GIF, `docs/results/run_summary.json` |
