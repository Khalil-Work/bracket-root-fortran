# Verification — what the tests establish, and the limits

`tests/test_roots.f90` is an assertion suite run by `make test`; it exits
with status 1 (`error stop 1`) if any check fails, so CI gates on it.
Every check prints `PASS`/`FAIL` with a specific message
(96 checks total in the reference run; all passing).

## Reference roots used

| equation | bracket | reference root | source |
|---|---|---|---|
| `x^3 - 2x - 5 = 0` | `[2, 3]` | `2.0945514815423265` | classic (Newton/Wallis) test root |
| `cos x - x = 0` | `[0, 1]` | `0.7390851332151607` | Dottie number |
| `sin x - 1/2 = 0` | `[0, 1.5]` | `pi/6 = 0.5235987755982989...` | **exact analytic root**, `sin(pi/6) = 1/2` |
| `x·e^x - 1 = 0` | `[0, 2]` | `W(1) = 0.5671432904097838...` | Omega constant (Lambert W at 1) |

The `sin x = 1/2` case is the required transcendental with a closed-form
analytic root; `x·e^x = 1` is an additional known-constant case.

## What each test group establishes

1. **Known roots (4 equations x 3 solvers).** Each solver returns
   `flag == OK` and a root within `1e-12` of the reference, at tolerances
   `xtol = ftol = 1e-14`. Establishes that all three methods actually solve
   representative smooth problems to well past single-precision accuracy.

2. **Invalid brackets (4 cases x 3 solvers).** Same-sign endpoints
   (`x^2 + 1`), reversed bracket (`a > b`), non-finite function value at an
   endpoint (`sqrt(x-4) - 1` evaluated at `x = 0` yields NaN), and a NaN
   endpoint. All return `INVALID_BRACKET` (code 1) — the solvers detect bad
   input instead of looping or returning garbage.

3. **Tiny budget.** With `maxiter = 3` and zero tolerances, bisection and the
   safeguarded secant return `MAX_ITERATIONS` (code 2) with
   `iterations == 3` and an iterate still inside the bracket. A zero budget
   returns code 2 with `iterations == 0`. Establishes that the budget flag is
   raised exactly and the reported iterate is always a bracket member.

4. **Machine-precision stagnation.** A synthetic step function
   (`f(x) = -1` for `x < 1`, `+1` otherwise, bracket `[0.5, 1.5]`, unattainable
   residual tolerance) forces bisection and the safeguarded secant down to
   adjacent floating-point numbers, where no representable bracket shrink is
   possible. Both return `NOT_CONVERGED` (code 3). Establishes that the
   stagnation exit fires instead of looping forever.

5. **Bisection contraction law.** With zero tolerances and a 25-iteration
   budget on the cubic, every consecutive pair of recorded history rows
   satisfies `halfwidth(k)/halfwidth(k-1) = 1/2` within `1e-9` (the slack is
   floating-point rounding of the midpoint), and every recorded iterate lies
   inside its recorded bracket. Establishes the per-iteration halving
   guarantee empirically, and that the history rows are self-consistent.

6. **Illinois vs regula falsi on a hard convex case.** `e^x - 2 = 0` on
   `[0, 5]` (root `ln 2`) with tolerance `1e-12`: plain regula falsi
   stagnates one-sidedally and needs **465** iterations; Illinois needs
   **12** (both `flag == OK`, Illinois root within `1e-12` of `ln 2`, and
   regula falsi uses at least 3x the Illinois count). The safeguarded
   secant also converges here. Establishes that the Illinois modification
   removes the false-position stagnation the model doc describes.

7. **Endpoint root.** `f(x) = x` on `[-1, 0]` returns `root = 0`,
   `iterations = 0`, `flag == OK` — the validation short-circuit for an
   endpoint that is exactly a root.

8. **Flag naming.** `flag_message` maps all four codes to their documented
   names (used by the driver's summary CSV).

## Limits of the verification

* **Scalar bracketed problems only.** Every claim concerns one-dimensional
  root finding with a sign-changing bracket. Nothing here applies to
  multivariate systems or unbounded root search.
* **No polynomial rootfinding claims.** The cubic is used as a generic
  smooth function; no companion-matrix, deflation, or all-roots machinery is
  present or verified.
* **Sign-change requirement.** Even-multiplicity roots (no sign change)
  are rejected as invalid brackets by design; tangent roots are out of scope.
* **Platform assumptions.** IEEE-754 double precision with round-to-nearest
  and quiet-NaN propagation (`sqrt` of a negative returns NaN rather than
  trapping). Verified with gfortran 16.2.0 (WinLibs, Windows) and gated on
  ubuntu gfortran in CI; other compilers are untested.
* **Iteration counts are tolerance-dependent.** All quoted counts use the
  tolerances stated with them (`1e-14` driver, `1e-12` convex case); they are
  measurements of this implementation, not universal constants of the
  algorithms.
* **Stagnation test is synthetic.** The `NOT_CONVERGED` path is exercised
  with a discontinuous step function; no claim is made about rates on
  discontinuous problems generally.
