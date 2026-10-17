# Model — bracketed scalar root solvers

All methods solve the scalar problem

> given a continuous `f : R -> R` and a bracket `[a, b]` with `f(a)·f(b) < 0`,
> find `r` in `(a, b)` with `f(r) = 0`.

Every solver maintains a current bracket `[a_k, b_k]` that is guaranteed to
contain a sign change, produces one new point `x_k` per iteration, and replaces
one endpoint with it. All arithmetic is IEEE double precision.

## Symbols

| symbol | meaning |
|---|---|
| `f` | user-supplied pure scalar function (double precision) |
| `a, b` | initial bracket endpoints, `a < b`, `f(a)·f(b) <= 0` |
| `a_k, b_k` | bracket after iteration `k` |
| `x_k` | point probed during iteration `k` (becomes an endpoint) |
| `w_k` | bracket width `b_k - a_k`; half-width `(b_k - a_k)/2` is the bracket error bound |
| `f_k` | `f(x_k)` |
| `v_a, v_b` | *stored* endpoint values used by regula falsi / Illinois (may be modified; see below) |
| `xtol` | bracket stopping tolerance: stop when `w_k/2 <= xtol` |
| `ftol` | residual stopping tolerance: stop when `|f_k| <= ftol` |
| `maxiter` | iteration budget (`>= 1` expected; tolerances `>= 0` expected) |

## Methods

### 1. Bisection

Iteration formula:

```
x_k = a_{k-1} + (b_{k-1} - a_{k-1}) / 2
```

Bracket update: replace `a` if `f(x_k)` has the same sign as `f(a)`,
otherwise replace `b`.

* **Guarantee:** `w_k = w_0 / 2^k` up to floating-point rounding, so
  `|x_k - r| <= w_0 / 2^k`. Linear convergence, rate exactly 1/2.
* **Caveats:** slow (about 3.3 iterations per decimal digit); needs nothing
  from `f` beyond continuity and a sign change.
* The test suite verifies the `1/2` contraction ratio directly on the
  recorded histories (within `1e-9`, the rounding slack).

### 2. Regula falsi (false position) — baseline

Iteration formula (zero of the secant line through the bracket endpoints):

```
x_k = (v_a·b_{k-1} - v_b·a_{k-1}) / (v_a - v_b)
```

with stored values `v_a = f(a)`, `v_b = f(b)` (true values in this method).

* **Behaviour:** superlinear *while the iterates alternate sides* (the
  proposal is then a secant step through two fresh points), but on convex or
  concave `f` one endpoint is retained forever; the bracket width stagnates
  and the convergence degrades to linear with a problem-dependent factor.
* This method is provided as a comparison baseline only; measured below at
  465 iterations on `e^x - 2 = 0` over `[0, 5]` at tolerance `1e-12`.

### 3. Illinois (regula falsi + Illinois modification)

Same proposal formula as regula falsi, but the *stored* endpoint values are
algorithmic state:

> when the same endpoint is retained on two (or more) **consecutive**
> iterations, its stored value is halved:
> replacing the `a`-side for a second time sets `v_b <- v_b / 2`
> (symmetrically `v_a <- v_a / 2` for the `b`-side).

* **Why it works:** halving the retained value doubles its weight in the
  proposal, pulling the next point across the root and breaking the one-sided
  crawl. In alternating phases no halving occurs, the stored values stay
  exact, and each proposal is a secant step through two fresh points.
* **Guarantees/caveats:** no a-priori bound better than bisection's; in
  practice superlinear with asymptotic order about 1.44. Stored endpoint
  values may differ from the true residuals `f(a_k)`, `f(b_k)` — only the
  reported iterate residuals are exact. Measured 12 iterations on the hard
  convex case above (vs 465 for plain regula falsi).

### 4. Safeguarded secant (original design)

Proposal from the secant line through the **two most recent iterates**
`x_k, x_{k-1}` (regardless of which endpoints they became):

```
s_k = x_k - f_k · (x_k - x_{k-1}) / (f_k - f_{k-1})
```

Safeguard rules (our own; no published Brent implementation was consulted
or copied):

* **S1 — containment:** `s_k` is accepted only if it is finite and lies
  *strictly inside* the current bracket `[a_k, b_k]`. Every accepted step
  replaces an endpoint with an interior point, so the bracket shrinks.
* **S2 — forced bisection:** if the bracket failed to at least halve in two
  consecutive iterations, the next point is the midpoint
  `a_k + (b_k - a_k)/2`, which always halves the bracket.

Consequences:

* **Worst-case guarantee:** the bracket halves at least every third
  iteration, so width-driven stopping needs at most
  `3·ceil(log2(w_0 / (2·xtol))) + O(1)` iterations — a small constant
  multiple of the bisection count, independent of how badly the secant
  proposals behave.
* **Typical behaviour:** superlinear (secant-grade) whenever the secant
  proposals stay inside the bracket; measured 6–20 iterations on the test
  equations vs 44–47 for bisection.
* The first proposal is the false-position point of the initial bracket
  (the initial secant pair is `(a, f(a))`, `(b, f(b))`).

## Stopping criteria and error flags

The kernel stops with `flag` set as follows:

| code | name | raised when |
|---|---|---|
| 0 | `OK` | converged: `|f(x_k)| <= ftol` (covers the exact-zero case) or `w_k/2 <= xtol`; an endpoint that is exactly a root is returned immediately with 0 iterations |
| 1 | `INVALID_BRACKET` | non-finite `a` or `b`; `a >= b`; non-finite `f(a)` or `f(b)`; `f(a)·f(b) > 0` (same nonzero sign); or a non-finite `f` value encountered during iteration |
| 2 | `MAX_ITERATIONS` | the iteration budget was exhausted before either stopping test held |
| 3 | `NOT_CONVERGED` | machine-precision stagnation: an iteration failed to strictly shrink the bracket (no representable progress) while the tolerances were still unmet |

On any non-`OK` exit the result still reports the last iterate, its true
residual, and the number of iterations used, so callers can inspect how far
the solver got. Tolerances are expected to be non-negative and `maxiter >= 1`;
violating budgets simply returns `MAX_ITERATIONS`.

## What is *not* claimed

* Only **scalar, bracketed** problems with a sign change are solved. Roots of
  even multiplicity (no sign change) are outside scope.
* No polynomial-specific theory (no companion matrices, no root counts).
* The history records the bracket *after* each iteration's update; iteration
  counts are the number of function-probing iterations.
