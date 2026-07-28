# polyfrac

[![install](https://img.shields.io/badge/install-from%20GitHub-blue)](https://github.com/nickharris808/polyfrac#install)
[![CI](https://img.shields.io/badge/ci-passing-brightgreen)](https://github.com/nickharris808/polyfrac/actions/workflows/ci.yml)
[![tests](https://img.shields.io/badge/tests-66%20passing-brightgreen)](tests/)
[![python](https://img.shields.io/badge/python-3.9%2B-blue)](pyproject.toml)
[![license](https://img.shields.io/badge/license-MIT-green)](LICENSE)
![deps](https://img.shields.io/badge/dependencies-none-brightgreen)

**Exact polynomial and rational-function arithmetic over ℚ, with Sturm-sequence real-root counting.
No floating point anywhere.**

## Why this exists

Testing a parameter at chosen values cannot establish that a property holds *everywhere* in a range —
a crossing may sit between two samples, and no sampling density fixes that. Interval arithmetic gives
you a conservative bound; a numerical root-finder gives you an approximate answer with an error term.
Neither gives you an **integer**.

Sturm's theorem does: the exact count of distinct real roots in an interval, computed by counting
sign changes over rational arithmetic. Zero roots plus a known sign at one endpoint settles the sign
for the whole continuum. That is the difference between "we tested 10,000 points" and "there is no
crossing", and it is what you need if the answer goes into a specification.

Sampling a parameter at chosen points cannot tell you that a property holds *everywhere* in a range —
a crossing may sit between two samples. `polyfrac` counts the real roots in an interval **exactly**,
as an integer, using Sturm's theorem over rational arithmetic. If the count is zero and the sign at one
endpoint is known, the sign is settled for the whole interval.

Zero dependencies. Pure standard library (`fractions`). ~390 lines you can read in one sitting.

## Install

```
# from GitHub (PyPI release pending)
pip install "polyfrac @ git+https://github.com/nickharris808/polyfrac.git"
```

> `pip install polyfrac` does not work yet — the package is not on PyPI. Install from GitHub as
> shown above. The distribution builds and is `twine check`-clean, with no unpublished
> dependencies, so it is ready to upload whenever that happens.

## 30-second quickstart

```python
from fractions import Fraction
from polyfrac import Poly, count_roots, positive_on

p = Poly.from_roots([1, 2, 3])         # (x-1)(x-2)(x-3)

count_roots(p, 0, 4)                   # 3   — exactly three real roots in (0, 4]
count_roots(p, 0, Fraction(3, 2))      # 1   — exactly one, namely x = 1
count_roots(p, 4, 10)                  # 0   — none out here

positive_on(Poly([-1, 1]), 2, 5)       # certificate that x - 1 > 0 on (2, 5]
# {'endpoint_a': '2', 'endpoint_b': '5', 'P_a': '1', 'P_b': '4',
#  'n_roots_in_interval': 0, 'positive_on_interval': True}
```

Intervals are half-open, `(a, b]` — left-exclusive, right-inclusive — which is the convention Sturm's
theorem states most cleanly. Repeated roots count **once**: `count_roots` counts *distinct* roots.

## Worked example — certifying a bound over a whole interval

Suppose a guarded procedure has failure probability `p/500` and an unguarded baseline has `3p/5`,
where `p ∈ (0, 1]` is an environment parameter you do not control. Does the guarded procedure meet a
target of `1/20` for **every** `p`, or only for the ones you happened to test?

```python
from fractions import Fraction
from polyfrac import Poly, count_roots, positive_on

target   = Fraction(1, 20)
guarded  = Poly([0, Fraction(1, 500)])     # p/500
baseline = Poly([0, Fraction(3, 5)])       # 3p/5

# Does target - guarded(p) stay positive on the whole interval?
positive_on(Poly([target]) - guarded, 0, 1)["positive_on_interval"]
# True  -> the guarded procedure meets the target for EVERY p in (0, 1]

# The baseline does not. Where exactly does it cross?
count_roots(Poly([target]) - baseline, 0, 1)      # 1  -> exactly one crossing
(Poly([target]) - baseline).eval(Fraction(1, 12)) # 0  -> the crossing is exactly p = 1/12
```

That last line is the point: the crossing is reported as the **exact rational `1/12`**, not as
`0.0833333...`. You can put it in a specification.

## Rational functions and linear systems

`PolyFrac` is a ratio of two `Poly`, kept in lowest terms, and `gauss_solve` solves linear systems
whose entries are rational functions of the parameter. This is what you need to solve a parameterised
Markov chain symbolically — the hitting probability comes out as an exact `N(p)/D(p)` rather than as a
number for one chosen `p`.

```python
from polyfrac import Poly, PolyFrac, gauss_solve

one = PolyFrac(Poly([1]), Poly([1]))
P   = PolyFrac(Poly([0, 1]), Poly([1]))    # the parameter p
(h,) = gauss_solve([[one]], [P])           # solve 1*h = p
# h is exactly p, as a rational function
```

## API

| Name | What it does |
|---|---|
| `Poly(coeffs)` | Dense univariate polynomial, `Fraction` coefficients, ascending order |
| `Poly.from_roots(rs)` | The monic polynomial with exactly those roots |
| `Poly.const`, `Poly.x` | Constant and identity constructors |
| `+ - * neg scale divmod gcd derivative eval degree is_zero` | Exact polynomial arithmetic |
| `PolyFrac(num, den)` | Rational function, auto-reduced |
| `gauss_solve(A, b)` | Gaussian elimination over the field of rational functions |
| `sturm_chain(P)` | The squarefree Sturm chain |
| `count_roots(P, a, b)` | Exact integer count of **distinct** real roots in `(a, b]` |
| `positive_on(P, a, b)` | Certificate dict: endpoints, values, root count, verdict |

## Why not SymPy?

SymPy will do all of this and much more. `polyfrac` exists for the case where you want the root count
and nothing else, with no dependency, no import cost, and a source file short enough to audit before
you trust it in a certification path.

## Honest scope

**Exactness is enforced at the boundary, not merely intended.** Passing a `float` raises
`InexactInput` rather than being converted, because `0.1` is not one tenth — it is
`3602879701896397/36028797018963968`, and exact arithmetic on it produces an exact answer to a
question you did not ask. Use `Fraction(1, 10)`, the string `"0.1"`, or the explicit
`Poly.from_floats()` when a float really is the value you mean.

```python
>>> Poly([0.1])
InexactInput: refusing the float 0.1: it is a binary approximation ...
>>> Poly(["0.1"]).c[0]
Fraction(1, 10)
```

**What it proves.** An exact count of *distinct* real roots in the half-open interval `(a, b]`, and
an exact sign certificate over an interval. No sampling, no tolerance, no numerical root finding —
the count comes from sign changes in a Sturm chain over exact rationals.

**What it does not prove.**

- Univariate only, over ℚ. For multivariate sign conditions you want cylindrical algebraic
  decomposition, which this package does not implement.
- Multiplicities are not reported. `Poly.from_roots([2, 2, 2])` has *one* distinct root.
- Irrational and complex roots are outside the interval arithmetic entirely; only real roots at
  rational-bounded intervals are counted.
- The interval convention is half-open `(a, b]` — the left endpoint is excluded and the right is
  included. Off-by-one here is silent, so it is asserted in the test suite.
- `gauss_solve` requires a non-singular system; a singular one raises rather than returning a value.
- The zero polynomial raises. It is zero at every point, so no finite root count exists, and
  returning `0` would be a confident wrong answer.

## Troubleshooting

**`InexactInput: refusing the float 0.1`.** Deliberate. `0.1` is not one tenth in binary, so exact
arithmetic on it answers a question you did not ask. Use `Fraction(1, 10)`, the string `"0.1"`, or
`Poly.from_floats([...])` if the binary value really is what you mean. This applies to interval
bounds too — a float `a` or `b` silently shifts what is being certified.

**`ValueError: need a < b`.** The interval is empty or reversed. Bounds are not sorted for you,
because swapping them would answer a different question than the one asked.

**`ValueError: the zero polynomial is zero at every point`.** It has infinitely many roots in any
interval, so there is no finite count to return. `0` would be a confident wrong answer.

**`count_roots` returned fewer roots than I expected.** It counts *distinct* roots.
`Poly.from_roots([2, 2, 2])` has one. Multiplicities are not reported.

**A root exactly at an endpoint is or is not counted.** The interval is half-open `(a, b]` — left
excluded, right included. `count_roots(Poly.from_roots([2]), 2, 5)` is `0`;
`count_roots(Poly.from_roots([2]), 0, 2)` is `1`.

**`positive_on` says `False` but the polynomial looks positive.** It certifies over the whole
interval. One root inside is enough to refuse, even if both endpoints are positive. Check
`n_roots_in_interval` in the returned dict.

**`gauss_solve` raised instead of returning.** The system is singular over the field of rational
functions, so there is no unique solution to return.

## Performance

Measured, not estimated: `count_roots` on a degree-40 polynomial with 40 distinct roots takes about
**8.5 ms**; degree 20 about **2.1 ms**; degree 10 about **0.6 ms**. Cost is dominated by the Sturm
chain, which is quadratic in the degree, and by `Fraction` growth in the coefficients. Nothing here
has needed optimising — if you hit a case that is slow, it is worth reporting rather than working
around.

## Tests

```
pip install -e ".[test]" && pytest
```

66 tests, including the worked example above and explicit checks that the arithmetic is exact
(`(1/3) * 3 == 1`, not `0.9999999999999999`).

## Where this came from

`polyfrac` was extracted from a formal-methods system that certifies reliability orderings over a
continuum of environment parameters — the worked example above is a stripped-down version of what it
does for real. The arithmetic is here under MIT. The certification pipeline built on it, and the
models it runs against, are the commercial offering.

## The portfolio

Five small, independently useful tools built around one idea: **a verdict you cannot check is not a verdict.**

| | |
|---|---|
| [`minicheck`](https://github.com/nickharris808/minicheck) | An explicit-state model checker in ~390 lines, with a CLI. Shortest counterexamples, no required dependencies. |
| [`protocol-bench`](https://github.com/nickharris808/protocol-bench) | 15 published IEEE 802.11 / 3GPP procedures with ground truth. A claimed detection must **replay**. |
| [`minicheck-mcp`](https://github.com/nickharris808/minicheck-mcp) | The checker as an **MCP server** — let an agent verify a state machine instead of guessing. |
| [`polyfrac`](https://github.com/nickharris808/polyfrac) ← *you are here* | Exact polynomial + rational-function arithmetic over ℚ with Sturm real-root counting. Zero deps. |
| [`failclosed`](https://github.com/nickharris808/failclosed) | Default-deny ASGI middleware: a gated endpoint succeeds only on an affirmative verdict. |
| [`protocol-bench-action`](https://github.com/nickharris808/protocol-bench-action) | Score a submission in CI and fail the build if a claimed detection cannot be proved |

Try it in your browser: **[live demo](https://huggingface.co/spaces/nickh007/protocol-bench-demo)** · Ground-truth tasks: **[dataset](https://huggingface.co/datasets/nickh007/protocol-bench)**

### The commercial offering

These are the engine. What is **not** open source is what makes it useful at scale: the maintained
hazard-property corpora, composition analysis that finds hazards existing only when two components
are combined, the trust-model sensitivity sweep, and the evidence trail that makes a verdict auditable
after the fact. The tools above are MIT and stay that way.

## Licence

MIT. See `LICENSE`.
