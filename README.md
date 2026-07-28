# polyfrac

[![install](https://img.shields.io/badge/install-from%20GitHub-blue)](https://github.com/nickharris808/polyfrac#install)
[![CI](https://img.shields.io/badge/ci-passing-brightgreen)](https://github.com/nickharris808/polyfrac/actions/workflows/ci.yml)
[![tests](https://img.shields.io/badge/tests-18%20passing-brightgreen)](tests/)
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

Zero dependencies. Pure standard library (`fractions`). ~250 lines you can read in one sitting.

## Install

```
# from GitHub (PyPI release pending)
pip install "polyfrac @ git+https://github.com/nickharris808/polyfrac.git"
```

> `pip install polyfrac` will work once the PyPI release lands. The distribution is built and `twine check`-clean; publication is pending.

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

## Scope

Univariate, over ℚ. `count_roots` counts distinct real roots; multiplicities are not reported.
`gauss_solve` assumes the system is non-singular over the field of rational functions. For multivariate
sign conditions you want cylindrical algebraic decomposition, which this package does not implement.

## Tests

```
pip install -e ".[test]" && pytest
```

18 tests, including the worked example above and explicit checks that the arithmetic is exact
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
| [`minicheck`](https://github.com/nickharris808/minicheck) | An explicit-state model checker in ~560 lines. Shortest counterexamples, no required dependencies. |
| [`protocol-bench`](https://github.com/nickharris808/protocol-bench) | 15 published IEEE 802.11 / 3GPP procedures with ground truth. A claimed detection must **replay**. |
| [`minicheck-mcp`](https://github.com/nickharris808/minicheck-mcp) | The checker as an **MCP server** — let an agent verify a state machine instead of guessing. |
| [`polyfrac`](https://github.com/nickharris808/polyfrac) ← *you are here* | Exact polynomial + rational-function arithmetic over ℚ with Sturm real-root counting. Zero deps. |
| [`failclosed`](https://github.com/nickharris808/failclosed) | Default-deny ASGI middleware: a gated endpoint succeeds only on an affirmative verdict. |

Try it in your browser: **[live demo](https://huggingface.co/spaces/nickh007/protocol-bench-demo)** · Ground-truth tasks: **[dataset](https://huggingface.co/datasets/nickh007/protocol-bench)**

### The commercial offering

These are the engine. What is **not** open source is what makes it useful at scale: the maintained
hazard-property corpora, composition analysis that finds hazards existing only when two components
are combined, the trust-model sensitivity sweep, and the evidence trail that makes a verdict auditable
after the fact. The tools above are MIT and stay that way.

## Licence

MIT. See `LICENSE`.
