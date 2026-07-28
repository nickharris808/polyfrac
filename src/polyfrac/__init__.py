"""polyfrac — exact univariate polynomial and rational-function arithmetic over the rationals,
with Sturm-sequence real-root counting.

Everything is exact. There is no floating point anywhere in this package: coefficients are
``fractions.Fraction``, and root *counts* are integers obtained by counting sign changes in a
Sturm chain rather than by numerical root finding.

Quick start
-----------
>>> from polyfrac import Poly, count_roots, positive_on
>>> p = Poly.from_roots([1, 2, 3])          # (x-1)(x-2)(x-3)
>>> count_roots(p, 0, 4)                    # exactly three real roots in (0, 4]
3
>>> count_roots(p, 0, Fraction(3, 2))       # exactly one, namely x = 1
1
>>> positive_on(Poly([-1, 1]), 2, 5)["positive_on_interval"]   # x - 1 > 0 on (2, 5]
True
"""

from fractions import Fraction

from ._core import Poly, PolyFrac, count_roots, gauss_solve, positive_on, sturm_chain

__all__ = [
    "Poly",
    "PolyFrac",
    "gauss_solve",
    "sturm_chain",
    "count_roots",
    "positive_on",
    "Fraction",
    "__version__",
]
__version__ = "0.1.0"
