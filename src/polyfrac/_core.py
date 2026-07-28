"""Exact univariate polynomial + rational-function arithmetic over Q, with Sturm-sequence real-root
counting.

Purpose: solve linear systems whose entries are rational functions of a parameter p -- for example an
absorbing-chain hitting probability P(p) = N(p)/D(p) with exact integer-rational coefficients -- and
then CERTIFY sign conditions of the result over an interval by exact real-root counting, rather than
by sampling the parameter at chosen points.

Classes:
  Poly      -- dense univariate polynomial over Fraction coefficients.
  PolyFrac  -- a ratio of two Poly, kept in lowest terms.

Functions:
  gauss_solve  -- Gaussian elimination over the field of rational functions.
  sturm_chain  -- the (squarefree) Sturm chain of a polynomial.
  count_roots  -- exact count of DISTINCT real roots in a half-open interval (a, b].
  positive_on  -- certify a polynomial has no root in (a, b] and report its sign there.

No floating point is used anywhere in the ARITHMETIC, and none is accepted at the door either.

That second half used to be untrue. `Poly([0.1, 0.2])` was accepted and stored as
``Fraction(3602879701896397, 36028797018963968)`` — the exact value of the *binary approximation* of
0.1, not of 0.1. Every operation after that is exact, and every result is exactly wrong by however
much the input was already off. Exactness you can only rely on when the caller happens to have
supplied exact input is not a guarantee, so a float now raises `InexactInput`.

Use `Fraction(1, 10)`, the string `"0.1"`, or `Poly.from_floats([...])` if a float really is what
you have and you accept the conversion.
"""

from __future__ import annotations

from fractions import Fraction

__all__ = [
    "Poly",
    "PolyFrac",
    "InexactInput",
    "gauss_solve",
    "sturm_chain",
    "count_roots",
    "positive_on",
    "exact",
]


class InexactInput(TypeError):
    """A float (or other inexact value) was supplied where an exact rational was required.

    Raised rather than converted, because the conversion is silent and lossy: a binary float is
    almost never the decimal you wrote, so converting it produces an exact answer to a question
    slightly different from the one you asked.
    """


def exact(x) -> Fraction:
    """Coerce to `Fraction`, refusing anything that is already approximate.

    Accepts int, Fraction, Decimal, and str (``"0.1"`` and ``"1/10"`` both parse exactly).
    Rejects float and complex.
    """
    if isinstance(x, float):
        if x != x or x in (float("inf"), float("-inf")):
            # Non-finite: there is no rational to convert to at all, and building the usual
            # message below would itself raise (Fraction(inf) overflows).
            raise InexactInput(f"refusing {x!r}: it is not a finite value and has no rational representation")
        raise InexactInput(
            f"refusing the float {x!r}: it is a binary approximation, so converting it would give "
            f"an exact answer to the wrong question (it equals {Fraction(x)}). "
            f'Use Fraction({x!r}) if you really mean that value, or the string "{x}" / '
            f"Fraction(numerator, denominator) for the decimal you intended, or Poly.from_floats()."
        )
    if isinstance(x, complex):
        raise InexactInput(f"refusing the complex value {x!r}: this package is real-valued")
    try:
        return Fraction(x)
    except (TypeError, ValueError) as e:
        raise InexactInput(f"cannot interpret {x!r} as an exact rational: {e}") from e


class Poly:
    """Dense univariate polynomial over Fraction: coeffs[i] is the coefficient of x^i."""

    __slots__ = ("c",)

    def __init__(self, coeffs):
        c = [exact(x) for x in coeffs]
        while len(c) > 1 and c[-1] == 0:
            c.pop()
        self.c = c

    @staticmethod
    def const(x) -> Poly:
        return Poly([exact(x)])

    @staticmethod
    def from_floats(coeffs) -> Poly:
        """Build from floats, converting each to its EXACT binary value.

        The explicit opt-in for `Poly([...])`'s refusal. Every guarantee downstream still holds —
        the arithmetic is exact — but it is exact about the binary approximations you passed in, not
        about the decimals you probably meant. ``0.1`` becomes ``3602879701896397/36028797018963968``.
        Prefer strings or `Fraction` when the decimal is what matters.
        """
        return Poly([Fraction(float(x)) for x in coeffs])

    @staticmethod
    def from_roots(roots) -> Poly:
        """The monic polynomial with exactly the given roots: prod_i (x - r_i).

        Repeated entries give repeated roots. Note that `count_roots` counts DISTINCT roots, so
        `Poly.from_roots([2, 2])` has one distinct root in any interval containing 2.
        """
        out = Poly([Fraction(1)])
        for r in roots:
            out = out * Poly([-exact(r), Fraction(1)])
        return out

    @staticmethod
    def x() -> Poly:
        return Poly([0, 1])

    def degree(self) -> int:
        return 0 if self.is_zero() else len(self.c) - 1

    def is_zero(self) -> bool:
        return len(self.c) == 1 and self.c[0] == 0

    def __add__(self, o: Poly) -> Poly:
        n = max(len(self.c), len(o.c))
        return Poly([(self.c[i] if i < len(self.c) else 0) + (o.c[i] if i < len(o.c) else 0) for i in range(n)])

    def __sub__(self, o: Poly) -> Poly:
        n = max(len(self.c), len(o.c))
        return Poly([(self.c[i] if i < len(self.c) else 0) - (o.c[i] if i < len(o.c) else 0) for i in range(n)])

    def __mul__(self, o: Poly) -> Poly:
        out = [Fraction(0)] * (len(self.c) + len(o.c) - 1)
        for i, a in enumerate(self.c):
            if a:
                for j, b in enumerate(o.c):
                    out[i + j] += a * b
        return Poly(out)

    def __neg__(self) -> Poly:
        return Poly([-a for a in self.c])

    def scale(self, k) -> Poly:
        k = exact(k)
        return Poly([a * k for a in self.c])

    def divmod(self, o: Poly) -> tuple[Poly, Poly]:
        """Exact polynomial long division (quotient, remainder) over Q."""
        if o.is_zero():
            raise ZeroDivisionError("poly division by zero")
        r = [Fraction(x) for x in self.c]
        q = [Fraction(0)] * max(1, len(r) - len(o.c) + 1)
        dlead = o.c[-1]
        while len(r) >= len(o.c) and not (len(r) == 1 and r[0] == 0):
            k = r[-1] / dlead
            pos = len(r) - len(o.c)
            q[pos] = k
            for i, b in enumerate(o.c):
                r[pos + i] -= k * b
            while len(r) > 1 and r[-1] == 0:
                r.pop()
            if len(r) < len(o.c):
                break
        return Poly(q), Poly(r)

    def derivative(self) -> Poly:
        if self.degree() == 0:
            return Poly([0])
        return Poly([self.c[i] * i for i in range(1, len(self.c))])

    def eval(self, x) -> Fraction:
        x = exact(x)
        acc = Fraction(0)
        for a in reversed(self.c):
            acc = acc * x + a
        return acc

    def gcd(self, o: Poly) -> Poly:
        a, b = self, o
        while not b.is_zero():
            a, b = b, a.divmod(b)[1]
        if a.is_zero():
            return Poly([0])
        return a.scale(1 / a.c[-1])  # monic normalization

    def __eq__(self, o) -> bool:
        return isinstance(o, Poly) and self.c == o.c

    def __repr__(self) -> str:
        return (
            "Poly(" + " + ".join(f"{a}*x^{i}" for i, a in enumerate(self.c) if a) + ")"
            if not self.is_zero()
            else "Poly(0)"
        )


class PolyFrac:
    """Reduced rational function num/den over Q (den normalized monic-lead positive)."""

    __slots__ = ("num", "den")

    def __init__(self, num: Poly, den: Poly | None = None):
        den = den if den is not None else Poly([1])
        if den.is_zero():
            raise ZeroDivisionError("PolyFrac with zero denominator")
        g = num.gcd(den)
        if not g.is_zero() and g.degree() >= 0 and not (g.degree() == 0 and g.c[0] == 1):
            num = num.divmod(g)[0]
            den = den.divmod(g)[0]
        if den.c[-1] < 0:
            num, den = -num, -den
        self.num, self.den = num, den

    @staticmethod
    def const(x) -> PolyFrac:
        return PolyFrac(Poly.const(x))

    def __add__(self, o: PolyFrac) -> PolyFrac:
        return PolyFrac(self.num * o.den + o.num * self.den, self.den * o.den)

    def __sub__(self, o: PolyFrac) -> PolyFrac:
        return PolyFrac(self.num * o.den - o.num * self.den, self.den * o.den)

    def __mul__(self, o: PolyFrac) -> PolyFrac:
        return PolyFrac(self.num * o.num, self.den * o.den)

    def __truediv__(self, o: PolyFrac) -> PolyFrac:
        if o.num.is_zero():
            raise ZeroDivisionError("PolyFrac division by zero")
        return PolyFrac(self.num * o.den, self.den * o.num)

    def __neg__(self) -> PolyFrac:
        return PolyFrac(-self.num, self.den)

    def is_zero(self) -> bool:
        return self.num.is_zero()

    def eval(self, x) -> Fraction:
        d = self.den.eval(x)
        if d == 0:
            raise ZeroDivisionError(f"pole at {x}")
        return self.num.eval(x) / d

    def __eq__(self, o) -> bool:
        return isinstance(o, PolyFrac) and (self.num * o.den) == (o.num * self.den)

    def __repr__(self) -> str:
        return f"PolyFrac({self.num!r} / {self.den!r})"


def gauss_solve(A: list[list[PolyFrac]], b: list[PolyFrac]) -> list[PolyFrac]:
    """Solve A x = b exactly over the rational-function field (partial 'pivoting' = first nonzero).
    Raises if singular. n is tiny for the DTMCs (<= 6), so plain elimination is exact and fast."""
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = next((r for r in range(col, n) if not M[r][col].is_zero()), None)
        if piv is None:
            raise ValueError("singular system")
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [e / pv for e in M[col]]
        for r in range(n):
            if r != col and not M[r][col].is_zero():
                f = M[r][col]
                M[r] = [M[r][k] - f * M[col][k] for k in range(n + 1)]
    return [M[i][n] for i in range(n)]


def sturm_chain(P: Poly) -> list[Poly]:
    """Canonical Sturm chain of the squarefree part of P (dividing out gcd(P, P') keeps the same
    DISTINCT real roots, which is what root counting needs)."""
    g = P.gcd(P.derivative())
    if g.degree() > 0:
        P = P.divmod(g)[0]
    chain = [P, P.derivative()]
    while not chain[-1].is_zero() and chain[-1].degree() > 0:
        rem = chain[-2].divmod(chain[-1])[1]
        if rem.is_zero():
            break
        chain.append(-rem)
    return [q for q in chain if not q.is_zero()]


def _sign_changes(chain: list[Poly], x) -> int:
    signs = []
    for q in chain:
        v = q.eval(x)
        if v != 0:
            signs.append(1 if v > 0 else -1)
    return sum(1 for i in range(len(signs) - 1) if signs[i] != signs[i + 1])


def count_roots(P: Poly, a, b) -> int:
    """Number of DISTINCT real roots of P in (a, b] (Sturm's theorem; a < b, exact rationals).

    Raises `ValueError` for the zero polynomial. Sturm's theorem does not apply to it, and the
    sign-change arithmetic happens to yield 0 — a confident "no roots here" about a function that
    is zero at every point of the interval. Refusing is the only correct answer available, since
    the true count is not finite.
    """
    a, b = exact(a), exact(b)
    if not a < b:
        raise ValueError(f"need a < b, got a={a} b={b}")
    if P.is_zero():
        raise ValueError(
            "the zero polynomial is zero at every point, so it has infinitely many roots in any "
            "interval; Sturm's theorem does not apply and no finite count exists"
        )
    chain = sturm_chain(P)
    return _sign_changes(chain, a) - _sign_changes(chain, b)


def positive_on(P: Poly, a, b) -> dict:
    """EXACT certificate that P(x) > 0 for all x in [a, b]: P(a) > 0, P(b) > 0, and ZERO roots in
    (a, b] by Sturm — a sign change inside would require a root. Returns the full evidence.

    Raises for the zero polynomial rather than reporting ``positive_on_interval: False`` with a
    root count of 0. That verdict was safe in direction but its evidence was wrong, and this
    function's value is the evidence.
    """
    a, b = exact(a), exact(b)
    if P.is_zero():
        raise ValueError("the zero polynomial is never positive and has no finite root count")
    pa, pb = P.eval(a), P.eval(b)
    n_roots = count_roots(P, a, b)
    return {
        "endpoint_a": str(a),
        "endpoint_b": str(b),
        "P_a": str(pa),
        "P_b": str(pb),
        "n_roots_in_interval": n_roots,
        "positive_on_interval": bool(pa > 0 and pb > 0 and n_roots == 0),
    }
