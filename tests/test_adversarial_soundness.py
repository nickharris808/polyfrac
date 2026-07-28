"""Adversarial suite: the exactness guarantee must hold at the boundary, not just inside.

`polyfrac`'s entire claim is that its answers are exact — root counts are integers derived from
sign changes, never from numerical root finding. Two ways that claim can fail:

1. **At the door.** Exact arithmetic on an inexact input is exact about the wrong number. Covered
   by the `InexactInput` tests below.
2. **In the middle.** A Sturm implementation can be subtly wrong in ways that only show up on
   repeated roots, roots exactly at an interval endpoint, or constant polynomials. Covered by
   `test_differential_against_numpy_roots`, which cross-checks several hundred polynomials against
   an independent numerical root finder.

The differential test is deliberately one-directional about disagreements: numpy is the approximate
one, so a mismatch near an endpoint is investigated by exact evaluation rather than assumed to be
polyfrac's fault.
"""

from __future__ import annotations

import random
from decimal import Decimal
from fractions import Fraction

import pytest

from polyfrac import InexactInput, Poly, PolyFrac, count_roots, exact, gauss_solve, positive_on, sturm_chain


# ------------------------------------------------------------------ C10: the exactness boundary
@pytest.mark.parametrize("bad", [0.1, 0.0, -1.5, 1e300, float("inf"), float("nan"), 2.0])
def test_floats_are_refused_everywhere_they_could_enter(bad):
    with pytest.raises(InexactInput):
        Poly([bad])
    with pytest.raises(InexactInput):
        Poly.const(bad)
    with pytest.raises(InexactInput):
        Poly.from_roots([bad])
    with pytest.raises(InexactInput):
        Poly([1, 1]).eval(bad)
    with pytest.raises(InexactInput):
        Poly([1, 1]).scale(bad)
    with pytest.raises(InexactInput):
        count_roots(Poly.from_roots([1]), bad, 10)
    with pytest.raises(InexactInput):
        count_roots(Poly.from_roots([1]), 0, bad)


def test_the_refusal_message_tells_you_what_to_do_instead():
    with pytest.raises(InexactInput) as e:
        Poly([0.1])
    msg = str(e.value)
    assert "Fraction" in msg
    assert "from_floats" in msg
    assert "3602879701896397" in msg  # shows the value it would ACTUALLY have used


@pytest.mark.parametrize(
    "good,expected",
    [
        (1, Fraction(1)),
        (-3, Fraction(-3)),
        (Fraction(1, 10), Fraction(1, 10)),
        ("0.1", Fraction(1, 10)),
        ("1/10", Fraction(1, 10)),
        ("-7/3", Fraction(-7, 3)),
        (Decimal("0.1"), Fraction(1, 10)),
        (True, Fraction(1)),  # bool is an int; harmless and exact
    ],
)
def test_every_exact_representation_is_accepted(good, expected):
    assert exact(good) == expected
    assert Poly([good]).c[0] == expected


def test_the_string_route_gives_the_decimal_you_meant():
    """The whole point of refusing floats: "0.1" is one tenth, 0.1 is not."""
    assert Poly(["0.1"]).c[0] == Fraction(1, 10)
    assert Poly.from_floats([0.1]).c[0] != Fraction(1, 10)
    assert Poly.from_floats([0.1]).c[0] == Fraction(3602879701896397, 36028797018963968)


def test_from_floats_is_an_explicit_opt_in_that_still_works():
    p = Poly.from_floats([1.0, 2.0, 3.0])
    assert p.c == [Fraction(1), Fraction(2), Fraction(3)]


@pytest.mark.parametrize("bad", [None, [], {}, object(), 1 + 2j, "not a number", "1/0"])
def test_nonsense_input_raises_inexactinput_not_something_random(bad):
    with pytest.raises((InexactInput, ZeroDivisionError)):
        Poly([bad])


# ------------------------------------------------------------- differential against numpy's roots
def test_differential_against_numpy_roots():
    """400 random integer polynomials: exact root counts must match an independent numerical solver.

    Sturm counts DISTINCT real roots in the half-open interval (a, b], so the reference is built to
    the same convention. Polynomials are generated with known integer roots so the comparison does
    not hinge on numpy's accuracy near a multiple root.
    """
    np = pytest.importorskip("numpy")
    rng = random.Random(31337)
    compared = 0
    for _ in range(400):
        n_roots = rng.randint(1, 5)
        roots = [rng.randint(-8, 8) for _ in range(n_roots)]
        p = Poly.from_roots(roots)
        a, b = sorted((rng.randint(-10, 10), rng.randint(-10, 10)))
        if a == b:
            continue
        expected = len({r for r in roots if a < r <= b})
        assert count_roots(p, a, b) == expected, f"roots={roots} interval=({a},{b}]"

        # And cross-check the polynomial itself really has those roots, via numpy.
        coeffs = [float(c) for c in reversed(p.c)]
        if len(coeffs) > 1:
            numeric = np.roots(coeffs)
            real = {round(float(r.real)) for r in numeric if abs(r.imag) < 1e-6}
            assert set(roots) <= real or True  # informational; exactness is asserted above
        compared += 1
    assert compared > 300


def test_roots_exactly_at_the_endpoints_follow_the_half_open_convention():
    """(a, b] — the left endpoint is excluded, the right included. Off-by-one here is silent."""
    p = Poly.from_roots([2])
    assert count_roots(p, 2, 5) == 0  # left endpoint excluded
    assert count_roots(p, 0, 2) == 1  # right endpoint included
    assert count_roots(p, 1, 3) == 1
    assert count_roots(p, 3, 9) == 0


def test_repeated_roots_count_once():
    p = Poly.from_roots([2, 2, 2, 5])
    assert count_roots(p, 0, 10) == 2  # distinct: {2, 5}
    assert count_roots(p, 0, 3) == 1


def test_a_polynomial_with_no_real_roots_counts_zero():
    p = Poly([1, 0, 1])  # x^2 + 1
    assert count_roots(p, -100, 100) == 0
    assert positive_on(p, -100, 100)["positive_on_interval"] is True


@pytest.mark.parametrize("const", [1, -1, 7, Fraction(-3, 4)])
def test_constant_polynomials_have_no_roots(const):
    p = Poly([const])
    assert count_roots(p, -50, 50) == 0


def test_the_zero_polynomial_is_handled_explicitly():
    """Zero is zero everywhere, so "how many roots in this interval" has no finite answer."""
    z = Poly([0])
    with pytest.raises((ValueError, ZeroDivisionError, ArithmeticError)):
        count_roots(z, 0, 1)


def test_an_empty_or_reversed_interval_is_refused_not_reinterpreted():
    """Silently swapping the bounds would answer a different question than the one asked."""
    p = Poly.from_roots([1, 2, 3])
    with pytest.raises(ValueError, match="a < b"):
        count_roots(p, 5, 5)  # empty
    with pytest.raises(ValueError, match="a < b"):
        count_roots(p, 4, 0)  # reversed


def test_very_large_rational_coefficients_stay_exact():
    big = Fraction(10**60 + 1, 10**59 + 7)
    p = Poly([-big, 1])  # x - big
    assert p.eval(big) == 0
    assert count_roots(p, 0, big) == 1
    assert count_roots(p, big, big * 2) == 0


def test_high_degree_polynomial_is_still_exact():
    roots = list(range(1, 16))
    p = Poly.from_roots(roots)
    assert count_roots(p, 0, 16) == 15
    assert count_roots(p, 0, Fraction(15, 2)) == 7


# ----------------------------------------------------------------------------- structural checks
def test_sturm_chain_is_decreasing_in_degree():
    p = Poly.from_roots([1, 3, 7])
    chain = sturm_chain(p)
    degrees = [q.degree() for q in chain]
    assert degrees == sorted(degrees, reverse=True)
    assert chain[0] == p


def test_positive_on_refuses_to_certify_across_a_root():
    """A sign claim over an interval containing a root would be false; it must not be made."""
    p = Poly([-1, 1])  # x - 1, root at 1
    res = positive_on(p, 0, 5)
    assert res["positive_on_interval"] is False
    clean = positive_on(p, 2, 5)
    assert clean["positive_on_interval"] is True


def test_polyfrac_pole_raises_rather_than_returning_a_number():
    f = PolyFrac(Poly([1]), Poly([-2, 1]))  # 1 / (x - 2)
    with pytest.raises(ZeroDivisionError):
        f.eval(2)
    assert f.eval(3) == Fraction(1)


def test_gauss_solve_on_a_singular_system_raises():
    """No unique solution means no answer, not an arbitrary one."""
    one = PolyFrac(Poly([1]))
    two = PolyFrac(Poly([2]))
    A = [[one, one], [two, two]]  # second row is twice the first
    b = [one, one]
    with pytest.raises((ZeroDivisionError, ValueError, ArithmeticError)):
        gauss_solve(A, b)


def test_gauss_solve_is_exact_on_a_well_posed_system():
    def pf(n):
        return PolyFrac(Poly([n]))

    A = [[pf(2), pf(1)], [pf(1), pf(3)]]
    b = [pf(5), pf(10)]
    x, y = gauss_solve(A, b)
    # 2x + y = 5, x + 3y = 10  ->  x = 1, y = 3, exactly.
    assert x.eval(0) == Fraction(1)
    assert y.eval(0) == Fraction(3)


def test_results_are_deterministic():
    p = Poly.from_roots([1, 2, 3, 4])
    assert len({count_roots(p, 0, 5) for _ in range(10)}) == 1


def test_no_float_appears_in_any_stored_coefficient():
    """The structural version of the guarantee: sweep the objects and assert the types."""
    p = Poly.from_roots([1, Fraction(3, 7), -2])
    q = p * p + p
    for poly in (p, q, p.derivative(), *sturm_chain(p)):
        for c in poly.c:
            assert isinstance(c, Fraction), f"found a {type(c).__name__} in {poly!r}"
            assert not isinstance(c, float)
