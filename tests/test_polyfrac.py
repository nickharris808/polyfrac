"""Tests for polyfrac.

The whole point of this package is exactness, so the tests are written to fail if any
floating-point contamination or off-by-one in the Sturm sign counting creeps in.
"""

from fractions import Fraction as F

import pytest

from polyfrac import Poly, PolyFrac, count_roots, gauss_solve, positive_on, sturm_chain


# --------------------------------------------------------------------------- Poly basics
def test_from_roots_builds_the_monic_polynomial():
    p = Poly.from_roots([1, 2, 3])
    assert p.degree() == 3
    for r in (1, 2, 3):
        assert p.eval(r) == 0
    assert p.eval(0) == -6  # (0-1)(0-2)(0-3)
    assert p.eval(4) == 6


def test_arithmetic_is_exact_over_rationals():
    third = Poly([F(1, 3)])
    assert (third * Poly([3])).eval(0) == 1  # exactly 1, not 0.9999999999999999
    p = Poly([F(1, 7)]) + Poly([F(2, 7)])
    assert p.eval(0) == F(3, 7)


def test_divmod_and_gcd():
    p = Poly.from_roots([1, 2])
    q, r = p.divmod(Poly([-1, 1]))  # divide by (x - 1)
    assert r.is_zero()
    assert q == Poly([-2, 1])  # x - 2
    g = Poly.from_roots([1, 2]).gcd(Poly.from_roots([2, 3]))
    assert g.eval(2) == 0 and g.degree() == 1  # shared root x = 2


def test_derivative():
    assert Poly([0, 0, 1]).derivative() == Poly([0, 2])  # d/dx x^2 = 2x


# --------------------------------------------------------------------------- Sturm counting
def test_counts_distinct_roots_in_an_interval():
    p = Poly.from_roots([1, 2, 3])
    assert count_roots(p, 0, 4) == 3
    assert count_roots(p, 0, F(3, 2)) == 1
    assert count_roots(p, F(3, 2), F(5, 2)) == 1
    assert count_roots(p, 4, 10) == 0


def test_half_open_convention_is_left_exclusive_right_inclusive():
    """Sturm's theorem as implemented counts roots in (a, b]."""
    p = Poly([-1, 1])  # single root at x = 1
    assert count_roots(p, 1, 2) == 0  # left endpoint excluded
    assert count_roots(p, 0, 1) == 1  # right endpoint included


def test_repeated_roots_count_once():
    p = Poly.from_roots([2, 2, 2])
    assert count_roots(p, 0, 5) == 1  # DISTINCT roots


def test_no_real_roots():
    p = Poly([1, 0, 1])  # x^2 + 1
    assert count_roots(p, -100, 100) == 0


def test_rational_endpoints_are_honoured_exactly():
    """A root at 1/3 must be found by an endpoint at 1/3 and missed by one just below."""
    p = Poly([F(-1, 3), 1])  # x - 1/3
    assert count_roots(p, 0, F(1, 3)) == 1
    assert count_roots(p, 0, F(33, 100)) == 0


def test_sturm_chain_is_squarefree_and_terminates():
    chain = sturm_chain(Poly.from_roots([1, 1, 2]))
    assert len(chain) >= 2
    assert chain[0].degree() >= chain[-1].degree()


def test_requires_ordered_interval():
    with pytest.raises(ValueError):
        count_roots(Poly([-1, 1]), 5, 2)


# --------------------------------------------------------------------------- positive_on
def test_positive_certificate_is_full_evidence():
    cert = positive_on(Poly([-1, 1]), 2, 5)  # x - 1 on (2, 5]
    assert cert["positive_on_interval"] is True
    assert cert["n_roots_in_interval"] == 0
    assert cert["P_a"] == "1" and cert["P_b"] == "4"


def test_positive_certificate_rejects_a_sign_change():
    cert = positive_on(Poly([-1, 1]), 0, 5)  # root at 1 lies inside
    assert cert["positive_on_interval"] is False
    assert cert["n_roots_in_interval"] == 1


def test_positive_certificate_rejects_a_negative_interval():
    cert = positive_on(Poly([-1, 1]), -5, 0)  # x - 1 < 0 throughout
    assert cert["positive_on_interval"] is False


# --------------------------------------------------------------------------- rational functions
def test_polyfrac_reduces_to_lowest_terms():
    f = PolyFrac(Poly.from_roots([1, 2]), Poly.from_roots([2]))  # (x-1)(x-2) / (x-2)
    assert f.num.degree() == 1  # reduced to (x-1)
    assert f.den.degree() == 0


def test_gauss_solve_over_rational_functions():
    """Solve a 2x2 system whose entries are rational functions of the parameter.

        [ x  1 ] [u]   [1]
        [ 0  1 ] [v] = [x]

    so v = x and u = (1 - x) / x.
    """
    one = PolyFrac(Poly([1]), Poly([1]))
    zero = PolyFrac(Poly([0]), Poly([1]))
    X = PolyFrac(Poly([0, 1]), Poly([1]))
    A = [[X, one], [zero, one]]
    b = [one, X]
    u, v = gauss_solve(A, b)
    # v == x
    assert v.num.eval(3) / v.den.eval(3) == 3
    # u == (1 - x)/x  ->  at x = 3 that is -2/3
    assert u.num.eval(3) / u.den.eval(3) == F(-2, 3)


def test_gauss_solve_matches_a_known_absorbing_chain():
    """A two-state absorbing chain: from s0, fail with prob p, else absorb.
    Hitting probability of the failure state is exactly p."""
    one = PolyFrac(Poly([1]), Poly([1]))
    P = PolyFrac(Poly([0, 1]), Poly([1]))  # the parameter p
    # (1)*h = p   =>   h = p
    (h,) = gauss_solve([[one]], [P])
    for val in (F(1, 10), F(1, 2), F(9, 10)):
        assert h.num.eval(val) / h.den.eval(val) == val


# --------------------------------------------------------------------------- the worked example
def test_worked_example_certifies_a_guard_over_a_whole_interval():
    """The README example: a guarded procedure with failure probability p/500 meets a target of
    1/20 for EVERY p in (0, 1], while an unguarded baseline of 3p/5 does not."""
    target = F(1, 20)
    guarded = Poly([0, F(1, 500)])  # p/500
    baseline = Poly([0, F(3, 5)])  # 3p/5

    # target - guarded(p) > 0 on (0, 1]  =>  guarded always meets the target
    margin = Poly([target]) - guarded
    assert positive_on(margin, 0, 1)["positive_on_interval"] is True

    # the baseline crosses the target inside the interval: exactly one root
    baseline_margin = Poly([target]) - baseline
    assert count_roots(baseline_margin, 0, 1) == 1
    # and the crossing is the exact rational 1/12
    assert baseline_margin.eval(F(1, 12)) == 0
