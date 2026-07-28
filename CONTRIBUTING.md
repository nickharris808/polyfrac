# Contributing to polyfrac

Small, focused, and dependency-free is the point of this package. That shapes what changes are easy to
accept.

## Ground rules

1. **No dependencies.** Standard library only. A pull request that adds a runtime dependency will be
   declined regardless of merit — depend on `polyfrac` from your own package instead.
2. **No floating point.** Every coefficient, endpoint, and returned value is exact. If you need a
   `float` anywhere in the computation path, this is the wrong package for that change.
3. **Univariate scope.** Multivariate sign conditions need cylindrical algebraic decomposition, which
   is deliberately out of scope.

## Getting set up

```
python -m venv .venv && . .venv/bin/activate
pip install -e ".[test]"
pytest
```

## Pull requests

- Add a test that fails before your change and passes after. Tests live in `tests/`.
- Keep the public API in `__all__` explicit; anything not listed there is internal.
- If you change interval semantics, say so loudly — `(a, b]` is left-exclusive and right-inclusive, and
  code depends on that.
- Sign-off by [DCO](https://developercertificate.org/) (`git commit -s`). There is no CLA.

## Reporting a wrong answer

An incorrect root count is the most serious possible bug here. If you find one, please include the
polynomial coefficients as exact fractions and the interval endpoints, so it can go straight into the
test suite.
