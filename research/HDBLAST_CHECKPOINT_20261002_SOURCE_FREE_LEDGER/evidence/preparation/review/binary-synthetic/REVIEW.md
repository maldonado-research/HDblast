# Synthetic exact binary-rational conversion review

**PASS: 103/103 checks in normal Python and 103/103 in optimized Python.** This review used constructed scalars only. It read no scientific archives, saved arrays, or numerical results, and executed no physical calculation or external call. It changes no registration, acceptance threshold, or scientific failure status.

The pinned runtime is `/workspace/hdblast-cloud-setup/venv-frw/bin/python`, NumPy 2.2.6 and mpmath 1.3.0. Its `numpy.longdouble` has 128 storage bits but **64 significand bits** (`nmant = 63`). Storage padding is therefore unsuitable for checking numerical equality.

The successor of one has the exact ratio `(9223372036854775809, 9223372036854775808)`, representing `1 + 2^-63`. The native `as_integer_ratio()` retains this extra bit. Constructing the result in mpmath's default 53-bit context loses it. An explicit context precision of at least `numpy.finfo(numpy.longdouble).nmant + 1` bits is required; the existing 50 and 70 decimal digit contexts exceed that requirement.

Use the following conversion contract:

```python
# x must still be a numpy.longdouble scalar; ctx.prec must be >= 64 here.
if not np.isfinite(x):
    raise ValueError("non-finite source")
n, d = x.as_integer_ratio()
y = ctx.mpf(int(n)) / ctx.mpf(int(d))
```

The denominator is a positive power of two. The source has at most 64 significant bits, so the integer constructors and division are exact at the required precision even for extreme exponents. Do not first use `float(x)`, `complex(x)`, `float64`, or display strings. Maintain sufficient precision for later arithmetic as well; exact input conversion does not make later finite-precision operations exact.

For a `numpy.clongdouble`, obtain `.real` and `.imag` directly and convert each as a longdouble. Both components retain their native type in the pinned runtime. Construct `ctx.mpc(real_mpf, imaginary_mpf)` only after their conversion, with the same sufficient context precision.

The check compares each original integer ratio with the exact finite mpmath representation. Its numerical roundtrip uses the small integer mantissa and `numpy.ldexp`, then compares `as_integer_ratio()` again. This avoids decimal formatting, byte padding comparisons, and binary64 intermediates. The verification reads mpmath's `_mpf_` four-tuple; that detail is valid for the pinned mpmath 1.3.0 runtime and should remain confined to the verifier.

Both signed zeros become mathematical zero in mpmath. Preserve a separate `negative_zero = (n == 0 and numpy.signbit(x))` flag if roundtrip validation requires the source zero sign. The roundtrip in this review restores that flag. This sign convention is an explicit representation choice, not a change to any physical equation.

The checks cover 26 scalar cases at 64, 166 and 236 bits, including values beyond binary64, maximum finite magnitude, the smallest normal, largest and smallest subnormals, signed zero, and a broad exponent grid. They also cover complex components, finite guards, inadequate precision, forbidden source downcasts, and two exact synthetic dyadic arithmetic identities. Exceptions and checks remain active under `python -O`.

A modest benchmark performed 4,000 synthetic conversions, including exact-ratio verification and extreme exponents, in approximately 0.90 seconds per mode. This is indicative conversion overhead only; it does not establish performance for a saved archive or a scientific calculation.

There is no conversion blocker in this pinned runtime. Adoption must keep the source type, finite guard, explicit precision, exact-ratio verifier, and optional signed-zero side channel. The review provides conversion infrastructure evidence only.

Evidence files:

- `review_binary_conversion.py`
- `SYNTHETIC_BINARY_CONVERSION_normal.json`
- `SYNTHETIC_BINARY_CONVERSION_optimized.json`

Source SHA-256: `1c28141fff46e72dafbd4c2304f6f044612f5b5074d7c3d11dc3f7446b61be81`.
