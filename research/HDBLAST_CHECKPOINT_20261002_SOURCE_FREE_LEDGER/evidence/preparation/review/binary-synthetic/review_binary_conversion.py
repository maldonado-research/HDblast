#!/usr/bin/env python3
"""Synthetic-only checks of lossless NumPy longdouble -> mpmath conversion.

No scientific archives, arrays, results, or external services are read.  All
inputs below are constructed synthetic scalars.  Source type checks prevent
an accidental binary64 conversion before the exact integer ratio is obtained.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time

import mpmath as mp
import numpy as np


OUTPUT = Path(__file__).resolve().parent
SOURCE_BITS = int(np.finfo(np.longdouble).nmant) + 1


@dataclass(frozen=True)
class ExactReal:
    value: object
    negative_zero: bool
    numerator: int
    denominator: int


def mp_ratio(value) -> tuple[int, int]:
    """Check exact mpf value using the representation of pinned mpmath 1.3.0.

    The four-tuple is used for verification; the conversion itself uses only
    mpf integer constructors and division.  It must represent a finite dyadic.
    """
    sign, mantissa, exponent, bitcount = value._mpf_
    mantissa, exponent, bitcount = int(mantissa), int(exponent), int(bitcount)
    if mantissa == 0 and exponent != 0:
        raise ValueError("Non-finite mpmath value")
    if mantissa.bit_length() != bitcount:
        raise ValueError("Malformed mpmath mantissa")
    if mantissa == 0:
        return (0, 1)
    numerator = -mantissa if sign else mantissa
    if exponent >= 0:
        return (numerator << exponent, 1)
    return (numerator, 1 << -exponent)


def exact_real(value, context) -> ExactReal:
    """Convert a finite longdouble without decimal formatting or downcasting."""
    if not isinstance(value, np.longdouble):
        raise TypeError("Source must remain a numpy.longdouble scalar")
    if not bool(np.isfinite(value)):
        raise ValueError("Source must be finite")
    if int(context.prec) < SOURCE_BITS:
        raise ValueError("Context precision is below source significand precision")
    numerator, denominator = value.as_integer_ratio()
    numerator, denominator = int(numerator), int(denominator)
    if denominator <= 0 or denominator & (denominator - 1):
        raise ValueError("Source ratio must have a positive power-of-two denominator")
    converted = context.mpf(numerator) / context.mpf(denominator)
    if mp_ratio(converted) != (numerator, denominator):
        raise ValueError("Exact ratio was not preserved")
    negative_zero = numerator == 0 and bool(np.signbit(value))
    return ExactReal(converted, negative_zero, numerator, denominator)


def restore_longdouble(exact: ExactReal) -> np.longdouble:
    """Reconstruct via a <=64-bit mantissa and ldexp, then compare ratios.

    No bytes comparison is made: x87 longdouble arrays can have padding bytes.
    No float(), complex(), decimal representation, or binary64 intermediate
    is used.  The input is an exact conversion returned by exact_real.
    """
    sign, mantissa, exponent, bitcount = exact.value._mpf_
    if int(bitcount) > SOURCE_BITS:
        raise ValueError("Value is not representable by source significand")
    with np.errstate(under="ignore"):
        restored = np.ldexp(np.longdouble(int(mantissa)), int(exponent))
    if sign:
        restored = -restored
    if exact.numerator == 0:
        restored = np.copysign(restored, np.longdouble(-1 if exact.negative_zero else 1))
    if restored.as_integer_ratio() != (exact.numerator, exact.denominator):
        raise ValueError("Restored source ratio differs")
    if exact.numerator == 0 and bool(np.signbit(restored)) != exact.negative_zero:
        raise ValueError("Signed-zero side channel was not preserved")
    return restored


def make_synthetic_complex(real, imaginary) -> np.clongdouble:
    value = np.empty((), dtype=np.clongdouble)
    value.real[()] = np.longdouble(real)
    value.imag[()] = np.longdouble(imaginary)
    return value[()]


def exact_complex(value, context):
    if not isinstance(value, np.clongdouble):
        raise TypeError("Source must remain a numpy.clongdouble scalar")
    real = exact_real(value.real, context)
    imaginary = exact_real(value.imag, context)
    converted = context.mpc(real.value, imaginary.value)
    if mp_ratio(converted.real) != (real.numerator, real.denominator):
        raise ValueError("Complex real ratio differs")
    if mp_ratio(converted.imag) != (imaginary.numerator, imaginary.denominator):
        raise ValueError("Complex imaginary ratio differs")
    return converted, real, imaginary


def run_checks() -> dict:
    checks = []

    def check(label, condition, detail=None):
        record = {"name": label, "passed": bool(condition)}
        if detail is not None:
            record["detail"] = detail
        checks.append(record)
        if not condition:
            raise RuntimeError(f"Check failed: {label}")

    def reject(label, function, expected_exception):
        try:
            function()
        except expected_exception:
            check(label, True)
        else:
            check(label, False)

    context = mp.mp.clone()
    one = np.longdouble(1)
    zero = np.longdouble(0)
    extra_bit = np.nextafter(one, np.longdouble(2))
    check("Pinned significand has 64 bits", SOURCE_BITS == 64)
    check("Beyond-binary64 exact ratio", extra_bit.as_integer_ratio() == ((1 << 63) + 1, 1 << 63))
    context.prec = 53
    numerator, denominator = extra_bit.as_integer_ratio()
    truncated = context.mpf(numerator) / context.mpf(denominator)
    check("Default 53-bit context demonstrably loses source bit", mp_ratio(truncated) != (numerator, denominator))
    reject("Converter refuses inadequate context", lambda: exact_real(extra_bit, context), ValueError)

    info = np.finfo(np.longdouble)
    scalars = {
        "positive_zero": zero,
        "negative_zero": np.longdouble("-0"),
        "one": one,
        "negative_one": -one,
        "successor_of_one": extra_bit,
        "predecessor_of_one": np.nextafter(one, zero),
        "maximum_64bit_integer": np.longdouble((1 << 64) - 1),
        "tiny_normal": info.tiny,
        "largest_subnormal": np.nextafter(info.tiny, zero),
        "smallest_subnormal": info.smallest_subnormal,
        "negative_smallest_subnormal": -info.smallest_subnormal,
        "second_subnormal": np.nextafter(info.smallest_subnormal, one),
        "maximum_finite": info.max,
        "negative_maximum_finite": -info.max,
        "eps": info.eps,
        "epsneg": info.epsneg,
    }
    for exponent in (-16445, -16382, -1024, -63, -1, 0, 1, 63, 1024, 16320):
        with np.errstate(under="ignore"):
            scalars[f"power_two_{exponent}"] = np.ldexp(one, exponent)

    scalar_counts = {}
    for precision in (64, 166, 236):
        context.prec = precision
        for label, scalar in scalars.items():
            result = exact_real(scalar, context)
            restored = restore_longdouble(result)
            check(f"Exact ratio and numerical roundtrip: {label}; precision {precision}",
                  mp_ratio(result.value) == scalar.as_integer_ratio()
                  and restored.as_integer_ratio() == scalar.as_integer_ratio())
        scalar_counts[str(precision)] = len(scalars)

    context.dps = 50
    for label, bad_value in (("NaN", np.longdouble("nan")), ("positive infinity", np.longdouble("inf")), ("negative infinity", np.longdouble("-inf"))):
        reject(f"Finite guard rejects {label}", lambda bad_value=bad_value: exact_real(bad_value, context), ValueError)
    reject("Source guard rejects Python float", lambda: exact_real(1.0, context), TypeError)
    reject("Source guard rejects binary64 scalar", lambda: exact_real(np.float64(1), context), TypeError)

    complex_values = {
        "both_components_beyond_binary64": make_synthetic_complex(extra_bit, np.nextafter(-one, np.longdouble(-2))),
        "extreme_components": make_synthetic_complex(info.max, -info.smallest_subnormal),
        "both_signed_zero": make_synthetic_complex(np.longdouble("-0"), np.longdouble("-0")),
    }
    for label, value in complex_values.items():
        check(f"Complex component types stay longdouble: {label}",
              isinstance(value.real, np.longdouble) and isinstance(value.imag, np.longdouble))
        converted, real, imaginary = exact_complex(value, context)
        check(f"Complex exact component ratios: {label}",
              mp_ratio(converted.real) == value.real.as_integer_ratio()
              and mp_ratio(converted.imag) == value.imag.as_integer_ratio())
        restored = make_synthetic_complex(restore_longdouble(real), restore_longdouble(imaginary))
        check(f"Complex component numerical roundtrip with zero signs: {label}",
              restored.real.as_integer_ratio() == value.real.as_integer_ratio()
              and restored.imag.as_integer_ratio() == value.imag.as_integer_ratio()
              and bool(np.signbit(restored.real)) == bool(np.signbit(value.real))
              and bool(np.signbit(restored.imag)) == bool(np.signbit(value.imag)))
    for label, bad_value in (("nonfinite real", make_synthetic_complex("nan", one)), ("nonfinite imaginary", make_synthetic_complex(one, "inf"))):
        reject(f"Complex finite guard rejects {label}", lambda bad_value=bad_value: exact_complex(bad_value, context), ValueError)
    reject("Complex source guard rejects Python complex", lambda: exact_complex(complex(1, 1), context), TypeError)
    reject("Complex source guard rejects binary64 complex scalar", lambda: exact_complex(np.complex128(1 + 1j), context), TypeError)

    # This is an arithmetic identity on constructed values, not a physical
    # residual or a test against the scientific acceptance thresholds.
    context.dps = 50
    converted = exact_real(extra_bit, context).value
    check("Synthetic arithmetic retains source perturbation", converted - context.mpf(1) == context.ldexp(context.mpf(1), -63))
    check("Synthetic dyadic multiplication identity", converted * context.ldexp(context.mpf(1), 63) == context.mpf((1 << 63) + 1))
    check("Signed zero normalizes in mpf but side channel records sign", mp_ratio(exact_real(np.longdouble("-0"), context).value) == (0, 1) and exact_real(np.longdouble("-0"), context).negative_zero)

    # Modest indicative performance measurement; no scientific input is read.
    iterations = 4000
    benchmark_values = tuple(scalars.values())
    start = time.perf_counter()
    for index in range(iterations):
        exact_real(benchmark_values[index % len(benchmark_values)], context)
    elapsed = time.perf_counter() - start
    return {
        "schema_version": 1,
        "status": "PASS_SYNTHETIC_BINARY_RATIONAL_CONVERSION_GUARDS",
        "prepared_utc": datetime.now(timezone.utc).isoformat(),
        "scope": {
            "synthetic_only": True,
            "saved_scientific_arrays_or_results_loaded": False,
            "physical_calculations_executed": False,
            "external_calls": False,
            "source_or_repository_files_changed": False,
            "scientific_status_changed": False,
        },
        "runtime": {
            "python": sys.version,
            "executable": sys.executable,
            "numpy": np.__version__,
            "mpmath": mp.__version__,
            "longdouble_storage_bits": info.bits,
            "longdouble_significand_bits": SOURCE_BITS,
            "longdouble_minexp": info.minexp,
            "longdouble_maxexp": info.maxexp,
            "optimizer_enabled": not __debug__,
        },
        "checks_passed": sum(record["passed"] for record in checks),
        "checks_total": len(checks),
        "checks": checks,
        "scalar_cases_per_precision": scalar_counts,
        "benchmark": {
            "iterations": iterations,
            "precision_dps": context.dps,
            "elapsed_seconds": elapsed,
            "conversions_per_second": iterations / elapsed,
            "synthetic_only": True,
            "includes_exact_ratio_verification": True,
            "interpretation": "Indicative conversion overhead; no archive-wide performance extrapolation.",
        },
        "requirements": [
            "Keep source scalar as numpy.longdouble; split numpy.clongdouble into native .real and .imag components.",
            "Reject NaN and infinities before calling as_integer_ratio.",
            "Use at least numpy.finfo(longdouble).nmant + 1 bits of mpmath precision, with existing 50/70 decimal digit contexts satisfying the 64-bit requirement.",
            "Construct mpf from exact integer numerator and denominator; do not pass the source through float(), complex(), float64, decimal display formatting, or padding bytes.",
            "Verify exact mpf rational representation against the original as_integer_ratio pair; numerical roundtrip can use a small exact mantissa and numpy.ldexp.",
            "Retain a separate negative-zero flag if bit-sign-sensitive roundtrip validation is required; mpf treats both signed zeros as mathematical zero.",
            "Treat this as conversion infrastructure evidence only, with no change to original registration, scientific thresholds, or failure receipts.",
        ],
        "blockers": [],
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def main():
    report = run_checks()
    suffix = "optimized" if not __debug__ else "normal"
    destination = OUTPUT / f"SYNTHETIC_BINARY_CONVERSION_{suffix}.json"
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "checks_passed": report["checks_passed"], "checks_total": report["checks_total"], "report": str(destination), "source_sha256": report["source_sha256"], "benchmark": report["benchmark"]}, indent=2))


if __name__ == "__main__":
    main()
