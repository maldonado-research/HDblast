"""Exact value conversion; storage padding is deliberately not an arithmetic input."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from mpmath.libmp import from_man_exp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@dataclass(frozen=True)
class BinaryValue:
    numerator: int
    exponent: int
    negative_zero: bool = False

    def mpf_tuple(self):
        return from_man_exp(self.numerator, self.exponent)

    def longdouble(self):
        if self.numerator == 0:
            return np.copysign(np.longdouble(0), np.longdouble(-1 if self.negative_zero else 1))
        return np.ldexp(np.longdouble(self.numerator), self.exponent)


def exact_binary(value):
    original = np.asarray(value)
    require(original.ndim == 0 and original.dtype == np.dtype(np.longdouble),
            "Native long-double scalar required before exact conversion")
    value = np.longdouble(value)
    require(np.isfinite(value), "Nonfinite retained real component")
    if value == 0:
        result = BinaryValue(0, 0, bool(np.signbit(value)))
    else:
        numerator, denominator = value.as_integer_ratio()
        require(denominator > 0 and denominator & (denominator - 1) == 0,
                "Retained value is not a binary rational")
        low_bit = abs(numerator) & -abs(numerator)
        shift = low_bit.bit_length() - 1
        result = BinaryValue(numerator >> shift, shift - (denominator.bit_length() - 1))
        require(abs(result.numerator).bit_length() <= np.finfo(np.longdouble).nmant + 1,
                "Exact mantissa exceeds the retained dtype")
    restored = result.longdouble()
    require(restored == value and np.signbit(restored) == np.signbit(value),
            "Exact binary-rational long-double roundtrip failed")
    return result


def exact_complex(value):
    original = np.asarray(value)
    require(original.ndim == 0 and original.dtype == np.dtype(np.clongdouble),
            "Native complex-long-double scalar required before exact conversion")
    return exact_binary(value.real), exact_binary(value.imag)
