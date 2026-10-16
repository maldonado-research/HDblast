#!/usr/bin/env python3
"""Exact arithmetic audit of source-only integer tick interval operations."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
NUMERIC = ROOT.parent/'numeric'


def need(ok, text):
    if not ok:
        raise ValueError(text)


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    need(not args.output.exists(), 'refuse overwrite')
    fast = load(NUMERIC/'fast_dyadic.py', 'independent_tick_engine')
    core = load(NUMERIC/'interval_core.py', 'independent_fraction_reference')
    G = core.GRID
    need(G == fast.GRID == 2**112, 'identical fixed grid')
    points = [-G-1, -G, -G+1, -2, -1, 0, 1, 2, G-1, G, G+1]
    intervals = [(F(points[i], G), F(points[j], G)) for i in range(len(points)) for j in range(i, len(points))]
    intervals += [(F(1, 3), F(7, 11)), (-F(7, 11), -F(1, 3)), (F(-11, 3), F(17, 7))]
    checks = 0
    def compare(actual, expected, label):
        nonlocal checks
        need(actual.row() == expected.row(), label)
        checks += 1
    for lo, hi in intervals:
        a, b = fast.D(lo, hi), core.I(lo, hi)
        compare(a, b, 'dyadic rational constructor')
        compare(fast.D(b), b, 'explicit Fraction interval import')
        compare(fast.D(a), b, 'dyadic copy constructor')
        compare(-a, -b, 'dyadic exact negation')
        for scalar in (0, 3, -7, F(1, 3), -F(7, 11)):
            for actual, expected in ((scalar+a, scalar+b), (a+scalar, b+scalar),
                                     (scalar-a, scalar-b), (a-scalar, b-scalar),
                                     (scalar*a, scalar*b), (a*scalar, b*scalar),
                                     (a+b, b+b), (a-b, b-b), (a*b, b*b)):
                compare(actual, expected, 'integer Fraction or imported interval scalar coercion')
            if scalar:
                compare(a/scalar, b/scalar, 'Fraction denominator coercion')
            if not b.lo <= 0 <= b.hi:
                compare(scalar/a, scalar/b, 'Fraction numerator reflected division')
        for exponent in range(7):
            compare(a**exponent, b**exponent, 'dyadic fixed rounding power')
        if not b.lo <= 0 <= b.hi:
            compare(a.inv(), b.inv(), 'dyadic signed reciprocal')
            compare(a**(-3), b**(-3), 'dyadic negative power')
        for lo2, hi2 in intervals:
            c, d = fast.D(lo2, hi2), core.I(lo2, hi2)
            compare(a+c, b+d, 'dyadic signed addition')
            compare(a-c, b-d, 'dyadic signed subtraction')
            compare(a*c, b*d, 'dyadic four-endpoint multiplication')
            if not d.lo <= 0 <= d.hi:
                compare(a/c, b/d, 'dyadic signed division')
    for lo, hi in ((-1, 1), (0, 1), (-1, 0), (0, 0)):
        try:
            fast.D(lo, hi).inv()
        except ValueError as exc:
            need(str(exc) == 'zero denominator', 'dyadic denominator failure reason')
            checks += 1
        else:
            raise ValueError('dyadic reciprocal accepted zero')
    paths = [NUMERIC/'fast_dyadic.py', NUMERIC/'interval_core.py', Path(__file__)]
    out = {'status': 'PASS_INDEPENDENT_INTEGER_DYADIC_REVIEW', 'exact_endpoint_differential_checks': checks,
           'all_point_proof': 'Intervals have endpoints l/G,h/G. Addition and negation are exact. Product enclosure is floor(min four products/G),ceil(max four products/G) in ticks. For either same-sign denominator, reciprocal endpoints are floor(G^2/h),ceil(G^2/l). Signed floor division and ceil(n/d)=-floor(-n/d) are exact for positive or negative d. Integer-power recursion has precisely the Fraction reference rounding order. Arbitrary-size integers avoid overflow.',
           'coercion_scope': 'D imports an interval object through exact lo/hi Fraction endpoints; D arithmetic accepts integer/Fraction scalars on either side and imported interval objects on its right. The source cell adapter explicitly converts every metric/scalar/background interval to D before evaluation, so no reference-I-left/D-right arithmetic is required.',
           'input_sha256': {str(path.relative_to(ROOT.parent)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}}
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True)+'\n')
    print(out['status'])


if __name__ == '__main__':
    main()
