#!/usr/bin/env python3
"""Fake formal-series checks only. Never calls registered_source_model."""
import argparse
import ast
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from math import factorial
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
SPEC = importlib.util.spec_from_file_location('reviewed_source_definitions_only', BASE / 'source_models.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def need(test, name):
    if not test:
        raise RuntimeError(name)


def multiply(a, b, n):
    out = [Q(0)] * (n + 1)
    for i, av in enumerate(a):
        for j, bv in enumerate(b):
            if i + j <= n:
                out[i + j] += av * bv
    return tuple(out)


def derivative(a):
    return tuple(j * a[j] for j in range(1, len(a)))


def fake_power_exponential(f, n):
    # Independent finite power expansion: [x^<=n] sum f(x)^j/j!.
    # f(0) is set to zero, so powers beyond n contribute no coefficients.
    exponent = (Q(0),) + tuple(f[1:])
    power = (Q(1),) + (Q(0),) * n
    out = [Q(0)] * (n + 1)
    for j in range(n + 1):
        for k, av in enumerate(power):
            out[k] += av / factorial(j)
        power = multiply(power, exponent, n)
    return tuple(out)


def check_fake_geometry_source_expression(h, center, degree):
    # Execute only the reviewed pure geometry/product coefficient statements
    # using artificial h and an intentionally nonregistered positive center.
    # No authorization, exp, source identifier, or source builder is called.
    tree = ast.parse((BASE / 'source_models.py').read_text())
    builder = next(v for v in tree.body if isinstance(v, ast.FunctionDef)
                   and v.name == 'registered_source_model')
    selected = []
    allowed_names = {'ls', 'l2', 'gamma', 'work_gamma'}
    for node in builder.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id in allowed_names for t in node.targets
        ):
            selected.append(node)
        elif isinstance(node, ast.For) and any(
            isinstance(v, ast.Call) and isinstance(v.func, ast.Attribute)
            and isinstance(v.func.value, ast.Name) and v.func.value.id == 'gamma'
            and v.func.attr == 'append' for v in ast.walk(node)
        ):
            selected.append(node)
    need(len(selected) == 5, 'exact geometry/product statement inventory')
    program = ast.fix_missing_locations(ast.Module(body=selected, type_ignores=[]))
    env = {'Q': Q, 'center': center, 'degree': degree, 'h': h}
    exec(compile(program, '<fake formal geometry/source expressions>', 'exec'), env)
    geometry = tuple(-Q((-1) ** n) / center ** (n + 1) for n in range(degree + 1))
    need(env['ls'] == geometry, 'variable negative inverse-geometry coefficients')
    need(env['l2'] == multiply(geometry, geometry, degree), 'exact L squared convolution')
    first = multiply(multiply(geometry, geometry, degree), h, degree)
    second = multiply(geometry, derivative(h), degree)
    dh2 = derivative(derivative(h))
    expected = tuple(4 * first[j] - 2 * second[j] - dh2[j]
                     for j in range(degree + 1))
    need(tuple(env['gamma']) == expected, 'g second derivative/index/scaling formula')
    need(env['work_gamma'] == multiply(geometry, expected, degree), 'Lg convolution formula')
    return 5


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    checks = 0
    n = 24
    fake_denominators = ((Q(2), Q(3, 7), Q(-1, 5)),
                         (Q(-3, 2), Q(5, 11)),
                         (Q(7), Q(0), Q(1, 13), Q(-2, 17)))
    for d in fake_denominators:
        inv = MODULE.reciprocal_coefficients(d, n)
        need(multiply(d, inv, n) == (Q(1),) + (Q(0),) * n,
             'fake rational reciprocal identity')
        checks += 1
    for sign in (-1, 1):
        f = (Q(5, 13), Q(sign, 3), Q(2, 7), Q(-sign, 11)) + (Q(0),) * 23
        beta = MODULE.normalized_exp_coefficients(f, 26)
        need(beta == fake_power_exponential(f, 26), 'independent fake finite exp-power algebra')
        checks += 1
        # Artificial even-like and signed-like forms; neither is the bump.
        for h in (beta, tuple(Q(2, 9) * beta[j] + (beta[j - 1] if j else 0)
                             for j in range(27))):
            checks += check_fake_geometry_source_expression(h, Q(13, 7), 24)
    for q in (Q(-1), Q(-2, 7), Q(-1, 3), Q(0)):
        lo20, hi20 = MODULE.exp_rational_interval(q, 20)
        lo200, hi200 = MODULE.exp_rational_interval(q, 200)
        need(Q(0) < lo20 <= lo200 <= hi200 <= hi20,
             'fabricated nested rational alternating exponential enclosures')
        checks += 1
    for q in (Q(1, 2), Q(-3, 2), -0.25):
        try:
            MODULE.exp_rational_interval(q)
        except ValueError:
            checks += 1
        else:
            raise RuntimeError('unsupported scalar exp input accepted')
    for q in (Q(-17, 13), Q(0), Q(13, 17)):
        chosen = MODULE.dyadic_point(q, 19)
        need(chosen <= q < chosen + Q(1, 2 ** 19), 'signed exact dyadic floor')
        checks += 1
    previous_even = Q(2737816576, 62462907)
    previous_signed = Q(2321190208, 62462907)
    need(Q(8, 27) * previous_even < 32 and Q(8, 27) * previous_signed < 32,
         'Lg M32 from independently proved source analytic domain')
    checks += 1
    ratio = Q(1, 16)
    need(32 * ratio ** 25 / (1 - ratio) == Q(1, 15 * 2 ** 91),
         'Lg degree24 analytic tail normalization')
    checks += 1
    tree = ast.parse((BASE / 'source_models.py').read_text())
    need(all(isinstance(v, (ast.Expr, ast.Import, ast.ImportFrom, ast.FunctionDef))
             for v in tree.body), 'source model has no module-level source evaluation')
    checks += 1
    receipt = {'status': 'PASS_FAKE_FORMAL_SOURCE_MODEL_REVIEW', 'checks': checks,
               'registered_source_model_calls': 0, 'registered_identifiers_used': 0,
               'physical_source_evaluations': 0, 'scientific_input_arrays': 0,
               'source_model_sha256': hashlib.sha256((BASE / 'source_models.py').read_bytes()).hexdigest(),
               'analytic_bound_Lg': 32, 'Lg_tail_exact': str(Q(1, 15 * 2 ** 91))}
    encoded = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
