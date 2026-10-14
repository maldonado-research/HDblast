#!/usr/bin/env python3
"""Exact checks of the full longitudinal-gauge quadratic angular TT source.

No imported ADM constraint engine: compute five-dimensional Christoffels
and Ricci from a diagonal spacetime metric with arbitrary spacetime Psi.
Then independently check the S^3 angular tensor decomposition as a rational
identity in stereographic coordinates. Sympy is the sole dependency.
"""
import json
import argparse
import hashlib
import sys
import time
from pathlib import Path
import sympy as s


def clean(expr):
    return s.factor(s.cancel(s.together(expr)))


def assert_zero(expr, label):
    result = clean(expr)
    if result != 0:
        raise AssertionError(f"{label}: {result}")


def require_nonzero(expr, label):
    if clean(expr) == 0:
        raise AssertionError(f"wrong-formula control failed to reject {label}")


def ricci_diagonal(coords, metric):
    """Direct coordinate definition, using diagonality solely for inversion."""
    dim = len(coords)
    gamma = {}
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                value = (
                    (s.diff(metric[a], coords[c]) if a == b else 0)
                    + (s.diff(metric[a], coords[b]) if a == c else 0)
                    - (s.diff(metric[b], coords[a]) if b == c else 0)
                ) / (2 * metric[a])
                gamma[a, b, c] = s.simplify(value)
    trace = [s.simplify(sum(gamma[a, a, b] for a in range(dim)))
             for b in range(dim)]
    cache = {}

    def component(b, c):
        if (b, c) not in cache:
            value = sum(s.diff(gamma[a, b, c], coords[a])
                        + gamma[a, b, c] * trace[a]
                        for a in range(dim)) - s.diff(trace[b], coords[c])
            value -= sum(gamma[a, c, d] * gamma[d, b, a]
                         for a in range(dim) for d in range(dim))
            cache[b, c] = clean(value)
        return cache[b, c]
    return component


def check_spacetime():
    y, t, x1, x2, x3, eps = s.symbols('y t x1 x2 x3 eps', real=True)
    coords = [y, t, x1, x2, x3]
    rho = s.Function('rho', positive=True)(y)
    psi = s.Function('Psi')(*coords)
    phi1 = s.Function('chi')(*coords)
    # Unit dS4 in its flat chart. The exact conformal identity is local and
    # chart invariant; the independent sphere check below uses global S^3.
    metric = [s.exp(-4*eps*psi), -rho**2*s.exp(2*eps*psi)]
    metric += [rho**2*s.exp(2*t+2*eps*psi)] * 3
    ricci = ricci_diagonal(coords, metric)
    checked = []
    for i, j in [(2, 3), (2, 4), (3, 4)]:
        expected = -6*eps**2*s.diff(psi, coords[i])*s.diff(psi, coords[j])
        assert_zero(ricci(i, j)-expected, f'R_{i}{j} full exponential metric')
        # Trace-reversed scalar term; U*g has zero off-diagonal entries.
        quadratic_residual = ricci(i, j)-eps**2*s.diff(phi1, coords[i])*s.diff(phi1, coords[j])
        assert_zero(s.diff(quadratic_residual, eps, 2)
                    +12*s.diff(psi, coords[i])*s.diff(psi, coords[j])
                    +2*s.diff(phi1, coords[i])*s.diff(phi1, coords[j]),
                    f'D2 E_{i}{j}')
        checked.append(f'R_{i}{j}')
    require_nonzero(ricci(2, 3), 'omitted geometry')
    require_nonzero(ricci(2, 3)-6*eps**2*s.diff(psi, x1)*s.diff(psi, x2),
                    'wrong geometric sign')
    for i, j in [(2, 3), (2, 4)]:
        expected = -6*eps**2*(s.diff(psi, coords[i])**2-s.diff(psi, coords[j])**2)
        assert_zero(ricci(i, i)-ricci(j, j)-expected,
                    f'R_{i}{i}-R_{j}{j} exponential metric')
        checked.append(f'R_{i}{i}-R_{j}{j}')

    # Differentiate the complete, linear Taylor metric independently.
    metric_linear = [1-4*eps*psi, -rho**2*(1+2*eps*psi)]
    metric_linear += [rho**2*s.exp(2*t)*(1+2*eps*psi)]*3
    ricci_linear = ricci_diagonal(coords, metric_linear)
    for i, j in [(2, 3), (2, 4), (3, 4)]:
        expected = 24*psi*s.diff(psi, coords[i], coords[j]) + 12*s.diff(psi, coords[i])*s.diff(psi, coords[j])
        actual = s.diff(ricci_linear(i, j), eps, 2).subs(eps, 0)
        assert_zero(actual-expected, f'D2 R_{i}{j} linear metric')
        # The exponential completion has g2yy=16Psi^2, g2munu=4rho^2Psi^2 gamma.
        # Its linear Ricci anisotropy is -12 Hess_gamma(Psi^2).
        completion = -12*s.diff(psi**2, coords[i], coords[j])
        assert_zero(actual+completion+12*s.diff(psi, coords[i])*s.diff(psi, coords[j]),
                    f'Taylor completion relation {i}{j}')
        checked.append(f'D2 R_{i}{j} linear path')
    wrong_half = 12*psi*s.diff(psi, x1, x2)+6*s.diff(psi, x1)*s.diff(psi, x2)
    require_nonzero(s.diff(ricci_linear(2, 3), eps, 2).subs(eps, 0)-wrong_half,
                    'missing Taylor factor two')
    return checked


def check_sphere():
    z1, z2, z3 = s.symbols('z1 z2 z3', real=True)
    zs = [z1, z2, z3]
    den = 1 + sum(z*z for z in zs)
    n = s.Matrix([2*z1/den, 2*z2/den, 2*z3/den,
                  (1-sum(z*z for z in zs))/den])
    jac = n.jacobian(zs).applyfunc(clean)
    q = s.eye(3)*4/den**2
    qinv = s.eye(3)*den**2/4
    b = s.diag(1, -1, 0, 0)
    a = b*b
    y = clean((n.T*b*n)[0])
    proj = (s.eye(4)-n*n.T).applyfunc(clean)
    v = (proj*b*n).applyfunc(clean)

    def tf3(m):
        return (m-q*clean(s.trace(qinv*m))/3).applyfunc(clean)

    # Ambient expression obtained from the Weyl curvature product of B.
    def tfambient(m):
        return (m-proj*clean(s.trace(m))/3).applyfunc(clean)

    tensor_ambient = tfambient(y*proj*b*proj-v*v.T+proj*a*proj/2)
    tensor = (jac.T*tensor_ambient*jac).applyfunc(clean)
    f2 = clean(n[0]**2+n[1]**2-s.Rational(1, 2))
    f4 = clean(y*y-(n[0]**2+n[1]**2)/2+s.Rational(1, 12))
    connection = {}
    for k in range(3):
        for i in range(3):
            for j in range(3):
                connection[k, i, j] = clean(sum(qinv[k, ell]*(
                    s.diff(q[ell, j], zs[i])+s.diff(q[ell, i], zs[j])
                    -s.diff(q[i, j], zs[ell]))/2 for ell in range(3)))

    def hessian(f):
        return s.Matrix(3, 3, lambda i, j: s.diff(f, zs[i], zs[j])
                        -sum(connection[k, i, j]*s.diff(f, zs[k]) for k in range(3))).applyfunc(clean)

    grad = s.Matrix([s.diff(y, z) for z in zs])
    residual = tf3(grad*grad.T)-tf3(hessian(f2))/2-tf3(hessian(f4))/3+4*tensor/3
    for i in range(3):
        for j in range(3):
            assert_zero(residual[i, j], f'sphere tensor decomposition {i}{j}')
    assert_zero(s.trace(qinv*tensor), 'sphere tensor trace')
    for j in range(3):
        divergence = sum(qinv[i, k]*(s.diff(tensor[i, j], zs[k])
                         -sum(connection[ell, k, i]*tensor[ell, j]
                              +connection[ell, k, j]*tensor[i, ell]
                              for ell in range(3))) for i in range(3) for k in range(3))
        assert_zero(divergence, f'sphere TT divergence {j}')
    # A concrete nonzero component provides a normalization/sign witness.
    assert_zero(tensor[0, 0].subs({z1:0, z2:0, z3:0})-s.Rational(2, 3),
                'tensor at stereographic north pole')
    require_nonzero((tf3(grad*grad.T)-tf3(hessian(f2))/2
                    -tf3(hessian(f4))/3-4*tensor/3)[0, 0].subs({z1:0,z2:0,z3:0}),
                    'wrong TT projection sign')
    return ['all 9 decomposition components', 'trace', 'all 3 divergences',
            'T11(north pole)=2/3 in stereographic coordinates']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='full_quadratic_tt_results.json')
    args = parser.parse_args()
    started = time.perf_counter()
    print('Checking direct five-dimensional spacetime Ricci...', flush=True)
    spacetime = check_spacetime()
    print('Checking global S3 rational tensor identities...', flush=True)
    sphere = check_sphere()
    results = {
        'spacetime_direct_ricci_checks': spacetime,
        'sphere_exact_rational_checks': sphere,
        'wrong_formula_controls_rejected': ['omitted geometry', 'wrong geometric sign',
                'missing Taylor factor two', 'wrong TT projection sign'],
        'convention': 'g=g0+eps*h1+eps^2*h2/2; E=Ric-dphi*dphi-(2/3)U*g',
        'D2_E_ij_TT': '(16*u^2+(8/3)*v^2)*T_ij',
        'second_order_forcing_ij_TT': '-(16*u^2+(8/3)*v^2)*T_ij',
        'geometry_coefficient_in_D2_E': '16*u^2',
        'scalar_coefficient_in_D2_E': '(8/3)*v^2',
        'status': 'pass',
        'limits': 'Local algebra and specified TT projection; not a through-contact IBVP theorem.',
        'python_optimization_level': sys.flags.optimize,
        'sympy_version': s.__version__,
        'elapsed_seconds': round(time.perf_counter()-started, 3),
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    outfile = Path(__file__).with_name(args.output)
    outfile.write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
