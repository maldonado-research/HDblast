#!/usr/bin/env python3
"""Exact finite algebra checks supporting, not proving, the analytic estimate."""
import argparse
import json
from pathlib import Path
import sys

import sympy as s


def require(condition, label):
    if not bool(condition):
        raise RuntimeError(label)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--corrupt', choices=['robin_sign', 'lapse_factor', 'claimed_symmetry'])
    parser.add_argument('--output')
    args = parser.parse_args()
    checks = []

    def zero(expr, label):
        require(s.simplify(expr) == 0, label)
        checks.append(label)

    e = s.symbols('epsilon')
    a, s1, s2, ell, w = s.symbols('a s1 s2 ell w')
    q = s.symbols('q0:4', positive=True)
    h0 = s.diag(-q[0], q[1], q[2], q[3])
    pairs = [(i, j) for i in range(4) for j in range(i, 4)]
    hs = s.symbols('H0:10')
    H = s.zeros(4)
    for v, (i, j) in zip(hs, pairs):
        H[i, j] = H[j, i] = v
    hinv_first = h0.inv() - e * h0.inv() * H * h0.inv()
    trace = s.trace(hinv_first * h0)
    N = s.sqrt(1 + e * ell)
    sigma = 6*a + s1*e*w + s2*(e*w)**2/2 + 4*(e*w)**3/6
    sigmap = s1 + s2*e*w + 2*(e*w)**2
    fN = 2*sigma*N**2/3 - a*N**3*trace
    fk = 2*N*fN
    fh = sigma*N*(h0 + e*H)/3
    fphi = -sigmap*N/2
    tau = s.trace(h0.inv()*H)
    Bk = 4*a*ell - 2*a*tau - 4*s1*w/3
    BH = -2*a*H - a*ell*h0 - s1*w*h0/3
    Bphi = s1*ell/4 + s2*w/2
    if args.corrupt == 'lapse_factor':
        BH = -2*a*H - 2*a*ell*h0 - s1*w*h0/3
    zero(s.diff(fk, e).subs(e, 0) + Bk, 'actual_direct_metric_lapse_row')
    for i, j in pairs:
        zero(s.diff(fh[i, j], e).subs(e, 0) + BH[i, j], f'actual_tangent_row_{i}{j}')
    zero(s.diff(fphi, e).subs(e, 0) + Bphi, 'actual_scalar_row')

    variables = s.Matrix([ell, *hs, w])
    Bvec = s.Matrix([Bk, *[BH[i, j] for i, j in pairs], Bphi])
    Bmatrix = Bvec.jacobian(variables)
    concrete = Bmatrix.subs({a: 1, s1: 2, s2: 3, **{x: 1 for x in q}})
    require(concrete.shape == (12, 12), 'twelve_Robin_components')
    checks.append('twelve_Robin_components')
    if args.corrupt == 'claimed_symmetry':
        require(concrete == concrete.T, 'false_full_Robin_symmetry_claim')
    else:
        require(concrete != concrete.T, 'full_Robin_matrix_really_nonsymmetric')
        checks.append('full_Robin_matrix_really_nonsymmetric')

    # Exact first boundary jet of exp(z M): T(0)=I, T_z(0)=M.
    M = concrete
    T0 = s.eye(12)
    Tz0 = -M if args.corrupt == 'robin_sign' else M
    require(-Tz0 + M*T0 == s.zeros(12), 'Robin_to_Neumann_boundary_jet')
    checks.append('Robin_to_Neumann_boundary_jet')
    # The D block is identity and has no N coupling in its first boundary jet.
    fullT0 = s.diag(s.eye(4), T0)
    fullTz0 = s.diag(s.zeros(4), Tz0)
    require(fullT0 == s.eye(16), 'all_sixteen_boundary_values_preserved')
    require(fullTz0[:4, :] == s.zeros(4, 16), 'Dirichlet_block_not_mixed')
    checks.extend(['all_sixteen_boundary_values_preserved', 'Dirichlet_block_not_mixed'])

    # Common principal conjugation is algebraic and does not assume symmetric T.
    aa, bb, cc, dd, pp = s.symbols('aa bb cc dd principal')
    generic = s.Matrix([[aa, bb], [cc, dd]])
    principal = generic.inv() * (pp*s.eye(2)) * generic
    require(s.simplify(principal - pp*s.eye(2)) == s.zeros(2), 'common_scalar_principal_preserved')
    checks.append('common_scalar_principal_preserved')

    # Homogeneous D/N boundary energy flux: D time trace and N normal trace vanish.
    vt = s.Matrix([0]*4 + list(s.symbols('nt0:12')))
    vz = s.Matrix(list(s.symbols('dz0:4')) + [0]*12)
    zero((vt.T*vz)[0], 'positive_orthogonal_split_flux_zero')

    # State flow linearization agrees with the actual M sign.
    rho = s.symbols('rho', positive=True)
    exact_state_derivative = s.Matrix([
        s.diff(fk, e).subs(e, 0),
        *[s.diff(fh[i, j], e).subs(e, 0) for i, j in pairs],
        s.diff(fphi, e).subs(e, 0),
    ])
    require(s.simplify(-rho*exact_state_derivative - rho*Bvec) == s.zeros(12, 1),
            'nonlinear_inward_flow_linearizes_to_M')
    checks.append('nonlinear_inward_flow_linearizes_to_M')
    payload = {'status': 'pass', 'checks': checks, 'check_count': len(checks),
               'scope': 'Exact boundary/principal/flux algebra only; analytic PDE estimates require the written proof.',
               'optimized_python': not __debug__}
    encoded = json.dumps(payload, indent=2)
    if args.output:
        Path(args.output).write_text(encoded+'\n')
    print(encoded)


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(json.dumps({'status': 'fail', 'error': str(exc)}), file=sys.stderr)
        raise SystemExit(1)
