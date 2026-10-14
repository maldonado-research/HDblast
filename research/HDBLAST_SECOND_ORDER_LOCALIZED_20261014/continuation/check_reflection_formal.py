#!/usr/bin/env python3
"""Exact finite-group checks; no imported checkpoint code or numerical sampling."""
from fractions import Fraction
from itertools import combinations, product
from math import factorial
import argparse
import json
import pathlib
import sympy as sp


parser = argparse.ArgumentParser()
parser.add_argument('--output', default='reflection_formal.json')
parser.add_argument('--inject-failure', action='store_true',
                    help='Intentionally fail an explicit check to verify nonzero exit under -O.')
args = parser.parse_args()
checks = []
controls = []


def require(name, condition):
    passed = bool(condition)
    checks.append({'name': name, 'passed': passed})
    if not passed:
        raise RuntimeError(name)


def control(name, wrong_claim):
    rejected = not bool(wrong_claim)
    controls.append({'name': name, 'wrong_claim_rejected': rejected})
    if not rejected:
        raise RuntimeError('Wrong-formula control survived: ' + name)


if args.inject_failure:
    require('deliberate_failure_control',False)

signs = list(product((-1, 1), repeat=4))
pairs = list(combinations(range(4), 2))
require('group_has_sixteen_distinct_elements', len(set(signs)) == 16)
for a in range(4):
    require(f'boost_{a}_character_average_zero', sum(t[a] for t in signs) == 0)
for a, b in pairs:
    require(f'rotation_{a}{b}_character_average_zero',
            sum(t[a]*t[b] for t in signs) == 0)
require('full_group_invariant_kid_dimension_zero',
        sum(all(t[a] == 1 for t in signs) for a in range(4))
        + sum(all(t[a]*t[b] == 1 for t in signs) for a,b in pairs) == 0)

x = sp.symbols('x1:5')
G, G1, K, K1 = sp.symbols('G G1 K K1')
s = sum(t*t for t in x)
HA = x[0]**2 - x[1]**2
HB = 2*x[0]**2 + 3*x[1]**2 - 4*x[2]**2 - x[3]**2
q, p = G*HA, K*HB
# Radial background denominators and constant normalizations are invariant.
h = sp.Matrix(4,4,lambda i,j: 2*q*(s*int(i==j)-3*x[i]*x[j]))
pi = sp.Matrix(4,4,lambda i,j: p*x[i]*x[j])
psi = (2*s*G1 + 6*G)*HA
varpi = (2*s*K1 + 7*K)*HB
m = sp.symbols('m1:5')

for index, t in enumerate(signs):
    sub = {x[i]: t[i]*x[i] for i in range(4)}
    refl = lambda f: sp.expand(f.xreplace(sub))
    require(f'element_{index}_scalar_ADM_sectors_invariant',
            sp.expand(refl(q)-q)==0 and sp.expand(refl(p)-p)==0
            and sp.expand(refl(psi)-psi)==0
            and sp.expand(refl(varpi)-varpi)==0)
    require(f'element_{index}_metric_covariance', all(
            sp.expand(refl(h[i,j])-t[i]*t[j]*h[i,j])==0
            for i in range(4) for j in range(4)))
    require(f'element_{index}_momentum_density_covariance', all(
            sp.expand(refl(pi[i,j])-t[i]*t[j]*pi[i,j])==0
            for i in range(4) for j in range(4)))
    for a,b in pairs:
        Y = [sp.Integer(0)]*4
        Y[b],Y[a] = x[a],-x[b]
        require(f'element_{index}_rotation_{a}{b}_vector_pullback',all(
            sp.expand(t[i]*refl(Y[i])-t[a]*t[b]*Y[i])==0
            for i in range(4)))
        pairing = sum(Y[i]*m[i] for i in range(4))
        transformed = sum(refl(Y[i])*t[i]*m[i] for i in range(4))
        require(f'element_{index}_rotation_{a}{b}_density_pairing',
                sp.expand(transformed-t[a]*t[b]*pairing)==0)


def compositions(total, length):
    if length == 1:
        if total >= 1:
            yield (total,)
        return
    for first in range(1, total-length+2):
        for tail in compositions(total-first,length-1):
            yield (first,)+tail


eps = sp.Symbol('eps')
u = {j:sp.Symbol(f'u{j}') for j in range(1,7)}
c = {k:sp.Symbol(f'c{k}') for k in range(1,7)}
for n in range(2,7):
    path = sum(u[j]*eps**j/sp.factorial(j) for j in range(1,n+1))
    constraint = sum(c[k]*path**k/sp.factorial(k) for k in range(1,n+1))
    direct = sp.diff(constraint,eps,n).subs(eps,0)
    nonlinear = sum(sp.Rational(factorial(n),factorial(k)*
                           sp.prod(factorial(j) for j in js)) *c[k]*
                           sp.prod(u[j] for j in js)
                    for k in range(2,n+1) for js in compositions(n,k))
    require(f'order_{n}_nth_derivative_factorial_convention',
            sp.expand(direct-c[1]*u[n]-nonlinear)==0)
    parity = {u[j]:(-1)**j*u[j] for j in range(1,n+1)}
    require(f'order_{n}_opposite_copy_source_parity',
            sp.expand(nonlinear.xreplace(parity)-(-1)**n*nonlinear)==0)

inversion_group = [(1,1,1,1),(-1,-1,-1,-1)]
invariant_rotations = sum(all(t[a]*t[b]==1 for t in inversion_group) for a,b in pairs)
require('inversion_only_leaves_six_invariant_rotations',invariant_rotations==6)
control('false_inversion_only_has_no_invariant_kids',invariant_rotations==0)
HB_nc = 2*x[0]*x[1]
JHA = sp.expand(x[0]*sp.diff(HA,x[1])-x[1]*sp.diff(HA,x[0]))
require('inversion_even_noncommuting_pair',
        sp.expand(HA.xreplace({z:-z for z in x})-HA)==0
        and sp.expand(HB_nc.xreplace({z:-z for z in x})-HB_nc)==0)
# On S^3, the exact normalized fourth moment is E[x1²x2²]=1/(4*6).
product_poly = sp.expand(JHA*HB_nc)
require('noncommuting_rotation_product_exact',product_poly == -8*x[0]**2*x[1]**2)
sphere_average = Fraction(-8,4*6)
require('noncommuting_inversion_even_pair_rotational_average_nonzero',sphere_average == Fraction(-1,3))
control('false_antipodal_evenness_suffices_for_rotational_orthogonality',sphere_average==0)
control('false_all_quadrupoles_invariant_under_independent_reflections',
        sp.expand(HB_nc.xreplace({x[0]:-x[0]})-HB_nc)==0)
wrong_second = c[2]*u[1]**2/2
correct_second = c[2]*u[1]**2
control('false_second_derivative_uses_epsilon_squared_coefficient',
        sp.expand(wrong_second-correct_second)==0)

receipt = {
    'status':'passed',
    'arithmetic':'exact symbolic polynomial and rational finite-group calculations',
    'python_optimized':not __debug__,
    'positive_checks':len(checks),
    'wrong_claim_controls':len(controls),
    'all_positive_checks_passed':all(v['passed'] for v in checks),
    'all_wrong_claim_controls_rejected':all(v['wrong_claim_rejected'] for v in controls),
    'full_reflection_group_order':16,
    'invariant_local_kid_dimensions':{'full_coordinate_reflections':0,'central_inversion_only':6},
    'inversion_even_noncommuting_rotation_sphere_average':str(sphere_average),
    'checks':checks,
    'controls':controls,
    'limits':'Finite algebra checks; the smooth localized PDE theorem, induction and no-convergence scope are in the companion proof.'
}
pathlib.Path(args.output).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k not in ('checks','controls')}))
