#!/usr/bin/env python3
"""Exact local coordinate identities, using only rational Laurent polynomials.

No PDE solve, numerical evolution, interval enclosure, or continuation claim.
Run with: python3 -I -S -B -O verify_proper_clock.py
The documented domain has p=F'(u)>0, q=F'(v)>0, L=exp(B_b)>0.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json


NAMES = ('p', 'q', 'L', 'w', 'r', 's', 'Bt', 'Bz', 'At', 'Az',
         'ft', 'fz', 'E', 'sig', 'sigp', 'Ft', 'Att', 'ftt')
NV = len(NAMES)
ZERO = (0,) * NV


class Laurent:
    def __init__(self, terms=None):
        self.terms = {e: Q(c) for e, c in (terms or {}).items() if c}

    @staticmethod
    def constant(v):
        return Laurent({ZERO: Q(v)})

    @staticmethod
    def variable(name):
        e = list(ZERO)
        e[NAMES.index(name)] = 1
        return Laurent({tuple(e): Q(1)})

    def __add__(self, other):
        if not isinstance(other, Laurent):
            other = Laurent.constant(other)
        terms = self.terms.copy()
        for e, c in other.terms.items():
            terms[e] = terms.get(e, Q(0)) + c
        return Laurent(terms)

    __radd__ = __add__

    def __neg__(self):
        return Laurent({e: -c for e, c in self.terms.items()})

    def __sub__(self, other):
        return self + (-other if isinstance(other, Laurent) else -Q(other))

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if not isinstance(other, Laurent):
            other = Laurent.constant(other)
        terms = {}
        for e, c in self.terms.items():
            for f, d in other.terms.items():
                k = tuple(a + b for a, b in zip(e, f))
                terms[k] = terms.get(k, Q(0)) + c * d
        return Laurent(terms)

    __rmul__ = __mul__

    def __pow__(self, n):
        if not isinstance(n, int):
            raise TypeError('Only integer Laurent powers are supported')
        if n < 0:
            if len(self.terms) != 1:
                raise ValueError('Negative powers require one nonzero monomial')
            e, c = next(iter(self.terms.items()))
            return Laurent({tuple(k * n for k in e): c ** n})
        ans = Laurent.constant(1)
        for _ in range(n):
            ans = ans * self
        return ans

    def rename(self, mapping):
        terms = {}
        for e, c in self.terms.items():
            f = [0] * NV
            for name, exponent in zip(NAMES, e):
                f[NAMES.index(mapping.get(name, name))] += exponent
            key = tuple(f)
            terms[key] = terms.get(key, Q(0)) + c
        return Laurent(terms)

    def substitute_monomials(self, mapping):
        ans = Laurent.constant(0)
        for e, c in self.terms.items():
            term = Laurent.constant(c)
            for name, exponent in zip(NAMES, e):
                base = mapping.get(name, Laurent.variable(name))
                term = term * base**exponent
            ans = ans + term
        return ans

    def zero(self):
        return not self.terms


V = {name: Laurent.variable(name) for name in NAMES}
p, q, L, w, r, s = [V[n] for n in ('p', 'q', 'L', 'w', 'r', 's')]
Bt, Bz, At, Az, ft, fz, E, sig, sigp, Ft = [V[n] for n in
    ('Bt', 'Bz', 'At', 'Az', 'ft', 'fz', 'E', 'sig', 'sigp', 'Ft')]
Att, ftt = V['Att'], V['ftt']
checks = []


def equal(name, lhs, rhs=0):
    residual = lhs - rhs
    if not residual.zero():
        raise RuntimeError('Exact identity failed: ' + name)
    checks.append({'name': name, 'result': 'EXACT_ZERO', 'remaining_terms': 0})


def nonzero(name, expr):
    if expr.zero():
        raise RuntimeError('Expected nonidentity vanished: ' + name)
    checks.append({'name': name, 'result': 'EXACT_NONZERO',
                   'remaining_terms': len(expr.terms)})


# J maps (dt,dz) to (dT,dZ); inverse maps (dT,dZ) to (dt,dz).
half = Q(1, 2)
J = [[(p + q)*half, (q - p)*half],
     [(q - p)*half, (p + q)*half]]
Ji = [[(p + q)*p**-1*q**-1*half, (p - q)*p**-1*q**-1*half],
      [(p - q)*p**-1*q**-1*half, (p + q)*p**-1*q**-1*half]]
equal('Jacobian determinant', J[0][0]*J[1][1] - J[0][1]*J[1][0], p*q)
for order, left, right in [('J_Jinv', J, Ji), ('Jinv_J', Ji, J)]:
    for i in range(2):
        for j in range(2):
            equal(order + '_%d%d' % (i, j),
                  sum(left[i][k]*right[k][j] for k in range(2)), int(i == j))

# Metric coefficients after division by exp(2B).
for i in range(2):
    for j in range(2):
        expected = (-1 if i == 0 else 1)*p**-1*q**-1 if i == j else 0
        equal('metric_pullback_%d%d' % (i, j),
              -Ji[0][i]*Ji[0][j] + Ji[1][i]*Ji[1][j], expected)

# Same F on both null coordinates implies U=V=F(t) on z=0.
equal('shell_Z_equals_zero', (Ft - Ft)*half)
equal('shell_T_equals_Ft', (Ft + Ft)*half, Ft)
shell = {'p': 'L', 'q': 'L', 'r': 'w', 's': 'w'}
for i in range(2):
    for j in range(2):
        equal('shell_inverse_J_%d%d' % (i, j),
              Ji[i][j].rename(shell), L**-1 if i == j else 0)

def dT(t_derivative, z_derivative):
    return Ji[0][0]*t_derivative + Ji[1][0]*z_derivative

def dZ(t_derivative, z_derivative):
    return Ji[0][1]*t_derivative + Ji[1][1]*z_derivative

# A and phi are scalar functions. B has the logarithmic conformal shift.
new_A_T, new_A_Z = dT(At, Az).rename(shell), dZ(At, Az).rename(shell)
new_phi_T, new_phi_Z = dT(ft, fz).rename(shell), dZ(ft, fz).rename(shell)
new_B_t = Bt - (r+s)*half
new_B_z = Bz + (r-s)*half
new_B_T = dT(new_B_t, new_B_z).rename(shell)
new_B_Z = dZ(new_B_t, new_B_z).rename(shell)
for name, got, want in [('A_T', new_A_T, At*L**-1),
                      ('A_Z', new_A_Z, Az*L**-1),
                      ('phi_T', new_phi_T, ft*L**-1),
                      ('phi_Z', new_phi_Z, fz*L**-1),
                      ('B_T', new_B_T, (Bt-w)*L**-1),
                      ('B_Z', new_B_Z, Bz*L**-1)]:
    equal('shell_' + name, got, want)

# E denotes exp(-B_b); exp(-Bnew_b)=L E on the shell for any F'>0.
equal('proper_time_one_form', (E**-1*L**-1)*L, E**-1)
equal('induced_metric_TT_pullback', -(E**-1*L**-1)**2,
      -E**-2*L**-2)
equal('induced_Hubble_invariance', L*E*new_A_T, E*At)
equal('scalar_proper_velocity_invariance', L*E*new_phi_T, E*ft)
new_A_TT = (Att - w*At)*L**-2
new_phi_TT = (ftt - w*ft)*L**-2
equal('Hubble_proper_derivative_invariance',
      (L*E)**2*(new_A_TT - new_B_T*new_A_T), E**2*(Att-Bt*At))
equal('scalar_proper_acceleration_invariance',
      (L*E)**2*(new_phi_TT - new_B_T*new_phi_T), E**2*(ftt-Bt*ft))
equal('normal_vector_t_component', L*E*Ji[0][1].rename(shell), 0)
equal('normal_vector_z_component', L*E*Ji[1][1].rename(shell), E)
equal('normal_scalar_derivative_invariance', L*E*new_phi_Z, E*fz)
equal('spatial_extrinsic_curvature_invariance', L*E*new_A_Z, E*Az)
equal('time_extrinsic_curvature_invariance', L*E*new_B_Z, E*Bz)
equal('covariant_extrinsic_TT_pullback', -E**-1*L**-1*new_B_Z,
      -E**-1*Bz*L**-2)
equal('extrinsic_trace_invariance', L*E*(new_B_Z+3*new_A_Z), E*(Bz+3*Az))

# Registered normal n=exp(-B) d_z; pure-tension junction residuals.
equal('normal_A_junction_residual', L*E*new_A_Z - sig*Q(1, 6),
      E*Az - sig*Q(1, 6))
equal('normal_B_junction_residual', L*E*new_B_Z - sig*Q(1, 6),
      E*Bz - sig*Q(1, 6))
equal('normal_phi_junction_residual', L*E*new_phi_Z + sigp*half,
      E*fz + sigp*half)
equal('coordinate_A_junction_residual', new_A_Z-sig*E**-1*L**-1*Q(1, 6),
      L**-1*(Az-sig*E**-1*Q(1, 6)))
equal('coordinate_B_junction_residual', new_B_Z-sig*E**-1*L**-1*Q(1, 6),
      L**-1*(Bz-sig*E**-1*Q(1, 6)))
equal('coordinate_phi_junction_residual', new_phi_Z+sigp*E**-1*L**-1*half,
      L**-1*(fz+sigp*E**-1*half))

# The clock condition adds L=exp(B_b)=E^-1 and w=dlogL/dt=B_t.
equal('clock_B_T_vanishes', new_B_T.rename({'w': 'Bt'}))
equal('clock_lapse_equals_one',
      (E**-1*L**-1).substitute_monomials({'E': L**-1}), 1)
equal('clock_lapse_direct', L*L**-1, 1)
equal('clock_H_is_A_T',
      (L*E*new_A_T - new_A_T).substitute_monomials({'E': L**-1}))

# Guard against invalid stronger statements: an arbitrary map does not fix B_T,
# and coordinate normal derivatives alone are not invariant.
nonzero('arbitrary_map_does_not_set_B_T_zero', new_B_T)
nonzero('coordinate_phi_Z_is_not_invariant', new_phi_Z - fz)

# Mass dimensions are separately verified with exact rational arithmetic.
dim_kappa5 = Q(-3, 2)
dim_registered_phi = Q(0)
dim_canonical_Phi = dim_registered_phi - dim_kappa5
dim_chi = Q(1)
dim_gb = Q(-1, 2)
dimensions = {
    'kappa5': dim_kappa5,
    'registered_phi': dim_registered_phi,
    'canonical_Phi': dim_canonical_Phi,
    'brane_chi': dim_chi,
    'brane_gb': dim_gb,
    'bulk_scalar_kinetic_density': 2*(dim_canonical_Phi+1),
    'brane_interaction_density': 2*dim_gb+2*dim_canonical_Phi+2*dim_chi,
    'production_q': dim_gb+dim_canonical_Phi+1,
}
for key, want in [('canonical_Phi', Q(3, 2)), ('bulk_scalar_kinetic_density', Q(5)),
                  ('brane_interaction_density', Q(4)), ('production_q', Q(2))]:
    if dimensions[key] != want:
        raise RuntimeError('Mass dimension failed: ' + key)
    checks.append({'name': 'dimension_' + key, 'result': 'EXACT_RATIONAL',
                   'mass_power': str(want)})

out = {
    'status': 'PASS_EXACT_LOCAL_COORDINATE_IDENTITIES',
    'check_count': len(checks),
    'method': 'Sparse Laurent polynomials with fractions.Fraction; exact zero coefficients',
    'reference_convention': {
        'file': 'Chat13/derive_evolution_equations.py',
        'sha256': '87082ca89e06f63b1cc588df9a054388b1b1dcbb22d11f00d51ef42a4d52a6fc',
        'normal': 'exp(-B) partial_z; bulk z<0',
        'junctions': 'n(A)=n(B)=sigma/6; n(phi)=-sigma_phi/2',
        'read_only_source': True,
        'reference_required_for_replay': False,
    },
    'assumptions': ['Common C3 map F for both null coordinates',
                    'Fprime positive throughout the local mapped neighborhood',
                    'Metric and scalar at least as smooth as required by the junction expressions',
                    'Fixed original brane z=0; orientation n=exp(-B) partial_z',
                    'Pure-tension junction convention inherited from Chat13',
                    'Clock specialization Fprime(t)=exp(B_b(t)) while lapse is positive'],
    'not_proved': ['Global coordinate coverage or late-time extension',
                   'Well-posedness or stability of a new PDE/gauge implementation',
                   'Numerical convergence, constraint control, or endpoint fate',
                   'Physical units of the registered numerical length scale'],
    'dimensions': {key: str(value) for key, value in dimensions.items()},
    'checks': checks,
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
path = Path(__file__).with_name('PROPER_CLOCK_EXACT_CHECKS.json')
path.write_text(json.dumps(out, indent=2, sort_keys=True) + '\n')
print(json.dumps({'status': out['status'], 'check_count': len(checks),
                  'receipt': path.name}, sort_keys=True))
