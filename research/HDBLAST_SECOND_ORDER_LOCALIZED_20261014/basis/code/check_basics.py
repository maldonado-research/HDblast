"""Basic exact checks: background equations and first-order constraints."""
import json, sys, time
import sympy as sp
from engine import X, O, F, U, register, d, reduce_bg, is_zero, dS, bg_rules, S2
from fields import reconstruction, constraints, harmonic_quadratic

out = {}
t0 = time.time()
rules, s = bg_rules()
# 1. integrability of the background: dU/ds = U'(phi) Phi'
U0 = rules[U[0]]
U1 = rules[U[1]]
chk = sp.simplify(dS(U0, rules) - U1 * F[1])
out['background_dU_consistency'] = (chk == 0)
print('dU/ds = U1 Phi_s :', chk == 0, flush=True)

# 2. background Hamiltonian constraint on T=0
z = sp.zeros(4)
Hm, M = constraints(z, 0, z, 0)
ok, num = is_zero(Hm.a)
out['background_hamiltonian'] = ok
print('background H = 0 :', ok, num if not ok else '', flush=True)

# 3. first-order linearized constraints for the inherited reconstruction
G = register('G')
K = register('K')
a = sp.symbols('a11 a12 a13 a14 a22 a23 a24 a33 a34')
A = sp.Matrix([[a[0], a[1], a[2], a[3]],
               [a[1], a[4], a[5], a[6]],
               [a[2], a[5], a[7], a[8]],
               [a[3], a[6], a[8], -(a[0] + a[4] + a[7])]])
b = sp.symbols('b11 b12 b13 b14 b22 b23 b24 b33 b34')
B = sp.Matrix([[b[0], b[1], b[2], b[3]],
               [b[1], b[4], b[5], b[6]],
               [b[2], b[5], b[7], b[8]],
               [b[3], b[6], b[8], -(b[0] + b[4] + b[7])]])
# general radial profiles times degree-2 harmonics
q = G[0] * harmonic_quadratic(A)
p = K[0] * harmonic_quadratic(B)
v = reconstruction(q, p)
t1 = time.time()
Hm, M = constraints(v['h'], v['dphi'], v['dpi'], v['dvarpi'])
print('constraint jets built in %.1fs' % (time.time() - t1), flush=True)
ok, num = is_zero(Hm.b)
out['first_order_hamiltonian_general_quadrupole'] = ok
print('DH[v] = 0 :', ok, flush=True)
okM = []
for i in range(4):
    oki, numi = is_zero(M[i].b)
    okM.append(oki)
out['first_order_momentum_general_quadrupole'] = all(okM)
print('DM[v] = 0 :', okM, flush=True)
# also other angular degrees: l = 0, 1, 3 harmonic polynomials
for name, P in [('l0', sp.Integer(1)), ('l1', X[0] + 2 * X[2]),
                ('l3', X[0]**3 - 3 * X[0] * X[1]**2),
                ('l4', X[0]**4 - 6 * X[0]**2 * X[1]**2 + X[1]**4)]:
    v2 = reconstruction(G[0] * P, K[0] * P)
    H2, M2 = constraints(v2['h'], v2['dphi'], v2['dpi'], v2['dvarpi'])
    okh, _ = is_zero(H2.b)
    okm = all(is_zero(M2[i].b)[0] for i in range(4))
    out[f'first_order_{name}'] = okh and okm
    print(name, 'DH=0', okh, 'DM=0', okm, flush=True)
out['seconds'] = round(time.time() - t0, 1)
json.dump(out, open('../results/basics.json', 'w'), indent=2)
print(json.dumps(out, indent=2))
