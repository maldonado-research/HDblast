#!/usr/bin/env python3
"""Cheap deliberate controls for the constraint-control laboratory.

1. Operator identity: lab order-4 operators == frozen solver operators.
2. RHS / constraint identity with the frozen solver on real initial data.
3. Stencil order on a smooth test function with a prescribed shell slope
   (orders 2, 4, 6 including the Hermite ghost closure) and a wrong-ghost control.
4. Exact SymPy check of the damped constraint-propagation law, with
   wrong-sign and wrong-denominator negative controls.
5. Projection control: the projected data change by a truncation-size amount
   that decreases with h, and the zero-amplitude reference is left exactly zero.
Writes CONTROLS.json.
"""
import json, hashlib, math, time
from pathlib import Path
import numpy as np
import sympy as sp
import lab

HERE = Path(__file__).resolve().parent
out = dict(status='PASS', checks=[])
def check(name, ok, **info):
    out['checks'].append(dict(name=name, pass_=bool(ok), **info))
    if not ok: out['status'] = 'FAIL'
    print(name, ok, info, flush=True)

# ---------------------------------------------------------------- 1,2
s, st, meta = lab.initial_data(.01, 4e-4, 3., 2.)
L4 = lab.Lab(s, order=4)
dd = [abs(L4.D1 - s.D1).max() / abs(s.D1).max(), abs(L4.D2 - s.D2).max() / abs(s.D2).max(),
      abs(L4.g1 - s.g1).max() / abs(s.g1).max(), abs(L4.g2 - s.g2).max() / abs(s.g2).max()]
check('order-4 operators equal frozen solver to roundoff (relative < 1e-14)', max(dd) < 1e-14, relative_max_differences=[float(x) for x in dd])
zt, zpt, zppt = lab.R.grid(4e-4, 3., 2.)
fa, fb = lab.R.matrices(zpt, zppt), lab.frozen_matrices_lowmem(zpt, zppt)
dm = [abs(fa[i] - fb[i]).max() / abs(fa[i]).max() for i in (0, 1)] + [abs(fa[i] - fb[i]).max() / abs(fa[i]).max() for i in (2, 3)] + [abs(fa[4] - fb[4]).max()]
check('low-memory operator builder equals frozen R.matrices (relative < 1e-14)', max(dm) < 1e-14, relative_max_differences=[float(x) for x in dm])
rng = np.random.default_rng(7)
pert = st + 1e-4 * rng.standard_normal(st.shape) * (np.arange(s.n) >= 2)
for nm, v in [('initial', st), ('perturbed', pert)]:
    r1, r2 = L4.rhs(v), s.rhs(v)
    rel = float(abs(r1 - r2).max() / abs(r2).max())
    check(f'rhs equals frozen solver to roundoff ({nm} state, relative < 1e-13)', rel < 1e-13, relative_max_difference=rel)
_, Hs, Ms, _ = s.diagnostics(pert, 0.)
Hl, Ml, _ = L4.constraints(pert)
relc = float(max(abs(Hs - Hl).max() / abs(Hs).max(), abs(Ms - Ml).max() / abs(Ms).max()))
check('constraint arrays equal frozen diagnostics to roundoff (relative < 1e-13)', relc < 1e-13, relative_max_difference=relc)
s.ko = L4.ko = 0.3; Lk = lab.Lab(s, order=4, ko=0.3)
relk = float(abs(Lk.rhs(pert) - s.rhs(pert)).max() / abs(s.rhs(pert)).max()); s.ko = 0.
check('KO dissipation equals the frozen solver KO term to roundoff', relk < 1e-13, relative_max_difference=relk)

# ---------------------------------------------------------------- 3
def stencil_errors(order, h, stretch=2., Lz=1., ghost_scale=1.0):
    z, zp, zpp = R.grid(h, Lz, stretch)
    D1, D2, g1, g2, gh = lab.build_operators(zp, zpp, order)
    u = np.sin(3 * z) + np.exp(z); ux = 3 * np.cos(3 * z) + np.exp(z); uxx = -9 * np.sin(3 * z) + np.exp(z)
    slope = ux[-1] * ghost_scale
    e1 = D1 @ u + g1 * slope - ux; e2 = D2 @ u + g2 * slope - uxx
    r = order // 2; sl = slice(r + 1, None)          # outer zero-ghost rows excluded
    return float(abs(e1[sl]).max()), float(abs(e2[sl]).max()), float(abs(e2[-1]))
R = lab.R
orders = {}
for order in (2, 4, 6):
    hs = [0.04, 0.02, 0.01] if order == 6 else [0.02, 0.01, 0.005]
    errs = [stencil_errors(order, h) for h in hs]
    p1 = [math.log2(errs[i][0] / errs[i + 1][0]) for i in range(2)]
    p2 = [math.log2(errs[i][1] / errs[i + 1][1]) for i in range(2)]
    orders[order] = dict(h=hs, errors=errs, observed_D1_order=p1, observed_D2_order=p2)
    if order < 6:
        check(f'order-{order} stencil incl. shell ghost converges at order >= {order}-0.3',
              min(p1 + p2) > order - 0.3, observed_D1=p1, observed_D2=p2)
    else:
        check('order-6 first derivative converges at order >= 5.7', min(p1) > 5.7, observed_D1=p1)
        out['finding_order6_shell_D2'] = dict(observed_D2_order=p2, errors=errs, finding=(
            'the degree-7 Hermite shell closure reaches ~1e-9 error already at h=0.01 and then '
            'loses accuracy to roundoff amplification; it is not a usable 6th-order closure at the '
            'evolution spacings (h<=4e-4). Recorded as a negative finding, not a pass/fail control.'))
e6 = [stencil_errors(6, h) for h in (0.01, 0.005, 0.0025)]
out['order6_shell_D2_roundoff_floor'] = dict(h=[0.01, 0.005, 0.0025], shell_D2_error=[e[2] for e in e6],
    note='degree-7 Hermite ghost weights amplify roundoff ~eps*W/h^2; error rises as h decreases below ~0.01')
errs = [stencil_errors(4, h, ghost_scale=1.01) for h in [0.02, 0.01, 0.005]]
pw = [math.log2(errs[i][2] / errs[i + 1][2]) for i in range(2)]
check('negative control: 1% wrong shell slope destroys boundary D2 convergence', max(pw) < 0, shell_D2_error_log2_ratios=pw)

# ---------------------------------------------------------------- 4
t, z, kap = sp.symbols('t z kappa', real=True)
A, B, p = [sp.Function(q)(t, z) for q in ('A', 'B', 'phi')]
U = sp.Function('U')(p)
At, Az, Bt, Bz, pt, pz = [sp.diff(f, x) for f in (A, B, p) for x in (t, z)]
EA = sp.diff(A, z, 2) - 3 * At ** 2 + 3 * Az ** 2 + sp.Rational(2, 3) * sp.exp(2 * B) * U
EB = sp.diff(B, z, 2) + 3 * At ** 2 - 3 * Az ** 2 - pt ** 2 / 2 + pz ** 2 / 2 - sp.exp(2 * B) * U / 3
EP = sp.diff(p, z, 2) - 3 * At * pt + 3 * Az * pz - sp.exp(2 * B) * sp.diff(U, p)
H = -2 * sp.exp(2 * B) * U + 6 * At ** 2 + 6 * At * Bt - 12 * Az ** 2 + 6 * Az * Bz - 6 * sp.diff(A, z, 2) - pt ** 2 - pz ** 2
M = -3 * sp.diff(A, t, z) - 3 * At * Az + 3 * At * Bz + 3 * Az * Bt - pt * pz
Cp, Cm = H + 2 * M, H - 2 * M
def law(damp_B):
    ev = {sp.diff(A, t, 2, z): sp.diff(EA, z), sp.diff(A, t, 2): EA, sp.diff(B, t, 2): EB + damp_B, sp.diff(p, t, 2): EP}
    red = lambda e: sp.simplify(sp.expand(e.subs(ev, simultaneous=True).doit()))
    lp = red(sp.diff(sp.exp(3 * A) * Cp, t) - sp.diff(sp.exp(3 * A) * Cp, z) + kap * sp.exp(3 * A) * Cp)
    lm = red(sp.diff(sp.exp(3 * A) * Cm, t) + sp.diff(sp.exp(3 * A) * Cm, z)
             + kap * sp.exp(3 * A) * Cp * (At - Az) / (At + Az))
    return lp, lm
t0 = time.time()
lp, lm = law(-kap * Cp / (6 * (At + Az)))
check('damped law: (d_t-d_z)[e^3A C+] = -kappa e^3A C+', lp == 0)
check('damped law: (d_t+d_z)[e^3A C-] = -kappa e^3A C+ (A_t-A_z)/(A_t+A_z)', lm == 0)
lp2, _ = law(+kap * Cp / (6 * (At + Az)))
check('negative control: wrong damping sign violates the damped law', lp2 != 0)
lp3, _ = law(-kap * Cp / (6 * (At - Az)))
check('negative control: denominator A_t-A_z violates the damped law', lp3 != 0)
out['sympy_seconds'] = time.time() - t0

# ---------------------------------------------------------------- 5
proj = []
for h in [4e-4, 2e-4]:
    s2, st2, _ = lab.initial_data(.01, h, 3., 2.)
    Lh = lab.Lab(s2)
    r0, _, _ = Lh.diagnostics(st2, 0.)
    v, log = lab.project_initial_hamiltonian(Lh, st2)
    r1, _, M1 = Lh.diagnostics(v, 0.)
    proj.append(dict(h=h, H_before=r0['H_max'], H_after=r1['H_max'], M_after=r1['M_max'],
                     Cw_before=r0['Cplus_w_max'], Cw_after=r1['Cplus_w_max'],
                     max_abs_change_a=float(abs(v[0] - st2[0]).max()), newton=log))
check('projection reduces discrete initial H to roundoff and keeps M exactly zero',
      all(p['H_after'] < 1e-12 and p['M_after'] == 0 for p in proj), rows=[{k: p[k] for k in p if k != 'newton'} for p in proj])
check('projection change decreases under refinement (converges to the continuum data)',
      proj[1]['max_abs_change_a'] < proj[0]['max_abs_change_a'] / 4,
      ratio=proj[0]['max_abs_change_a'] / proj[1]['max_abs_change_a'])
s0, st0, _ = lab.initial_data(0., 4e-4, 3., 2.)
L0 = lab.Lab(s0); v0, _ = lab.project_initial_hamiltonian(L0, st0)
check('zero-amplitude control: projection leaves the zero perturbation exactly zero', not np.any(v0))

out['stencil_orders'] = orders
out['projection_rows'] = proj
out['sha256'] = {f: hashlib.sha256((HERE / f).read_bytes()).hexdigest() for f in ['lab.py', 'controls.py']}
(HERE / 'CONTROLS.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
print(out['status'])
