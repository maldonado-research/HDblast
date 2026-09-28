#!/usr/bin/env python3
"""Run one balanced-seed evolution with a chosen constraint-control option set.
Writes <out>.json (records, settings, hashes) and <out>.npz (snapshots)."""
import argparse, hashlib, json, time, os
from pathlib import Path
import numpy as np
import lab

p = argparse.ArgumentParser()
p.add_argument('--epsilon', type=float, default=.01)
p.add_argument('--hmin', type=float, default=2e-4)
p.add_argument('--L', type=float, default=3.)
p.add_argument('--stretch', type=float, default=2.)
p.add_argument('--order', type=int, default=4)
p.add_argument('--kappa', type=float, default=0.)
p.add_argument('--damp-form', default='cplus')
p.add_argument('--ko', type=float, default=0.)
p.add_argument('--damp-off', type=float, nargs=2, default=None, help='kappa switched off for z>-d0, full for z<-d1')
p.add_argument('--xtol', type=float, default=2e-14)
p.add_argument('--tf', type=float, default=1.)
p.add_argument('--cfl', type=float, default=.4)
p.add_argument('--gauss', type=float, nargs=3, default=None, help='smooth calibration data: add AMP*exp(-((z-Z0)/W)^2) to f (use with --epsilon 0 --project)')
p.add_argument('--project', action='store_true', help='discrete Hamiltonian projection of initial data')
p.add_argument('--out', required=True)
a = p.parse_args()
out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
t0 = time.monotonic()
s, state, meta = lab.initial_data(a.epsilon, a.hmin, a.L, a.stretch, a.xtol)
if a.tf >= .85 * a.L: raise ValueError('tf reaches taper causal window')
L = lab.Lab(s, order=a.order, kappa=a.kappa, ko=a.ko, damp_form=a.damp_form, damp_off=a.damp_off)
proj_log = None
if a.gauss:
    amp, z0, w = a.gauss
    state = state.copy(); state[2] += amp * np.exp(-((s.z - z0) / w) ** 2); state[:, :2] = 0
if a.project:
    raw_row, _, _ = L.diagnostics(state, 0.)
    state, proj_log = lab.project_initial_hamiltonian(L, state)
    proj_log = dict(before=raw_row, iterations=proj_log)
t_setup = time.monotonic() - t0
init_row, H0, M0 = L.diagnostics(state, 0.)
print(json.dumps({'setup_seconds': t_setup, 'nodes': L.n, 'initial': init_row}), flush=True)
res, snaps, final = lab.evolve(L, state, tf=a.tf, cfl=a.cfl, progress=lambda r: print(json.dumps(r), flush=True))
here = Path(lab.__file__).resolve().parent
hashes = {f: hashlib.sha256((here / f).read_bytes()).hexdigest() for f in
          ['lab.py', 'run_case.py', 'src/registered_solver.py', 'src/evolve_balanced.py', 'src/balanced_constraint_seed.py', 'src/constraint_seed.py']}
doc = dict(status='FLOATING_POINT_CONSTRAINT_CONTROL_RUN', settings=vars(a), nodes=L.n, dt=res['dt'], steps=res['steps'],
           stop_reason=res['reason'], initial=init_row, projection=proj_log, snapshots={str(k): v['row'] for k, v in snaps.items()},
           records=res['rows'], seed_meta=dict(analytic_boundary_slope_minus_target=meta['analytic_boundary_slope_minus_target'],
           initial_shell_scalar_displacement=meta['initial_shell_scalar_displacement']),
           runtime_seconds=time.monotonic() - t0, sha256=hashes,
           scope='Classical matter-free floating-point run; no interval certificate; diagnostics defined in lab.py.')
out.with_suffix('.json').write_text(json.dumps(doc, indent=1) + '\n')
arrs = dict(z=L.z, background=np.array([L.rho, L.phi, L.hc, L.phz]), initial=state, H_initial=H0, M_initial=M0, final=final)
for k, v in snaps.items():
    tag = f'{k:.2f}'
    arrs['state_' + tag] = v['state']; arrs['H_' + tag] = v['H']; arrs['M_' + tag] = v['M']
np.savez_compressed(str(out) + '.npz', **arrs)
print(json.dumps({'done': str(out), 'reason': res['reason'], 'runtime': doc['runtime_seconds'],
                  'snap': doc['snapshots']}), flush=True)
