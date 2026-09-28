"""Step 1: reproduce the static backgrounds from the (read-only) package scripts.

+1 branch  : static_branch/solve_plus_branch.py :: solve(delta)          (22 Sept 2026 checkpoint)
original   : frozen/registered_solver.py :: shell_background(tdet)        (unstable shell, cone near phi=-1)
Bytecode writing is disabled so that nothing is created inside the package folders.
Each background is then re-integrated and root-polished by gi_core.Background (independent code path).
Output: BACKGROUNDS.json
"""
import sys, json, time, hashlib, math
sys.dont_write_bytecode = True
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
PKG = Path('/home/user/unified-theory-maldonado/new-files/latest-work/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922')
sys.path.insert(0, str(PKG/'static_branch')); sys.path.insert(0, str(PKG/'frozen')); sys.path.insert(0, str(HERE))
import solve_plus_branch as spb
import registered_solver as rs
from gi_core import Background, C_REG

STATED = dict(phi_b=0.9999159473169134, rho_b=129.9247628497, H2=5.92401479433e-5, H_over_H0=0.606721732)
RHO0_ARCHIVE = 78.82817714224423
DELTAS = [0.0003, 0.001, 0.003, 0.01, 0.03, 0.1]
t0 = time.time()
out = dict(c=C_REG, package_c=spb.C, package_hashes={
    'solve_plus_branch.py': hashlib.sha256((PKG/'static_branch/solve_plus_branch.py').read_bytes()).hexdigest(),
    'registered_solver.py': hashlib.sha256((PKG/'frozen/registered_solver.py').read_bytes()).hexdigest()},
    plus=[], original=[])
assert spb.C == C_REG == rs.C
for d in DELTAS:
    row, _ = spb.solve(d)
    bg = Background(+1, d, e_h=row['eta_h'], y_b=row['y_b'])
    sv = bg.shell_values()
    rec = dict(delta=d, package=dict(phi_b=row['phi_b'], rho_b=row['rho_b'], H2=row['H2'], eta_h=row['eta_h'], y_b=row['y_b'],
                                     junction_residual=row['junction_residual']),
               gi_core=dict(e_h=bg.e_h, y_b=bg.y_b, polish=bg.polish_info, **{k: float(v) for k, v in sv.items()}))
    rec['rel_diff_rho_b'] = abs(sv['rho_b']/row['rho_b'] - 1)
    rec['diff_phi_b'] = abs(sv['phi_b'] - row['phi_b'])
    out['plus'].append(rec); print(json.dumps(rec), flush=True)
r1 = [r for r in out['plus'] if r['delta'] == 0.001][0]['package']
out['stated_digits_check'] = dict(
    phi_b=dict(got=r1['phi_b'], stated=STATED['phi_b'], ok=abs(r1['phi_b'] - STATED['phi_b']) < 5e-16),
    rho_b=dict(got=r1['rho_b'], stated=STATED['rho_b'], ok=abs(r1['rho_b'] - STATED['rho_b']) < 5e-11),
    H2=dict(got=r1['H2'], stated=STATED['H2'], ok=abs(r1['H2'] - STATED['H2']) < 5e-17),
    H_over_H0=dict(got=math.sqrt(r1['H2'])*RHO0_ARCHIVE, stated=STATED['H_over_H0'],
                   ok=abs(math.sqrt(r1['H2'])*RHO0_ARCHIVE - STATED['H_over_H0']) < 5e-10))
print(out['stated_digits_check'], flush=True)
for d in [0.001, 0.003, 0.01]:
    samp, meta = rs.shell_background(d)
    bg = Background(-1, d, e_h=meta['psi_h'], y_b=meta['yb'])
    sv = bg.shell_values()
    rec = dict(delta=d, package=meta, gi_core=dict(e_h=bg.e_h, y_b=bg.y_b, polish=bg.polish_info, **{k: float(v) for k, v in sv.items()}))
    rec['rel_diff_rho_b'] = abs(sv['rho_b']/meta['rho_b'] - 1)
    out['original'].append(rec); print(json.dumps(rec, default=float), flush=True)
out['original_rho0_matches_archive'] = abs(out['original'][0]['gi_core']['rho_b'] - RHO0_ARCHIVE)
out['runtime_s'] = time.time() - t0
out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(HERE/'BACKGROUNDS.json').write_text(json.dumps(out, indent=2, default=float) + '\n')
print('done', out['runtime_s'])
