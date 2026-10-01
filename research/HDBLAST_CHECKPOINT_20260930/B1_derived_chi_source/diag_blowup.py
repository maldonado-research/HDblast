#!/usr/bin/env python3
"""B1 diagnostic (informative, not part of the registered rule): where and when does the bulk evolution become non-finite?
Reruns one evolve_b1.py configuration unchanged (same Run class, same arguments) but with dense snapshots (every snap_dT in
coordinate time T), then locates the first bulk point with phi < -1 (past the W critical point phi = -1 of the bulk potential),
the global phi minimum and its Z position per snapshot, and the shell values.  Output: diag/<tag>_blowup.json.
Usage: python3 diag_blowup.py <tag> <snap_dT> -- <evolve_b1.py arguments>"""
import sys, os, json, math, hashlib
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import evolve_b1 as E
import shell_matter as SM

def main():
    tag = sys.argv[1]; sdt = float(sys.argv[2]); args = sys.argv[sys.argv.index('--') + 1:]
    import argparse
    ap = argparse.ArgumentParser()
    for k, t, dflt in [('delta', float, 0.1), ('dc', float, 1e-2), ('G', float, 0.0), ('phistar', float, 0.5), ('y', float, 1.0),
                       ('b', float, 0.0), ('Y', float, 0.0), ('Gt', float, 10.0), ('dzf', float, 1e-3), ('dzc', float, 4e-3),
                       ('L', float, 16.0), ('kappa', float, 0.0), ('project', float, None), ('Tf', float, 14.0), ('xc', float, math.inf)]:
        ap.add_argument('--' + k, type=t, default=dflt)
    ap.add_argument('--source', default='chi'); ap.add_argument('--dstar', action='store_true'); ap.add_argument('--xc_from', default=None)
    a = ap.parse_args(args)
    M8 = json.load(open('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/M8_QUADRATIC_TENSION_TUNING.json'))
    d = M8['d_star']['0.1'] if a.dstar else 0.0
    xc = a.xc
    if a.xc_from:
        xs = json.load(open(a.xc_from))['xc_suggest']; xc = math.inf if xs is None else xs
    if a.source == 'friction': matter = SM.Friction(a.Y); mp = dict(source='friction', Y=a.Y)
    else:
        matter = (lambda rb: SM.ChiGas(a.G, a.phistar, a.y, a.b, rb)); mp = dict(source='chi', G=a.G, phistar=a.phistar, y=a.y, b=a.b)
    run = E.Run(delta=a.delta, d=d, dc=a.dc, matter=matter, matter_params=mp, dzf=a.dzf, dzc=a.dzc, L=a.L, xc=xc, kappa=a.kappa,
                project=a.project, log=lambda s: print(s, flush=True))
    out_dir = HERE/'diag'; out_dir.mkdir(exist_ok=True)
    run.evolve(T_final=a.Tf, tag=tag, out_dir=str(out_dir), snap_dT=sdt)
    z = np.load(out_dir/(tag + '_timeseries.npz')); s = json.load(open(out_dir/(tag + '_summary.json')))
    Z = z['z'][::4] if len(z['z']) != len(z['snap_0.000_phi']) else z['z']
    Ts = sorted({k.split('_')[1] for k in z.files if k.startswith('snap_')}, key=float)
    rows = []
    for T in Ts:
        ph = z['snap_%s_phi' % T]; B = z['snap_%s_B' % T]
        fin = np.isfinite(ph)
        i = int(np.nanargmin(ph)) if fin.any() else -1
        below = np.flatnonzero(ph < -1.0)
        rows.append(dict(T=float(T), phi_min=float(ph[i]) if i >= 0 else None, Z_phi_min=float(Z[i]) if i >= 0 else None,
                         Z_range_phi_below_m1=[float(Z[below].min()), float(Z[below].max())] if len(below) else None,
                         n_nonfinite=int((~fin).sum()), B_max=float(np.nanmax(B)), Z_B_max=float(Z[int(np.nanargmax(B))]),
                         phi_shell=float(ph[-1])))
    cols = list(z['cols']); rec = z['rec']
    res = dict(status='numerical diagnostic (informative)', tag=tag, args=args, stop_reason=s['stop_reason'], T_end=s['T_end'], H0tau_end=s['H0tau_end'],
               snapshots=rows, code_sha256=hashlib.sha256((HERE/'evolve_b1.py').read_bytes()).hexdigest())
    (out_dir/(tag + '_blowup.json')).write_text(json.dumps(res, indent=1) + '\n')
    for r in rows: print(r)

if __name__ == '__main__':
    main()
