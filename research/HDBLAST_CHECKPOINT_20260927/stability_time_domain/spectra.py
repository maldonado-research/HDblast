"""Discrete-operator spectra of the registered 1+1 PDE system linearised about a static shell.

Dense mode: full eigen-decomposition of the reduced linear operator (frozen outer nodes removed),
every eigenpair classified (constraint residual, gauge-template residual, shell observables).
Shift mode: shift-invert Arnoldi near chosen targets for finer grids.

Usage examples
  python3 spectra.py --bg plus --hmin 4e-4 --mode dense --out results/spectra/plus_h4_L6
  python3 spectra.py --bg plus --hmin 1e-4 --mode shift --shifts 0.5,0,-1,-1.45 --out results/spectra/plus_h1_shift
Control option --s20-scale multiplies sigma''(phi_b) in the scalar junction (wrong-formula control).
"""
import argparse, json, time, os
from pathlib import Path
import numpy as np, scipy.linalg as sl
from scipy.sparse.linalg import eigs
import td_core as T


def build(bg_name, hmin, L, stretch, tdet, s20_scale=1.0, precise=True):
    bg = T.plus_background(tdet) if bg_name == 'plus' else T.original_background(tdet)
    s = T.GeneralSolver(bg, tdet=tdet, hmin=hmin, L=L, stretch=stretch, precise_eta=precise)
    if s20_scale != 1.0: s.s20 = s.s20*s20_scale
    return s


def summarize_modes(s, w, vecs_full, top=40):
    rows = []
    for k in range(len(w)):
        d = T.classify(s, w[k], vecs_full[:, k]); d['label'] = T.label(d)
        lam = w[k]; d['mu2'] = [float((-lam*(lam + 3)).real), float((-lam*(lam + 3)).imag)]
        d['on_dS_line'] = bool(abs(lam.imag) > 1e-6 and abs(lam.real + 1.5) < 1e-6)
        rows.append(d)
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--bg', choices=['plus', 'original'], required=True)
    p.add_argument('--hmin', type=float, default=4e-4); p.add_argument('--L', type=float, default=6.)
    p.add_argument('--stretch', type=float, default=.05); p.add_argument('--tdet', type=float, default=1e-3)
    p.add_argument('--mode', choices=['dense', 'shift'], default='dense')
    p.add_argument('--shifts', default='0.5,0,-1,-1.45'); p.add_argument('--k', type=int, default=8)
    p.add_argument('--s20-scale', type=float, default=1.0)
    p.add_argument('--phi-poly', action='store_true', help='use the registered phi-polynomial potential derivatives (precision control)')
    p.add_argument('--out', required=True)
    a = p.parse_args(); Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    s = build(a.bg, a.hmin, a.L, a.stretch, a.tdet, a.s20_scale, not a.phi_poly)
    M, jerr = T.verify_matrix(s)
    if a.s20_scale == 1.0 and jerr > 1e-12: raise RuntimeError('matrix/RHS mismatch %g' % jerr)
    idx = T.free_index(s.n)
    out = dict(script='spectra.py', background=s.meta, bg=a.bg, hmin=a.hmin, L=a.L, stretch=a.stretch, nodes=s.n,
               tdet=a.tdet, s20_scale=a.s20_scale, matrix_vs_complex_step_rhs_rel_error=jerr,
               background_junction_residual_on_grid=s.junction_on_grid, mode=a.mode, pot_mode=s.pot_mode)
    if a.mode == 'dense':
        A = M.toarray()[np.ix_(idx, idx)]
        w, V = sl.eig(A)
        full = np.zeros((6*s.n, len(w)), complex); full[idx] = V
        # eigen-residuals
        res = np.linalg.norm(A@V - V*w, axis=0)/np.linalg.norm(V, axis=0)
        rows = summarize_modes(s, w, full)
        for r, e in zip(rows, res): r['eigen_residual'] = float(e)
        keep = [r for r in rows if r['imag'] >= -1e-9]
        keep.sort(key=lambda r: -r['real'])
        out['n_eigenvalues'] = len(w)
        out['max_real_all'] = keep[0]['real']
        out['top_by_real_part'] = keep[:40]
        # Physical-sector bookkeeping.
        line = [r for r in keep if r['on_dS_line']]
        off = [r for r in keep if not r['on_dS_line']]
        out['n_on_dS_line_upper_half'] = len(line)
        out['max_abs_real_plus_1p5_on_line'] = float(max(abs(r['real'] + 1.5) for r in line)) if line else None
        out['real_axis_eigenvalues'] = [r for r in off if abs(r['imag']) < 1e-6]
        out['off_line_complex_with_real_gt_minus3'] = [r for r in off if abs(r['imag']) >= 1e-6 and r['real'] > -3][:80]
        # Shell-supported modes above the continuum line: what would signal a bound state or instability.
        out['shell_supported_above_line'] = [r for r in keep if r['real'] > -1.5 + 1e-6 and r['shell_phi_rel'] > 1e-6]
        low = sorted([r for r in line if r['shell_phi_rel'] > 1e-3], key=lambda r: r['imag'])
        out['lowest_frequency_shell_supported_line_modes'] = low[:10]
        big = sorted(line, key=lambda r: -r['shell_phi_rel'])
        out['largest_shell_amplitude_line_modes'] = big[:10]
        np.savez_compressed(a.out + '_eigs.npz', w=w, eigen_residual=res)
    else:
        found = []
        Mr = M.tocsr()[idx][:, idx].tocsc()
        for sh in [complex(x) for x in a.shifts.split(',')]:
            try:
                vals, vr = eigs(Mr, k=a.k, sigma=sh, tol=1e-11, maxiter=5000)
            except Exception as ex:
                found.append(dict(shift=[sh.real, sh.imag], error=repr(ex))); continue
            modes = []
            for j in np.argsort(np.abs(vals - sh)):
                v = np.zeros(6*s.n, complex); v[idx] = vr[:, j]
                d = T.classify(s, vals[j], v); d['label'] = T.label(d)
                d['eigen_residual'] = float(np.linalg.norm(Mr@vr[:, j] - vals[j]*vr[:, j])/np.linalg.norm(vr[:, j]))
                d['mu2'] = [float((-vals[j]*(vals[j] + 3)).real), float((-vals[j]*(vals[j] + 3)).imag)]
                modes.append(d)
            found.append(dict(shift=[sh.real, sh.imag], modes=modes))
        out['shift_invert'] = found
    out['runtime_seconds'] = time.time() - t0
    T.dump(a.out + '.json', out)
    print(json.dumps(dict(out=a.out, nodes=s.n, runtime=out['runtime_seconds'], jerr=jerr)), flush=True)


if __name__ == '__main__':
    main()
