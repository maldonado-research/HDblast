#!/usr/bin/env python3
"""Independent recomputation of the workstream's tables from its saved runs (runs/**/*.json, *.npz).

Own code (does not import analyze.py). Compares against README values (typed in below from
README.md) and ANALYSIS.json, and adds audit diagnostics the workstream did not report:
  * field self-difference orders restricted to the complement of the taper's causal future
    and to a near-shell window, plus the location of the largest self-difference;
  * the S0 (background normalisation) factor between the shell and the front;
  * a grid-scale (odd-even) decomposition of the initial discrete Hamiltonian defect near the shell.
Writes recompute_saved.json.
"""
import json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
HS = ['4e-4', '2e-4', '1e-4', '5e-5']
VAR = dict(base='baseline_order4', proj='remedy_A_discrete_projection', k10='remedy_B_outgoing_damping_k10',
           k10p='remedy_A_plus_B', o6p='remedy_C_order6_projected')
README = {  # H_max table (README section 2), as printed
    'base': {'0.5': [2.45e-1, 4.49e-2, 7.55e-3, 6.38e-3], '1.0': [1.33, 3.01e-1, 3.36e-2, 2.78e-2]},
    'proj': {'0.5': [2.27e-1, 4.83e-2, 7.01e-3, 9.35e-4], '1.0': [1.25, 2.70e-1, 3.03e-2, 3.46e-3]},
    'k10': {'0.5': [5.96e-3, 1.16e-3, 1.78e-4, 6.50e-5], '1.0': [4.34e-2, 1.15e-2, 1.19e-3, 1.43e-4]},
    'k10p': {'0.5': [5.39e-3, 1.17e-3, 1.81e-4, 2.28e-5], '1.0': [4.34e-2, 1.15e-2, 1.19e-3, 1.44e-4]},
    'o6p': {'0.5': [6.88e-2, 9.73e-3, 1.21e-3, 8.41e-5], '1.0': [2.42e-1, 3.59e-2, 4.54e-3, 3.56e-4]}}
README_ORD = {'base': {'0.5': [2.45, 2.57, .24], '1.0': [2.14, 3.16, .27]},
              'proj': {'0.5': [2.23, 2.78, 2.91], '1.0': [2.22, 3.15, 3.13]},
              'k10': {'0.5': [2.36, 2.70, 1.45], '1.0': [1.92, 3.27, 3.06]},
              'k10p': {'0.5': [2.20, 2.69, 2.99], '1.0': [1.92, 3.27, 3.05]},
              'o6p': {'0.5': [2.82, 3.01, 3.85], '1.0': [2.75, 2.99, 3.67]}}
README_L2C = {'base': [.58, .63], 'proj': [3.48, 3.38], 'k10': [1.71, 3.10], 'k10p': [3.40, 3.10], 'o6p': [3.97, 3.91]}

def sig_match(x, y, digits=3):
    """True if y (printed, 3 significant digits) equals x rounded to 3 significant digits (+-1 in last digit)."""
    if y == 0: return x == 0
    e = math.floor(math.log10(abs(y))) - (digits - 1)
    return abs(x - y) <= 1.01 * 10 ** e

def restrict(zc, zf, arr, stride):
    rf = arr[..., ::-1][..., ::stride]; zr = zf[::-1][::stride]; m = min(len(zc), len(zr))
    assert np.abs(zr[:m] - zc[::-1][:m]).max() < 1e-12
    out = np.full(arr.shape[:-1] + (len(zc),), np.nan); out[..., ::-1][..., :m] = rf[..., :m]
    return out

def l2(z, u):
    # own trapezoid rule
    return float(np.sqrt(np.sum(0.5 * (u[..., 1:] ** 2 + u[..., :-1] ** 2) * np.diff(z), axis=-1))) if u.ndim == 1 else \
        np.sqrt(np.sum(0.5 * (u[..., 1:] ** 2 + u[..., :-1] ** 2) * np.diff(z), axis=-1))

out = dict(status='AUDIT_RECOMPUTATION', table_H_max={}, mismatches=[])
AN = json.loads((ROOT / 'ANALYSIS.json').read_text())
for v, name in VAR.items():
    rows = {}
    for t in ('0.5', '1.0'):
        vals = [json.loads((ROOT / f'runs/final/{v}_h{h}.json').read_text())['snapshots'][t]['H_max'] for h in HS]
        ords = [math.log2(vals[i] / vals[i + 1]) for i in range(3)]
        ok_v = [sig_match(x, y) for x, y in zip(vals, README[v][t])]
        ok_o = [abs(round(o, 2) - r) <= 0.011 for o, r in zip(ords, README_ORD[v][t])]
        rows[t] = dict(values=vals, orders=ords, readme_values_match=ok_v, readme_orders_match=ok_o)
        if not (all(ok_v) and all(ok_o)): out['mismatches'].append(dict(variant=v, t=t, values=vals, orders=ords))
    out['table_H_max'][v] = rows

# causal-window L2 orders and field self-differences with three windows
L = 3.
out['L2_causal'] = {}; out['field_selfdiff'] = {}
for v, name in VAR.items():
    npz = {h: np.load(ROOT / f'runs/final/{v}_h{h}.npz') for h in HS}
    l2o = []
    fs = {}
    for t, tk in (('0.50', 0.5), ('1.00', 1.0)):
        e = []
        for h in HS:
            A = npz[h]; z = A['z']; bg = A['background']; S0 = 1 + 6 * bg[2] ** 2 + bg[3] ** 2
            w = (z > -.85 * L + tk + .05) & (np.arange(len(z)) > 5)
            e.append(l2(z[w], np.abs(A['H_' + t])[w] / S0[w]))
        l2o.append(math.log2(e[2] / e[3]))
        zc = npz[HS[0]]['z']
        F = [restrict(zc, npz[h]['z'], npz[h]['state_' + t], 2 ** i) for i, h in enumerate(HS)]
        windows = {'original_mask_z>-2.4': (zc > -2.4) & (np.arange(len(zc)) > 5),
                   'causal_z>-0.85L+t+0.05': zc > -.85 * L + tk + .05,
                   'bulk_-1.5<z<-0.2': (zc > -1.5) & (zc < -.2),
                   'near_shell_z>-0.2': zc > -.2,
                   'taper_future_only_-2.4<z<-0.85L+t+0.05': (zc > -2.4) & (zc < -.85 * L + tk + .05)}
        W = {}
        for wn, m in windows.items():
            d = [l2(zc[m], (F[i] - F[i + 1])[:, m]) for i in range(3)]
            W[wn] = dict(L2=[x.tolist() for x in d],
                         order_last_pair=[math.log2(a / b) for a, b in zip(d[1], d[2])],
                         order_first_pairs=[math.log2(a / b) for a, b in zip(d[0], d[1])])
        m = windows['original_mask_z>-2.4']
        dpf = np.abs(F[2][5] - F[3][5]); dpf[~m] = 0
        W['pf_last_pair_argmax_z'] = float(zc[np.nanargmax(dpf)])
        W['pf_last_pair_max'] = float(np.nanmax(dpf))
        fs[t] = W
    out['L2_causal'][v] = dict(orders_last_pair_t05_t1=l2o, readme=README_L2C[v],
                               match=[abs(round(a, 2) - b) <= .011 for a, b in zip(l2o, README_L2C[v])])
    if not all(out['L2_causal'][v]['match']): out['mismatches'].append(dict(variant=v, L2_causal=l2o))
    out['field_selfdiff'][v] = fs
    # consistency with ANALYSIS.json (original mask)
    for t in ('0.50', '1.00'):
        mine = fs[t]['original_mask_z>-2.4']['order_last_pair']
        theirs = AN['variants'][name]['times'][t]['field_self_difference_order'][-1]
        if max(abs(a - b) for a, b in zip(mine, theirs)) > 1e-6:
            out['mismatches'].append(dict(variant=v, t=t, field_orders_mine=mine, field_orders_analysis=theirs))

# early outgoing density table (README section 3)
early = {}
for tag in ['final/base', 'final/proj', 'final/o6p']:
    for h in HS:
        d = json.loads((ROOT / f'runs/{tag}_h{h}.json').read_text())
        r = [x for x in d['records'] if abs(x['time'] - .05) < 1e-3][0]
        early[f'{tag}_h{h}'] = dict(initial=d['initial']['Cplus_w_max'], t005=r['Cplus_w_max'], time=r['time'])
for f in sorted((ROOT / 'runs/short').glob('*.json')):
    d = json.loads(f.read_text()); r = [x for x in d['records'] if abs(x['time'] - .05) < 1e-3][0]
    r2 = [x for x in d['records'] if abs(x['time'] - .1) < 1e-3]
    early['short/' + f.stem] = dict(initial=d['initial']['Cplus_w_max'], t005=r['Cplus_w_max'], t01=r2[0]['Cplus_w_max'] if r2 else None)
out['early_density'] = early

# S0 normalisation factor and transport-gain wording check
A = np.load(ROOT / 'runs/final/base_h5e-5.npz'); z = A['z']; bg = A['background']
rho, hc, phz = bg[0], bg[2], bg[3]; S0 = 1 + 6 * hc ** 2 + phz ** 2
def at(zz, arr): return float(np.interp(zz, z, arr))
d = json.loads((ROOT / 'runs/final/base_h5e-5.json').read_text())
zpk = d['snapshots']['0.5']['H_peak_z']
H0 = A['H_initial']; Hn0 = np.abs(H0) / S0
near = z > -0.02
out['normalisation_check'] = dict(
    S0_at_shell=float(S0[-1]), S0_at_minus0p003=at(-.003, S0), S0_at_front_t05=at(zpk, S0), front_z_t05=zpk,
    rho_ratio_cubed_e3t_gain_from_m0p003=float(math.exp(-1.5) * (at(-.003, rho) / at(-.003 - .5, rho)) ** 3),
    initial_max_H_over_S0_near_shell=float(Hn0[near].max()), initial_max_abs_H_near_shell=float(np.abs(H0[near]).max()),
    H_max_normalised_t05=d['snapshots']['0.5']['H_max'],
    implied_growth_of_normalised_H=float(d['snapshots']['0.5']['H_max'] / Hn0[near].max()),
    note='The conserved-density gain e^{-3t}[rho(z0)/rho(z0-t)]^3 applies to un-normalised H+2M; H/S0 additionally '
         'gains S0(z0)/S0(front).')

# grid-scale content of the initial defect near the shell, h=1e-4 and 5e-5 baseline
noise = {}
for h in ('1e-4', '5e-5'):
    A = np.load(ROOT / f'runs/final/base_h{h}.npz'); z = A['z']; H0 = A['H_initial']
    rho = A['background'][0]; wt = (rho / rho[-1]) ** 3
    m = (z > -0.02) & (z < -1e-3)
    u = (wt * H0)[m]
    smooth = np.convolve(u, [.25, .5, .25], mode='same')   # removes the odd-even (Nyquist) mode exactly
    rough = u - smooth
    noise[h] = dict(max_weighted_H=float(np.abs(u).max()), max_smooth_part=float(np.abs(smooth[1:-1]).max()),
                    max_nyquist_part=float(np.abs(rough[1:-1]).max()),
                    z_of_max=float(z[m][np.argmax(np.abs(u))]))
out['initial_defect_near_shell'] = noise
(HERE / 'recompute_saved.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(dict(mismatches=out['mismatches'], L2=out['L2_causal'], norm=out['normalisation_check'], noise=noise), indent=1))
for v in VAR:
    for t in ('0.50', '1.00'):
        W = out['field_selfdiff'][v][t]
        print(v, t, 'pf argmax z', round(W['pf_last_pair_argmax_z'], 4), W['pf_last_pair_max'])
        for wn in ('original_mask_z>-2.4', 'causal_z>-0.85L+t+0.05', 'bulk_-1.5<z<-0.2', 'near_shell_z>-0.2', 'taper_future_only_-2.4<z<-0.85L+t+0.05'):
            print('   %-42s' % wn, ['%.2f' % x for x in W[wn]['order_last_pair']], ['%.1e' % x for x in W[wn]['L2'][2]])
