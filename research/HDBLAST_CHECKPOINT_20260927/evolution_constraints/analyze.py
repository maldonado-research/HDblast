#!/usr/bin/env python3
"""Convergence analysis of the constraint-control runs (reads runs/, writes ANALYSIS.json).

For every variant with several resolutions: Hamiltonian/momentum residual norms at
t=0.25,0.5,0.75,1 and observed orders p=log2(E_h/E_{h/2}); pointwise self-differences
of H and of the six fields on exactly nested (shell-anchored) nodes; and the
constraint-transport gain G(t)=e^{-3t}[rho(z0)/rho(z0-t)]^3 of the background.
"""
import json, math, glob, hashlib
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
VARIANTS = {
    'baseline_order4': 'runs/final/base_h{h}',
    'remedy_A_discrete_projection': 'runs/final/proj_h{h}',
    'remedy_B_outgoing_damping_k10': 'runs/final/k10_h{h}',
    'remedy_A_plus_B': 'runs/final/k10p_h{h}',
    'remedy_C_order6_projected': 'runs/final/o6p_h{h}',
}
HS = ['4e-4', '2e-4', '1e-4', '5e-5']
TIMES = ['0.25', '0.50', '0.75', '1.00']

def nest(zc, zf, arr):
    """Values of a fine-grid array on the coarse nodes (grids share the shell node
    and every second fine node from the shell is a coarse node)."""
    rf = arr[..., ::-1][..., ::2]; zr = zf[::-1][::2]; m = min(len(zc), len(zr))
    if abs(zr[:m] - zc[::-1][:m]).max() > 1e-12: raise RuntimeError('grids not nested')
    out = np.full(arr.shape[:-1] + (len(zc),), np.nan); out[..., ::-1][..., :m] = rf[..., :m]
    return out

def restrict(zc, zf, arr, stride):
    """Fine-grid values on the coarse nodes; both grids share the shell node and
    every stride-th fine node counted from the shell is a coarse node."""
    rf = arr[..., ::-1][..., ::stride]; zr = zf[::-1][::stride]; m = min(len(zc), len(zr))
    if abs(zr[:m] - zc[::-1][:m]).max() > 1e-12: raise RuntimeError('grids not nested')
    out = np.full(arr.shape[:-1] + (len(zc),), np.nan); out[..., ::-1][..., :m] = rf[..., :m]
    return out

def l2(z, u, mask):
    zz, uu = z[mask], u[..., mask]
    return np.sqrt(np.trapezoid(uu * uu, zz, axis=-1))

result = dict(status='FLOATING_POINT_CONVERGENCE_ANALYSIS', variants={})
for name, pat in VARIANTS.items():
    runs = {}
    for h in HS:
        f = HERE / (pat.format(h=h) + '.json')
        if f.exists():
            d = json.loads(f.read_text())
            if d['stop_reason'] != 'tf': continue
            runs[h] = d
    if len(runs) < 2: continue
    hs = [h for h in HS if h in runs]
    V = dict(resolutions=hs, settings=runs[hs[0]]['settings'], nodes=[runs[h]['nodes'] for h in hs],
             runtime_seconds=[runs[h]['runtime_seconds'] for h in hs], times={})
    V['initial'] = {h: {k: runs[h]['initial'][k] for k in ('H_max', 'Cplus_w_max')} for h in hs}
    npz = {h: np.load(HERE / (pat.format(h=h) + '.npz')) for h in hs}
    for t in TIMES:
        T = {}
        for q in ('H_max', 'H_L2', 'M_max', 'M_L2', 'core_H', 'Cplus_w_max', 'Cminus_w_max', 'H_peak_z'):
            vals = [runs[h]['snapshots'][str(float(t))][q] for h in hs]
            T[q] = vals
            if q not in ('H_peak_z',):
                T[q + '_order'] = [math.log2(vals[i] / vals[i + 1]) if vals[i + 1] > 0 else None for i in range(len(vals) - 1)]
        # residual restricted to the complement of the taper's causal future
        cm = []
        for h in hs:
            A = npz[h]; z = A['z']; bg = A['background']; S0 = 1 + 6 * bg[2] ** 2 + bg[3] ** 2
            w = (z > max(-2.4, -.85 * 3. + float(t) + .05)) & (np.arange(len(z)) > 5)
            Hn = np.abs(A['H_' + t]) / S0
            cm.append((float(Hn[w].max()), float(np.sqrt(np.trapezoid(Hn[w] ** 2, z[w])))))
        T['H_max_causal'] = [c[0] for c in cm]; T['H_L2_causal'] = [c[1] for c in cm]
        for q in ('H_max_causal', 'H_L2_causal'):
            T[q + '_order'] = [math.log2(T[q][i] / T[q][i + 1]) for i in range(len(hs) - 1)]
        # self-convergence of H/S0 and of fields on the coarsest grid's nodes
        zc = npz[hs[0]]['z']; bgc = npz[hs[0]]['background']
        sc = 1 + 6 * bgc[2] ** 2 + bgc[3] ** 2
        mask = (zc > -2.4) & (np.arange(len(zc)) > 5)  # original mask; far-end NaN padding lies outside it
        Hs, Fs = [], []
        for h in hs:
            A = npz[h]; lev = hs.index(h)
            Hn = A['H_' + t]; St = A['state_' + t]; z = A['z']
            Hn = restrict(zc, z, Hn, 2 ** lev); St = restrict(zc, z, St, 2 ** lev)
            Hs.append(Hn / sc); Fs.append(St)
        dH = [float(l2(zc, Hs[i] - Hs[i + 1], mask)) for i in range(len(hs) - 1)]
        dF = [l2(zc, Fs[i] - Fs[i + 1], mask).tolist() for i in range(len(hs) - 1)]
        T['H_self_difference_L2'] = dH
        T['H_self_difference_order'] = [math.log2(dH[i] / dH[i + 1]) for i in range(len(dH) - 1)]
        T['field_self_difference_L2'] = dF
        T['field_self_difference_order'] = [[math.log2(a / b) if b > 0 else None for a, b in zip(dF[i], dF[i + 1])]
                                            for i in range(len(dF) - 1)]
        V['times'][t] = T
    result['variants'][name] = V

# constraint-transport gain along outgoing rays of the static background
f = HERE / 'runs/final/base_h2e-4.npz'
if f.exists():
    A = np.load(f); z = A['z']; rho = A['background'][0]
    gain = {}
    for z0 in (0., -0.003, -0.01):
        g = {}
        for t in (0.1, 0.25, 0.5, 0.75, 1.0):
            r0 = np.interp(z0, z, rho); r1 = np.interp(z0 - t, z, rho)
            g[str(t)] = float(math.exp(-3 * t) * (r0 / r1) ** 3)
        gain[str(z0)] = g
    result['transport_gain_unweighted'] = gain
    result['transport_gain_definition'] = 'e^{-3t}[rho(z0)/rho(z0-t)]^3 on the zero-bump reference; linear interpolation of the sampled background'

# early outgoing weighted constraint density (after the near-shell passage), from
# the t=0.05 record of every completed run (final and short diagnostic runs)
early = {}
for f in sorted(glob.glob(str(HERE / 'runs/final/*.json')) + glob.glob(str(HERE / 'runs/short/*.json'))):
    d = json.loads(Path(f).read_text())
    rec = [r for r in d['records'] if abs(r['time'] - 0.05) < 1e-3]   # nearest record (h=4e-4 has dt=1/6252)
    if not rec: continue
    st = d['settings']
    early[Path(f).parent.name + '/' + Path(f).stem] = dict(hmin=st['hmin'], order=st['order'], kappa=st['kappa'],
        ko=st['ko'], xtol=st['xtol'], project=st.get('project', False), cfl=st.get('cfl', .4),
        initial_Cplus_w_max=d['initial']['Cplus_w_max'], Cplus_w_max_t005=rec[0]['Cplus_w_max'],
        H_max_t005=rec[0]['H_max'], record_time=rec[0]['time'])
result['early_outgoing_density_t0.05'] = early

# projection correction size (projected minus sampled initial warp) per resolution
pc = {}
for h in HS:
    fa, fb = HERE / f'runs/final/proj_h{h}.npz', HERE / f'runs/final/base_h{h}.npz'
    if fa.exists() and fb.exists():
        A, B = np.load(fa), np.load(fb)
        pc[h] = dict(max_abs_delta_a=float(np.abs(A['initial'][0] - B['initial'][0]).max()),
                     delta_a_at_shell=float(A['initial'][0, -1] - B['initial'][0, -1]),
                     f_unchanged=bool(np.array_equal(A['initial'][2], B['initial'][2])))
ks = [h for h in HS if h in pc]
result['projection_correction'] = dict(rows=pc, orders=[math.log2(pc[ks[i]]['max_abs_delta_a'] / pc[ks[i + 1]]['max_abs_delta_a']) for i in range(len(ks) - 1)])

def calibration(folder):
    cal = {}
    for h in HS:
        f = HERE / f'{folder}/g_h{h}.npz'
        if not f.exists() or not (HERE / f'{folder}/g_h{h}.json').exists(): continue
        A = np.load(f); z = A['z']; bg = A['background']; S0 = 1 + 6 * bg[2] ** 2 + bg[3] ** 2
        row = {}
        for t in TIMES:
            w = (z > -.85 * 3. + float(t) + .05) & (z < -.05)
            Hn = np.abs(A['H_' + t]) / S0
            row[t] = dict(H_max=float(Hn[w].max()), H_L2=float(np.sqrt(np.trapezoid(Hn[w] ** 2, z[w]))), z_peak=float(z[w][np.argmax(Hn[w])]))
        w0 = (z > -2.4) & (z < -.05); Hn0 = np.abs(A['H_initial']) / S0
        row['initial'] = dict(H_max=float(Hn0[w0].max()), H_L2=float(np.sqrt(np.trapezoid(Hn0[w0] ** 2, z[w0]))))
        cal[h] = row
    kc = [h for h in HS if h in cal]
    calo = {}
    for t in TIMES + ['initial']:
        calo[t] = dict(H_max=[cal[h][t]['H_max'] for h in kc], H_L2=[cal[h][t]['H_L2'] for h in kc])
        for q in ('H_max', 'H_L2'):
            v = calo[t][q]; calo[t][q + '_order'] = [math.log2(v[i] / v[i + 1]) for i in range(len(v) - 1)]
    calf = {}
    if len(kc) >= 3:
        Z = {h: np.load(HERE / f'{folder}/g_h{h}.npz') for h in kc}
        zc = Z[kc[0]]['z']
        for t in TIMES:
            w = (zc > -.85 * 3. + float(t) + .05) & (zc < -.05)
            F = [restrict(zc, Z[h]['z'], Z[h]['state_' + t], 2 ** i) for i, h in enumerate(kc)]
            dF = [l2(zc, F[i] - F[i + 1], w) for i in range(len(kc) - 1)]
            calf[t] = dict(field_self_difference_L2=[d.tolist() for d in dF],
                           field_self_difference_order=[[math.log2(x / y) for x, y in zip(dF[i], dF[i + 1])] for i in range(len(dF) - 1)])
    res_fields = calf
    return res_fields, dict(resolutions=kc, by_time=calo,
        window='z > -0.85L + t + 0.05 and z < -0.05 (outside the taper causal future and away from the shell); L=3')



for key, folder, w in [('smooth_bulk_calibration', 'runs/calib3', 0.02), ('smooth_bulk_calibration_wide_pulse', 'runs/calib2', 0.05)]:
    fields, table = calibration(folder)
    table['pulse_width'] = w
    result[key] = table; result[key + '_fields'] = fields

# damping controls at h=4e-4 (exploratory runs, same code path for these options)
dc = {}
for tag, f in [('baseline', 'runs/final/base_h4e-4'), ('kappa=+10 uniform (no shell switch-off)', 'runs/damp/o4k10_h4e-4'),
               ('kappa=-5 uniform (anti-damping control)', 'runs/ctrl/o4km5_h4e-4'),
               ('naive H-only damping kappa=10 (expected unstable)', 'runs/ctrl/o4honly10_h4e-4'),
               ('kappa=+10 with shell switch-off', 'runs/final/k10_h4e-4'), ('kappa=+30 with shell switch-off', 'runs/damp/k30off_h4e-4')]:
    fj = HERE / (f + '.json')
    if not fj.exists(): continue
    d = json.loads(fj.read_text())
    dc[tag] = dict(file=f, stop_reason=d['stop_reason'], last_time=d['records'][-1]['time'],
                   snapshots={k: dict(H_max=r['H_max'], core_H=r['core_H'], H_peak_z=r['H_peak_z']) for k, r in d['snapshots'].items()})
b = dc.get('baseline'); pos = dc.get('kappa=+10 uniform (no shell switch-off)'); neg = dc.get('kappa=-5 uniform (anti-damping control)')
if b and pos and neg:
    dc['predicted_vs_measured_t0.25'] = dict(
        kappa_plus10=dict(measured_ratio=pos['snapshots']['0.25']['H_max'] / b['snapshots']['0.25']['H_max'], exp_minus_kappa_t=math.exp(-2.5)),
        kappa_minus5=dict(measured_ratio=neg['snapshots']['0.25']['H_max'] / b['snapshots']['0.25']['H_max'], exp_minus_kappa_t=math.exp(1.25)),
        kappa_minus5_t0p5=dict(measured_ratio=neg['snapshots']['0.5']['H_max'] / b['snapshots']['0.5']['H_max'], exp_minus_kappa_t=math.exp(2.5)))
result['damping_controls_h4e-4'] = dc

# uncompensated single-pulse data (runs/calib/gauss_*): corner-incompatible control
un = {}
for h in HS:
    f = HERE / f'runs/calib/gauss_h{h}.json'
    if f.exists():
        d = json.loads(f.read_text())
        un[h] = {k: dict(H_max=r['H_max'], H_peak_z=r['H_peak_z']) for k, r in d['snapshots'].items()}
result['uncompensated_single_pulse_control'] = dict(rows=un, description=(
    'epsilon=0 reference + one Gaussian scalar pulse at z=-0.3 (width 0.05, amplitude 0.01), warp from the '
    'discrete projection anchored at the far taper. The single pulse changes the shell data (no compensation), '
    'so the residual concentrates on the corner characteristic z=-t and grows under refinement.'))

# reproduction of the checkpoint wrapper by the lab code (same h=2e-4 wide run)
fo = HERE / 'runs/baseline_original/balanced_epsp01_wide_h2.json'
if fo.exists():
    o = json.loads(fo.read_text()); lb = json.loads((HERE / 'runs/final/base_h2e-4.json').read_text())
    ro = [r for r in o['records'] if abs(r['time'] - 1) < 1e-9][0]
    result['wrapper_reproduction'] = dict(original_copy_H_t1=ro['hamiltonian_max'], lab_H_t1=lb['snapshots']['1.0']['H_max'],
        archived_checkpoint_H_t1=0.30047427375135166,
        relative_difference_lab_vs_copy=abs(ro['hamiltonian_max'] - lb['snapshots']['1.0']['H_max']) / ro['hamiltonian_max'])

result['sha256'] = {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in ['lab.py', 'run_case.py', 'analyze.py', 'calib.py', 'controls.py']}
(HERE / 'ANALYSIS.json').write_text(json.dumps(result, indent=1) + '\n')
for name, V in result['variants'].items():
    print(name, V['resolutions'])
    for t, T in V['times'].items():
        print('  t', t, 'H_max', ['%.3e' % x for x in T['H_max']], 'ord', ['%.2f' % x for x in T['H_max_order']],
              '| causal', ['%.2e' % x for x in T['H_max_causal']], ['%.2f' % x for x in T['H_max_causal_order']],
              '| L2c ord', ['%.2f' % x for x in T['H_L2_causal_order']], '| selfdiff ord', ['%.2f' % x for x in T['H_self_difference_order']])
print(json.dumps({k: result[k] for k in ('projection_correction', 'smooth_bulk_calibration', 'smooth_bulk_calibration_fields', 'wrapper_reproduction') if k in result}, indent=1))
