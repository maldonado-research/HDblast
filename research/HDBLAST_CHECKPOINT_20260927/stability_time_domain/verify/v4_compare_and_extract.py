"""Audit V4: (a) compare the audit re-runs in verify/rerun/ with the audited results/;
(b) re-extract rates from the audited raw time series with independent code (ESPRIT pole
estimator, plain log-slope fits), recompute the Richardson numbers including the h=4e-4 point,
recompute line deviation and RK4 fully-discrete rates from the saved eigenvalue files;
(c) inspect the initial data actually used (effective weight of the 'shell-localised' bump).
Writes verify/v4_compare_and_extract.json.  No PDE solves."""
import json, glob
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent; RES = ROOT/'results'; RR = HERE/'rerun'
GN = 1.6571936312453648
out = {}


def esprit(y, dt, order, frac=0.4):
    N = len(y); Lp = int(N*frac)
    H = np.array([y[i:i + Lp] for i in range(N - Lp + 1)]).T
    U, S, Vh = np.linalg.svd(H, full_matrices=False)
    Us = U[:, :order]
    Phi = np.linalg.lstsq(Us[:-1], Us[1:], rcond=None)[0]
    z = np.linalg.eigvals(Phi)
    lam = np.log(z.astype(complex))/dt
    Z = np.vander(z, N, increasing=True).T
    amp = np.linalg.lstsq(Z, y.astype(complex), rcond=None)[0]
    o = np.argsort(-np.abs(amp))
    return [(float(lam[k].real), float(lam[k].imag), float(abs(amp[k]))) for k in o]


def slope(t, y, a, b):
    m = (t >= a) & (t <= b); return float(np.polyfit(t[m], np.log(y[m]), 1)[0])


# (a) re-run comparison -----------------------------------------------------------------
cmp = {}
for f in sorted(RR.glob('*.json')):
    name = f.stem; orig = RES/('spectra' if (name.startswith('plus_h') or 'shift' in name) else 'td')/(name + '.json')
    if not orig.exists(): cmp[name] = 'no original'; continue
    a = json.loads(f.read_text()); b = json.loads(orig.read_text())
    ent = {}
    if 'shift_invert' in a:
        def modes(d):
            return sorted({(round(m['real'], 7), round(m['imag'], 5)) for blk in d['shift_invert'] for m in blk.get('modes', [])
                           if m['eigen_residual'] < 1e-6 and m['imag'] > -1e-9}, key=lambda x: -x[0])
        ma, mb = modes(a), modes(b)
        ent['n_modes'] = [len(ma), len(mb)]
        common = [(x, min(mb, key=lambda y: abs(complex(*x) - complex(*y)))) for x in ma]
        ent['max_abs_diff_matched'] = float(max(abs(complex(*x) - complex(*y)) for x, y in common))
        ent['top_modes_rerun'] = ma[:6]
    elif a.get('mode') == 'dense':
        wa = np.load(str(f).replace('.json', '_eigs.npz'))['w']; wb = np.load(str(orig).replace('.json', '_eigs.npz'))['w']
        wa = wa[np.argsort(-wa.real)][:50]; wb = wb[np.argsort(-wb.real)][:50]
        ent['top50_max_abs_diff'] = float(max(min(abs(x - wb)) for x in wa))
        ent['max_real'] = [a['max_real_all'], b['max_real_all']]
        ent['max_dev_line'] = [a['max_abs_real_plus_1p5_on_line'], b['max_abs_real_plus_1p5_on_line']]
    else:
        za = np.load(str(f).replace('.json', '.npz')); zb = np.load(str(orig).replace('.json', '.npz'))
        ent['f_b_max_abs_diff_over_max'] = float(np.max(np.abs(za['f_b'] - zb['f_b']))/np.max(np.abs(zb['f_b'])))
        ent['E_max_rel_diff'] = float(np.max(np.abs(za['E']/zb['E'] - 1)))
        ent['beta'] = [a['initial_data']['beta'], b['initial_data']['beta']]
    cmp[name] = ent
out['rerun_vs_original'] = cmp

# (b) independent extraction from the audited raw time series -----------------------------
rates = {}
for f in sorted((RES/'td').glob('*.json')):
    d = json.loads(f.read_text()); z = np.load(str(f).replace('.json', '.npz')); t = z['t']; name = f.stem
    A = d['args']; te = t[-1]; dt = t[1] - t[0]
    r = dict(hmin=A['hmin'], bg=A['bg'], s20=A['s20_scale'])
    if A['bg'] == 'original' or A['s20_scale'] < 0:
        tp = 2.0 if A['bg'] == 'original' else 0.3*te
        m = t >= tp; y = z['f_b'][m]
        # decimate to keep the Hankel matrix small
        step = max(1, len(y)//1500); poles = esprit(y[::step], dt*step, 6)
        real = [p for p in poles if abs(p[1]) < 1e-6]
        r['esprit_dominant_real_pole'] = real[0][0] if real else None
        r['log_slope_last_1p5'] = slope(t, np.abs(z['f_b']), te - 1.5, te)
        if A['bg'] == 'original' and real: r['esprit_minus_GN'] = real[0][0] - GN
    else:
        r['E_half_slope_1_to_tf'] = slope(t, z['E'], 1.0, te)/2
        r['E_half_slope_3_to_tf'] = slope(t, z['E'], 3.0, te)/2
        r['E_near_half_slope_1_to_tf'] = slope(t, z['E_near'], 1.0, te)/2
        fb = z['f_b']; r['f_b_abs_max'] = float(np.max(np.abs(fb)))
        r['f_b_abs_max_t_gt_1'] = float(np.max(np.abs(fb[t > 1])))
        r['f_b_abs_max_t_gt_5'] = float(np.max(np.abs(fb[t > 5])))
        r['compensated_e1p5t_fb_max_t_gt_1'] = float(np.max(np.abs(fb[t > 1]*np.exp(1.5*t[t > 1]))))
    rates[name] = r
out['independent_extraction'] = rates
seq = [rates['plus_lin_h%s' % h]['E_half_slope_1_to_tf'] for h in ['4e-4', '2e-4', '1e-4', '5e-5']]
d1, d2, d3 = seq[1] - seq[0], seq[2] - seq[1], seq[3] - seq[2]
out['plus_rate_sequence_h_4_2_1_0p5e-4'] = seq
out['sequence_monotone_toward_-1.5'] = bool(abs(seq[0] + 1.5) > abs(seq[1] + 1.5) > abs(seq[2] + 1.5) > abs(seq[3] + 1.5))
p = np.log2(abs(d2/d3)); out['richardson_last_three'] = dict(order=float(p), extrap=float(seq[3] + d3/(2**p - 1)))
out['richardson_assuming_order_2'] = float(seq[3] + d3/3)
out['richardson_assuming_order_4'] = float(seq[3] + d3/15)

# dense spectra: line deviation, RK4 fully-discrete rate, positive semi-discrete modes
sp_ = {}
R4 = lambda q: 1 + q + q*q/2 + q**3/6 + q**4/24
for f in sorted((RES/'spectra').glob('*_eigs.npz')):
    name = f.name.replace('_eigs.npz', ''); d = json.loads((RES/'spectra'/(name + '.json')).read_text())
    w = np.load(f)['w']; dt = 0.4*d['hmin']
    rate = np.log(np.abs(R4(w*dt)))/dt
    up = w[w.imag > 1e-6]; line = up[np.abs(up.real + 1.5) < 1e-6]
    pos = w[(w.real > 1e-6)]
    sp_[name] = dict(n=len(w), rk4_max_rate=float(rate.max()), rk4_argmax=[float(w[rate.argmax()].real), float(w[rate.argmax()].imag)],
                     n_upper_on_line=int(len(line)), max_dev_on_line=float(np.max(np.abs(line.real + 1.5))) if len(line) else None,
                     n_positive_real_part=int(len(pos)),
                     positive_real_part_modes_abs_imag_gt_1000=sorted([[float(x.real), float(x.imag)] for x in pos if abs(x.imag) > 1000 and x.imag > 0], key=lambda q: -q[0])[:5],
                     positive_real_part_modes_abs_imag_lt_1000=sorted([[float(x.real), float(x.imag)] for x in pos if abs(x.imag) <= 1000 and x.imag >= 0], key=lambda q: -q[0])[:8])
out['dense_spectra_recomputed'] = sp_

# (c) initial data composition
ini = {}
for f in sorted((RES/'td').glob('*.json')):
    d = json.loads(f.read_text()); A = d['args']; beta = d['initial_data']['beta']
    ini[f.stem] = dict(bump1=[A['z1'], A['w1']], bump2=[A['z2'], A['w2']], beta=beta,
                       relative_weight_of_bump1_after_normalisation=float(1/abs(beta)) if abs(beta) > 1 else 1.0,
                       support_of_dominant_bump=[A['z2'] - A['w2'], A['z2'] + A['w2']] if abs(beta) > 1 else [A['z1'] - A['w1'], A['z1'] + A['w1']])
out['initial_data_composition'] = ini
(HERE/'v4_compare_and_extract.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(dict(cmp=cmp, seq=seq, mono=out['sequence_monotone_toward_-1.5'], rich=out['richardson_last_three'],
                      r2=out['richardson_assuming_order_2'], r4=out['richardson_assuming_order_4']), indent=1)[:6000])
for k, v in rates.items(): print(k, {kk: (round(vv, 8) if isinstance(vv, float) else vv) for kk, vv in v.items()})
for k, v in sp_.items(): print(k, v['rk4_max_rate'], v['max_dev_on_line'], v['n_positive_real_part'], v['positive_real_part_modes_abs_imag_gt_1000'][:2])
