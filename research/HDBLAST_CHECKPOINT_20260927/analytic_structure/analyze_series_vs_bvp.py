"""Compare exact series (SERIES_COEFFICIENTS.json) with high-precision BVP solutions (runs/BVP_SCAN_*.json).
Writes SERIES_VS_BVP.json.  Tests:
 (i)   scaled remainders r_n(delta) = (Q - S_n)/delta^{n+1} -> q_{n+1};
 (ii)  blind extraction: Richardson / polynomial extrapolation of (Q - S_2)/delta^3 on the geometric delta list
       (uses ONLY the previously known O(delta), O(delta^2) terms) -> estimates of q_3, q_4, q_5;
 (iii) wrong-coefficient control (q_3 -> q_3 (1+1e-4));  (iv) perturbed-parameter control (c -> c(1+1e-6));
 (v)   convergence: two solver settings (dps 50/order 50/h 0.125 vs dps 65/order 64/h 0.1).
"""
import json, hashlib
from pathlib import Path
import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
mp.mp.dps = 60
ser = json.loads((HERE / 'SERIES_COEFFICIENTS.json').read_text())
cs = sp.Symbol('c')

def coeffs(name, cval):
    return [mp.mpf(str(sp.N(sp.sympify(e).subs(cs, sp.Float(str(cval), 70)), 65))) for e in ser[name]]

def quantities(row, cval):
    d = mp.mpf(row['delta'])
    return dict(eta_b=mp.mpf(row['eta_b']), H2=mp.mpf(row['H2']),
                rho_b_factor=mp.mpf(row['rho_b']) * mp.sqrt(d * (1 + cval) / 27),
                eta_h_over_delta8=mp.mpf(row['eta_h']) / d**8,
                y_b_correction=mp.mpf(row['y_b']) - mp.mpf(9) / 2 * mp.log(4 / (3 * (1 + cval) * d)))

def S(q, d, n):
    return mp.fsum(q[i] * d**i for i in range(n + 1))

def fmt(x, n=12):
    return mp.nstr(x, n)

report = dict(status='NUMERICAL comparison of exact series with high-precision BVP', scans={})
names = ['eta_b', 'H2', 'rho_b_factor', 'eta_h_over_delta8', 'y_b_correction']
for scanfile in sorted((HERE / 'runs').glob('BVP_SCAN_*.json')):
    scan = json.loads(scanfile.read_text())
    cval = mp.mpf(scan['c'])
    qs = {n: coeffs(n, cval) for n in names}
    rows = scan['rows']
    entry = dict(c=scan['c'], dps=scan['dps'], taylor_order=scan['taylor_order'], h=scan['step_h'], per_delta=[])
    for row in rows:
        d = mp.mpf(row['delta']); Q = quantities(row, cval)
        rec = dict(delta=row['delta'], J2_residual=row['J2'], first_integral_residual=row['first_integral_residual'])
        for nme in names:
            q = qs[nme]; top = len(q) - 1
            rr = []
            for n in range(0, top):
                r = (Q[nme] - S(q, d, n)) / d**(n + 1)
                target = q[n + 1]
                rr.append(dict(n=n, scaled_remainder=fmt(r, 15), exact_next=fmt(target, 15),
                               rel_dev=fmt(abs(r - target) / abs(target), 6) if target != 0 else fmt(abs(r), 6)))
            last = (Q[nme] - S(q, d, top)) / d**(top + 1)
            rec[nme] = dict(value=fmt(Q[nme], 40), remainders=rr, remainder_after_all_known_terms_over_delta_power=fmt(last, 8),
                            abs_diff_from_full_series=fmt(abs(Q[nme] - S(q, d, top)), 6))
        entry['per_delta'].append(rec)
    # (ii) blind extraction on geometric list
    geo = [r for r in rows if mp.mpf(r['delta']) in [mp.mpf('0.0001') * 2**k for k in range(7)]]
    if len(geo) == 7 and cval != 0:
        blind = {}
        for nme in ['eta_b', 'H2']:
            q = qs[nme]
            xs = [mp.mpf(r['delta']) for r in geo]
            ys = [(quantities(r, cval)[nme] - S(q, mp.mpf(r['delta']), 2)) / mp.mpf(r['delta'])**3 for r in geo]
            V = mp.matrix([[x**j for j in range(7)] for x in xs]); coef = mp.lu_solve(V, mp.matrix(ys))
            # lower-degree fit (5 points) as a stability check
            V5 = mp.matrix([[x**j for j in range(5)] for x in xs[:5]]); coef5 = mp.lu_solve(V5, mp.matrix(ys[:5]))
            blind[nme] = {f'q{3 + j}': dict(extrapolated_7pt=fmt(coef[j], 20), extrapolated_5pt=fmt(coef5[j], 20), exact=fmt(q[3 + j], 20),
                                            rel_err_7pt=fmt(abs(coef[j] - q[3 + j]) / abs(q[3 + j]), 4)) for j in range(3)}
        entry['blind_extraction'] = blind
        # (iii) wrong-coefficient control and (iv) perturbed-parameter control at smallest delta
        r0 = geo[0]; d0 = mp.mpf(r0['delta']); Q0 = quantities(r0, cval)
        q = list(qs['H2']); q_wrong = list(q); q_wrong[3] = q[3] * (1 + mp.mpf('1e-4'))
        good = (Q0['H2'] - S(q, d0, 3)) / d0**4; bad = (Q0['H2'] - S(q_wrong, d0, 3)) / d0**4
        qp = coeffs('H2', cval * (1 + mp.mpf('1e-6')))
        pert = (Q0['H2'] - S(qp, d0, 3)) / d0**4
        entry['controls'] = dict(delta=r0['delta'], exact_q4=fmt(q[4], 15), correct_scaled_remainder=fmt(good, 15),
                                 wrong_q3_times_1p1e4_scaled_remainder=fmt(bad, 15), perturbed_c_1e6_scaled_remainder=fmt(pert, 15),
                                 wrong_rel_dev=fmt(abs(bad - q[4]) / abs(q[4]), 4), correct_rel_dev=fmt(abs(good - q[4]) / abs(q[4]), 4),
                                 perturbed_rel_dev=fmt(abs(pert - q[4]) / abs(q[4]), 4))
    report['scans'][scanfile.stem] = entry

# (v) convergence: compare reg and reg_conv at common deltas
if 'BVP_SCAN_reg' in report['scans'] and (HERE / 'runs' / 'BVP_SCAN_reg_conv.json').exists():
    A = {r['delta']: r for r in json.loads((HERE / 'runs' / 'BVP_SCAN_reg.json').read_text())['rows']}
    B = {r['delta']: r for r in json.loads((HERE / 'runs' / 'BVP_SCAN_reg_conv.json').read_text())['rows']}
    conv = []
    for dlt in sorted(set(A) & set(B), key=lambda s: mp.mpf(s)):
        conv.append({k: fmt(abs(mp.mpf(A[dlt][k]) - mp.mpf(B[dlt][k])) / abs(mp.mpf(B[dlt][k])), 4) for k in ['eta_b', 'H2', 'rho_b', 'eta_h', 'y_b']} | dict(delta=dlt))
    report['convergence_reg_vs_reg_conv_relative_differences'] = conv
report['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(HERE / 'SERIES_VS_BVP.json').write_text(json.dumps(report, indent=2) + '\n')

# compact console summary
for k, e in report['scans'].items():
    print('==', k, 'c=', e['c'][:12])
    for rec in e['per_delta']:
        print(' delta', rec['delta'], ' H2 rel_dev by order:', [x['rel_dev'] for x in rec['H2']['remainders']],
              ' eta_b:', [x['rel_dev'] for x in rec['eta_b']['remainders']])
    if 'blind_extraction' in e:
        print(json.dumps(e['blind_extraction'], indent=0)); print(json.dumps(e['controls'], indent=0))
if 'convergence_reg_vs_reg_conv_relative_differences' in report:
    print(json.dumps(report['convergence_reg_vs_reg_conv_relative_differences'], indent=0))
