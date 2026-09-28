#!/usr/bin/env python3
"""Deliberate controls for the chi mode solver.  Writes CONTROLS.json.

1. Calibration: linear crossing in flat space (h=0, a=1) must reproduce n_k = exp(-pi (k^2+mu0^2)/q).
2. Wrong-formula controls: the same numerics compared with exp(-pi k^2/(2q)) and exp(-2 pi k^2/q) must fail visibly.
3. Numerical convergence on the archived trajectory: adiabatic-window tolerance, ODE tolerance, momentum nodes,
   adiabatic order of the vacuum/particle basis, spline order.
4. Perturbed crossing point phi* +- 0.01; nonzero bare mass mu0 on the actual trajectory.
5. Continuation: occupation evaluated at the window end and at the data end must agree.
6. Asymptotic approach to the instant formula: relative number deviation times G should approach a constant.
"""
import json, math, sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ph_lib as L

PRIMARY = 'v2_t1e3_plus_c4_finecoarse'


def summary(case):
    return {'q': case['q'], 'N_comoving_numeric': case['analytic']['N_comoving_numeric'],
            'N_exact_formula': case['analytic']['N_exact_formula'],
            'rel_diff_number': case['analytic']['rel_diff_number'], 'max_abs_diff_nk': case['analytic']['max_abs_diff_nk'],
            'window': case['window'], 'wronskian_dev': case['wronskian_dev'], 'order_used': case['order_used']}


def main():
    t0 = time.time()
    out = {'status': 'CONTROLS', 'archive_sha256': L.sha256_file(L.ARCHIVE)}
    # 1-2 calibration in flat space
    cal = {}
    for v, G in ((0.5, 20.0), (0.5, 200.0), (0.58, 2000.0)):
        q = G*v
        for mu0 in (0.0, 0.5*math.sqrt(q)):
            tr = L.LinearToy(v, span=max(20.0, 3*math.sqrt(1/(G*v))*40))
            case = L.run_case(tr, G, 0.0, mu0=mu0, tol=1e-3)
            k = np.array(case['kappa']); n = np.array(case['n_k']); w = np.array(case['weights'])
            ana = np.exp(-math.pi*(k**2 + mu0**2)/q)
            good = ana > 1e-6
            wrong1 = np.exp(-math.pi*(k**2 + mu0**2)/(2*q)); wrong2 = np.exp(-2*math.pi*(k**2 + mu0**2)/q)
            def relN(x): return float(np.sum(w*k**2*n)/np.sum(w*k**2*x) - 1)
            cal['q=%g|mu0=%.4g' % (q, mu0)] = {
                'q': q, 'mu0': mu0, 'window': case['window'],
                'max_rel_diff_nk_where_exact_n>1e-6': float(np.max(np.abs(n[good]/ana[good] - 1))),
                'max_abs_diff_nk_all_nodes': float(np.max(np.abs(n - ana))),
                'rel_diff_number_correct_formula': relN(ana),
                'rel_diff_number_wrong_half_q': relN(wrong1), 'rel_diff_number_wrong_double_exponent': relN(wrong2),
                'N_numeric': float(np.sum(w*k**2*n))/(2*math.pi**2),
                'N_exact': q**1.5/(8*math.pi**3)*math.exp(-math.pi*mu0**2/q),
                'wronskian_dev': case['wronskian_dev']}
            print('cal', q, mu0, cal['q=%g|mu0=%.4g' % (q, mu0)]['max_rel_diff_nk_where_exact_n>1e-6'], flush=True)
    out['calibration_flat_linear'] = cal
    out['calibration_pass'] = all(c['max_rel_diff_nk_where_exact_n>1e-6'] < 1e-5 and abs(c['rel_diff_number_wrong_half_q']) > 0.5
                                  and abs(c['rel_diff_number_wrong_double_exponent']) > 0.5 for c in cal.values())

    # 3 convergence on the archived trajectory
    tr = L.Trajectory(PRIMARY); tr3 = L.Trajectory(PRIMARY, k=3)
    conv = {}
    for ps in (0.5, 0.75):
        for G in (100.0, 1e3, 1e4):
            base = L.run_case(tr, G, ps)
            variants = {'tol3e-4': dict(tol=3e-4), 'rtol1e-12': dict(rtol=1e-12), 'nk80': dict(nk=80),
                        'order0': dict(order=0), 'kfac60': dict(kfac=60.0)}
            row = {'base': summary(base)}
            Nb = base['analytic']['N_comoving_numeric']
            for name, kw in variants.items():
                c = L.run_case(tr, G, ps, **kw)
                row[name] = {'N_comoving_numeric': c['analytic']['N_comoving_numeric'],
                             'rel_change_vs_base': c['analytic']['N_comoving_numeric']/Nb - 1, 'window': c['window'],
                             'order_used': c['order_used']}
            c = L.run_case(tr3, G, ps)
            row['cubic_spline'] = {'N_comoving_numeric': c['analytic']['N_comoving_numeric'],
                                   'rel_change_vs_base': c['analytic']['N_comoving_numeric']/Nb - 1}
            # continuation: particle number re-evaluated at the data end
            if base['window']['s_b'] < tr.s_max - 1e-6:
                kaps = np.array(base['kappa'])
                r2 = L.evolve_modes(tr, kaps, G, ps, 0.0, base['window']['s_a'], tr.s_max, base['_lna_star'], order=2)
                N2 = float(np.sum(np.array(base['weights'])*kaps**2*r2['n']))/(2*math.pi**2)
                row['continuation_to_data_end'] = {'N_comoving': N2, 'rel_change_vs_base': N2/Nb - 1}
            conv['%s|%g' % (ps, G)] = row
            print('conv', ps, G, {k: v.get('rel_change_vs_base') for k, v in row.items() if k != 'base'}, flush=True)
    out['convergence'] = conv

    # 4 perturbed crossing point and bare mass
    pert = {}
    for ps in (0.49, 0.5, 0.51):
        c = L.run_case(tr, 1e3, ps)
        pert[str(ps)] = summary(c)
    q = pert['0.5']['q']; mu0 = 0.5*math.sqrt(q)
    c = L.run_case(tr, 1e3, 0.5, mu0=mu0)
    pert['mu0'] = {'mu0': mu0, 'N_comoving_numeric': c['analytic']['N_comoving_numeric'],
                   'ratio_to_mu0_zero': c['analytic']['N_comoving_numeric']/pert['0.5']['N_comoving_numeric'],
                   'instant_prediction_ratio': math.exp(-math.pi*mu0**2/q), 'rel_diff_number': c['analytic']['rel_diff_number']}
    out['perturbed'] = pert

    # 6 asymptotic scaling of the deviation from the instant formula
    asym = {}
    for ps in (0.5, 0.75):
        xs, ys = [], []
        for G in (1e3, 3e3, 1e4, 3e4, 1e5, 1e6):
            c = L.run_case(tr, G, ps)
            xs.append(G); ys.append(c['analytic']['rel_diff_number'])
        asym[str(ps)] = {'G': xs, 'rel_diff_number': ys, 'G_times_rel_diff': [x*y for x, y in zip(xs, ys)],
                         'loglog_slope': float(np.polyfit(np.log(xs), np.log(np.abs(ys)), 1)[0])}
    out['asymptotic_scaling'] = asym
    out['trajectory'] = {'lna_crosscheck_max_abs': tr.lna_crosscheck, 'spline_node_residual': tr.phi_fit_residual,
                         'max_rel_stored_vs_spline_velocity_on_crossing_interval':
                             float(np.max(np.abs(tr.bg(tr.s[(tr.s > 4.5) & (tr.s < 6.8)], 1)[:, 0]/tr.v_stored[(tr.s > 4.5) & (tr.s < 6.8)] - 1)))}
    out['runtime_s'] = round(time.time() - t0, 1)
    (HERE/'CONTROLS.json').write_text(json.dumps(out, indent=1, allow_nan=False, default=float) + '\n')
    print('done', out['runtime_s'], 'calibration_pass', out['calibration_pass'])


if __name__ == '__main__':
    main()
