#!/usr/bin/env python3
"""Analyze saved runs only; standard library, no evolution or eigenproblem replay."""
from pathlib import Path
import hashlib, json, math, statistics

ROOT = Path(__file__).resolve().parent
EXPECTED = -1.5 + math.sqrt(2.25 + 7.717871625176294)
checks = []

def check(name, condition):
    checks.append(dict(name=name, passed=bool(condition)))
    if not condition:
        raise RuntimeError(name)

def read(stem):
    p = ROOT / 'runs' / (stem + '.json')
    return json.loads(p.read_text()), hashlib.sha256(p.read_bytes()).hexdigest()

def fit(records, lo, hi):
    rows = [r for r in records if lo <= r['time'] <= hi]
    x = [r['time'] for r in rows]
    y = [math.log(abs(r['scalar_deviation'])) for r in rows]
    mx, my = statistics.mean(x), statistics.mean(y)
    slope = sum((a-mx)*(b-my) for a,b in zip(x,y))/sum((a-mx)**2 for a in x)
    return dict(window=[lo,hi], samples=len(x), slope=slope,
                max_log_residual=max(abs(b-my-slope*(a-mx)) for a,b in zip(x,y)))

spectra = []
for stem in ['spectrum_reg_h4','spectrum_reg_h2','spectrum_reg_h1','spectrum_reg_h05',
             'spectrum_reg_h1_L8','spectrum_t01_h4']:
    d, digest = read(stem)
    val = d['eigenvalues'][0]['real']
    expected = EXPECTED if d['background']['tdet'] == .001 else 1.6270213424531335
    spectra.append(dict(name=stem, sha256=digest, detuning=d['background']['tdet'],
                        hmin=d['hmin'], L=d['L'], nodes=d['nodes'], growth=val,
                        expected_growth=expected, error=val-expected,
                        eigen_residual=d['eigenvalues'][0]['absolute_eigen_residual'],
                        matrix_rhs_relative_error=d['matrix_rhs_relative_error']))
    check(stem+' matrix agrees with nonlinear RHS',d['matrix_rhs_relative_error']<2e-12)
check('registered selected rate approaches previous independent rate with refinement',
      all(abs(spectra[i+1]['error'])<abs(spectra[i]['error']) for i in range(3)))
check('finest selected rate is within 1e-7 of independent rate',abs(spectra[3]['error'])<1e-7)
domain_change = spectra[4]['growth']-spectra[2]['growth']
check('L6 versus L8 selected-rate change below 1e-8',abs(domain_change)<1e-8)
check('larger-detuning selected-rate control within 1e-7',abs(spectra[5]['error'])<1e-7)

evolutions = []
for stem in ['reg_mode_plus_h2','reg_mode_minus_h2','reg_mode_ko_h2','reg_mode_h1','reg_plus_h2','reg_plus_h1']:
    d, digest = read(stem)
    rec = d['records']
    mode = bool(d['mode_seed'])
    windows = [(0.5,1.5),(1,2.5)] if mode else [(2.5,3.5),(3,4.5)]
    row = dict(name=stem, sha256=digest, seed='linear eigenmode' if mode else 'shifted tension',
               hmin=d['hmin'], ko=d['ko'], amplitude=d['amplitude'] if mode else None,
               dc=d['dc'], t_end=d['t_end'], fits=[fit(rec,*w) for w in windows],
               max_normalized_H=max(r['hamiltonian_max'] for r in rec),
               max_normalized_M=max(r['momentum_max'] for r in rec),
               max_weighted_characteristic=max(r['weighted_characteristic_max'] for r in rec),
               initial=rec[0], final=rec[-1])
    evolutions.append(row)
    check(stem+' reached requested final time',d['stop_reason']=='tf' and d['t_end']==(2.5 if mode else 4.5))
    check(stem+' recovers expected growth within 1e-4',abs(row['fits'][-1]['slope']-EXPECTED)<1e-4)
check('fine mode H residual improves by more than factor8',evolutions[0]['max_normalized_H']>8*evolutions[3]['max_normalized_H'])
check('fine mode M residual improves by more than factor8',evolutions[0]['max_normalized_M']>8*evolutions[3]['max_normalized_M'])
ko_change=evolutions[2]['fits'][-1]['slope']-evolutions[0]['fits'][-1]['slope']
sign_change=evolutions[1]['fits'][-1]['slope']-evolutions[0]['fits'][-1]['slope']
check('small KO change in fitted growth',abs(ko_change)<1e-7)
check('opposite-sign small-mode growth agrees',abs(sign_change)<1e-6)
# A diagnostic deliberately records a failure of an intended future acceptance gate.
seed_refinement_fails=(evolutions[5]['max_normalized_H']>evolutions[4]['max_normalized_H']
                       and evolutions[5]['max_normalized_M']>evolutions[4]['max_normalized_M'])
check('negative finding: shifted-tension seed has worsening refined constraints',seed_refinement_fails)

result=dict(status='PASS_SAVED_CONTROL_ANALYSIS', expected_growth=EXPECTED, checks=checks,
            spectra=spectra, evolutions=evolutions, L6_L8_growth_change=domain_change,
            KO_growth_change=ko_change, opposite_sign_growth_change=sign_change,
            scope='Saved-run analysis only. Linear-limit controls, not full nonlinear fate, complete mode counting, interval certification, or evidence of a cosmological origin.',
            constraints='H,M are divided by 1+6Hc(z)^2+phi_z(z)^2 on z>−0.8L excluding first 6 nodes. This is a background normalization, not a relative perturbation-error bound. Raw arrays are in NPZ snapshots.',
            nonlinear_seed_gate='Neither evolved seed family is fully nonlinear compatible. The shifted-tension family also fails constraint refinement. Separate initial_data construction is not evolved here.')
(ROOT/'CONTROL_ANALYSIS.json').write_text(json.dumps(result,indent=2)+'\n')
lines=['# Completed registered-shell controls','',
       'Analysis of six spectra and six short nonlinear-PDE evolutions. All numbers are floating-point results. The expected rate was supplied before the new eigensolves; comparison is a code-validation test, not a blind physical prediction.','',
       '| Detuning | Shell spacing | L | Nodes | Selected growth | Difference from prior rate |',
       '|---:|---:|---:|---:|---:|---:|']
for r in spectra:
    lines.append(f"| {r['detuning']} | {r['hmin']} | {r['L']} | {r['nodes']} | {r['growth']:.12f} | {r['error']:.3e} |")
lines+=['','The finest selected eigenpair has absolute residual %.3e; the observed rate difference is not a rigorous error bound. Other near-shift eigenvalues are not independently validated.'%spectra[3]['eigen_residual'],'',
        '| Evolution | Spacing | Fitted growth | Fit time window | Maximum normalized H | Maximum normalized M |',
        '|---|---:|---:|---|---:|---:|']
for r in evolutions:
    f=r['fits'][-1]
    lines.append(f"| {r['name']} | {r['hmin']} | {f['slope']:.12f} | {f['window']} | {r['max_normalized_H']:.3e} | {r['max_normalized_M']:.3e} |")
lines+=['','Mode seeds use shell scalar amplitude ±10⁻⁸. They test the linear limit of the nonlinear equations, ending at coordinate time 2.5. Their full nonlinear junction mismatch starts at second order. Shifted-tension seeds use dc=10⁻⁷ and end at 4.5; their constraints worsen under refinement, so they do not pass nonlinear initial-data acceptance.',
        '',f'The L6→L8 selected-rate change is {domain_change:.3e}. The KO 0→0.02 fitted-rate change is {ko_change:.3e}; the opposite-sign fitted-rate difference is {sign_change:.3e}. These are bounded controls, not guarantees of late-time behavior.',
        '',result['constraints'],'',
        'Time and rates are dimensionless: one unit of the static coordinate time equals one initial-shell Hubble time. No observed cosmological time or energy scale has been calibrated.',
        '',f'All {len(checks)} analysis checks passed, including the check that preserves the negative shifted-seed finding. This script did not rerun the spectra or evolutions. Full fits, inputs, hashes and final rows are retained in CONTROL_ANALYSIS.json.','']
(ROOT/'CONTROL_ANALYSIS.md').write_text('\n'.join(lines))
print(json.dumps(dict(status=result['status'],checks=len(checks),finest_growth=spectra[3]['growth'],fine_evolution_growth=evolutions[3]['fits'][-1]['slope'])))
