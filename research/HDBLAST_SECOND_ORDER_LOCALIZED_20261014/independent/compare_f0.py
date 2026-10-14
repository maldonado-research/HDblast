"""Compare the independently measured l = 0 remainder with the original
engine's closed-form f0 (l0_weight.json) on the test background."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp, mpmath as mp
import indep_check as ic
here = os.path.dirname(os.path.abspath(__file__))
f0 = sp.sympify(json.load(open(os.path.join(here, '..', 'basis', 'results', 'l0_weight.json')))['f0'])
fn = sp.lambdify(ic.JET_SYMS, f0, 'mpmath')
d = json.load(open(os.path.join(here, 'indep_check.json')))
out = {}
for case, (G, K) in {'commuting_A_eq_B': (ic.Gf, ic.Kf), 'pure_velocity_G0': (0, ic.Kf)}.items():
    pk = ic.Packet(G, K)
    for sph in d['cases'][case]['spheres']:
        r = mp.mpf(sph['r']); sv = r * r
        eng = fn(*pk.jetvals(sv))
        meas = mp.mpf(sph['f0_from_remainder'])
        out[f'{case}_r{sph["r"]}'] = {'engine_f0': mp.nstr(eng, 20), 'independent_f0': mp.nstr(meas, 20),
                                      'rel_diff': mp.nstr(abs(eng - meas) / abs(eng), 5)}
print(json.dumps(out, indent=2))
json.dump(out, open(os.path.join(here, 'compare_f0.json'), 'w'), indent=2)
