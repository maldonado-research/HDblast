"""Where does the purely local (shell) expansion stop?  Exact probe of the first resonance.

The bulk formal series N(X,T), P(X,T) (from SERIES_COEFFICIENTS.json, total degree <= 8) is extended to degree 9.
At (m,j)=(2,7) the linear operator (p-14)(p+18), p=14m-2j, vanishes; the equation there reduces to
   0 * n_{2,7} + S_{2,7} = 0.
If S_{2,7} != 0 a secular term u X^2 T^7 = u alpha^2 e^{14u} is required, and the coefficient n_{2,7} (the
coefficient of alpha^2 e^{14u}) becomes a free constant fixed only by the global cone-regularity condition.
Also fits the numerically measured delta^9 remainder of eta_b and H^2 with and without a log(delta) term.
Writes RESONANCE_PROBE.json
"""
import json, hashlib
from fractions import Fraction as Fr
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
ser = json.loads((HERE / 'SERIES_COEFFICIENTS.json').read_text())
K = 9
N = {tuple(int(t) for t in k.split(',')): Fr(v) for k, v in ser['bulk_N'].items()}
P = {tuple(int(t) for t in k.split(',')): Fr(v) for k, v in ser['bulk_P'].items()}
def mul(a, b):
    o = {}
    for (m1, j1), v in a.items():
        for (m2, j2), w in b.items():
            if m1 + m2 + j1 + j2 <= K:
                o[(m1 + m2, j1 + j2)] = o.get((m1 + m2, j1 + j2), 0) + v * w
    return o
def add(*xs):
    o = {}
    for a in xs:
        for k, v in a.items(): o[k] = o.get(k, 0) + v
    return o
def sc(s, a): return {k: s * v for k, v in a.items()}
def D(a): return {(m, j): (14 * m - 2 * j) * v for (m, j), v in a.items()}
N2 = mul(N, N); N3 = mul(N2, N); N4 = mul(N3, N); N5 = mul(N4, N)
F = add(sc(252, N), sc(450, N2), sc(-54, N3), sc(-180, N4), sc(-36, N5))
E1 = add(mul(P, D(D(N))), sc(4, mul(add(P, D(P)), D(N))), sc(-1, mul(P, F)))
S27 = E1.get((2, 7), 0)
# lower-degree residuals must vanish (consistency of the imported series)
maxlow = max(abs(v) for (m, j), v in E1.items() if m + j <= 8)
kappa = -S27 / 32   # secular coefficient: (D^2+4D-252)(u X^2 T^7) = 32 X^2 T^7 + ... at leading order in T
out = dict(status='EXACT (rational arithmetic) resonance probe + NUMERICAL remainder fits',
           S_2_7=str(S27), S_2_7_float=float(S27), secular_coefficient_kappa=str(kappa), kappa_float=float(kappa),
           max_abs_lower_degree_residual=str(maxlow),
           interpretation=('S_2_7 != 0: the regular solution contains kappa*u*alpha^2 e^{14u}. At the shell this can be absorbed into the '
                           'locally solved growing-mode amplitude, but the constant n_2_7 (coefficient of alpha^2 e^{14u}) is fixed only by '
                           'cone regularity. Hence eta_b, H^2, rho_b are determined by the local shell analysis through delta^8; the first '
                           'globally determined coefficient enters at delta^9, and eta_h acquires a relative delta^8*log(delta) correction.'))
# numerical delta^9 remainders
mp.mp.dps = 60
import sympy as sp
cs = sp.Symbol('c')
fits = {}
for label in ['reg', 'cm04', 'cp13']:
    scan = json.loads((HERE / 'runs' / f'BVP_SCAN_{label}.json').read_text())
    cval = mp.mpf(scan['c'])
    geo = [r for r in scan['rows'] if mp.mpf(r['delta']) in [mp.mpf('0.0001') * 2**k for k in range(7)]]
    fits[label] = {}
    for name in ['eta_b', 'H2']:
        q = [mp.mpf(str(sp.N(sp.sympify(e).subs(cs, sp.Float(scan['c'], 70)), 65))) for e in ser[name]]
        xs = [mp.mpf(r['delta']) for r in geo]
        ys = [(mp.mpf(r[name]) - mp.fsum(q[i] * x**i for i in range(9))) / x**9 for r, x in zip(geo, xs)]
        # model 1: polynomial a0 + a1 d + a2 d^2 + a3 d^3 (4 smallest points after the first, which sits at the precision floor)
        pts = list(range(1, 6))
        V = mp.matrix([[xs[i]**k for k in range(5)] for i in pts]); c1 = mp.lu_solve(V, mp.matrix([ys[i] for i in pts]))
        # model 2: with log term b ln d
        V2 = mp.matrix([[1, mp.log(xs[i]), xs[i], xs[i]**2, xs[i]**3] for i in pts]); c2 = mp.lu_solve(V2, mp.matrix([ys[i] for i in pts]))
        fits[label][name] = dict(delta9_coefficient_poly_fit=mp.nstr(c1[0], 12), with_log_model_constant=mp.nstr(c2[0], 12),
                                 with_log_model_log_coefficient=mp.nstr(c2[1], 6),
                                 log_to_constant_ratio=mp.nstr(abs(c2[1] / c2[0]), 4),
                                 scaled_remainders=[mp.nstr(y, 12) for y in ys])
out['numerical_delta9_fits'] = fits
# eta_h: predicted relative delta^8 log(delta) term.  eta_h = (136/3) alpha, alpha = alpha_loc (1 - kappa alpha_loc u_b + const*alpha_loc),
# alpha_loc = X_b T_b^7, u_b = -(1/2) ln T_b  =>  eta_h/delta^8 contains  B delta^8 ln(delta),  B = (kappa/2)(3/136) e_0^2
etah = {}
for label in ['reg', 'cm04', 'cp13']:
    scan = json.loads((HERE / 'runs' / f'BVP_SCAN_{label}.json').read_text())
    geo = [r for r in scan['rows'] if mp.mpf(r['delta']) in [mp.mpf('0.0001') * 2**k for k in range(7)]]
    q = [mp.mpf(str(sp.N(sp.sympify(e).subs(cs, sp.Float(scan['c'], 70)), 65))) for e in ser['eta_h_over_delta8']]
    xs = [mp.mpf(r['delta']) for r in geo]
    ys = [(mp.mpf(r['eta_h']) / x**8 - mp.fsum(q[i] * x**i for i in range(8))) / x**8 for r, x in zip(geo, xs)]
    V2 = mp.matrix([[1, mp.log(x), x, x**2, x**3, x * mp.log(x)] for x in xs[:6]]); c2 = mp.lu_solve(V2, mp.matrix(ys[:6]))
    V1 = mp.matrix([[1, x, x**2, x**3, x**4, x**5] for x in xs[:6]]); c1 = mp.lu_solve(V1, mp.matrix(ys[:6]))
    Bpred = mp.mpf(kappa.numerator) / kappa.denominator / 2 * mp.mpf(3) / 136 * q[0]**2
    # residual of the pure-polynomial model on the 7th point (out-of-sample)
    pred7_poly = mp.fsum(c1[k] * xs[6]**k for k in range(6))
    pred7_log = c2[0] + c2[1] * mp.log(xs[6]) + c2[2] * xs[6] + c2[3] * xs[6]**2 + c2[4] * xs[6]**3 + c2[5] * xs[6] * mp.log(xs[6])
    etah[label] = dict(B_predicted=mp.nstr(Bpred, 10), B_fitted=mp.nstr(c2[1], 10), rel_diff=mp.nstr(abs(c2[1] - Bpred) / abs(Bpred), 4),
                       out_of_sample_err_log_model=mp.nstr(abs(pred7_log - ys[6]), 4), out_of_sample_err_poly_model=mp.nstr(abs(pred7_poly - ys[6]), 4),
                       scaled_remainders=[mp.nstr(y, 12) for y in ys])
out['eta_h_delta8_log_check'] = etah
out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(HERE / 'RESONANCE_PROBE.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({k: v for k, v in out.items() if k not in ('numerical_delta9_fits',)}, indent=1))
for l, v in fits.items():
    for n, w in v.items(): print(l, n, {k: w[k] for k in w if k != 'scaled_remainders'})
