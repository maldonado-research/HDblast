"""Independent numerical test (no original-engine code) of the l = 0 sector:
 (G) radial gauge data are annihilated by DH (diffeomorphism invariance);
 (W) int s^{3/2} Om^5 DH[alpha delta + a x x^T, Psi] ds
       == -(3/2) int Om sqrt(s) (3 alpha + s a) ds  for compactly supported data;
 (P) Psi-only data give a vanishing weighted integral (scalar cannot carry m).
Uses indep_check's numerical evaluator on the exact test background."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp, mpmath as mp
import indep_check as ic
S = ic.S
s1, s2 = mp.mpf('0.5'), mp.mpf('1.5')
def bumped(poly):
    J = ic.jets_fn(poly)
    def f(sv):
        if not (s1 < sv < s2):
            return [mp.mpf(0)] * 3
        b = mp.exp(-1 / ((sv - s1) * (s2 - sv)))
        # derivatives of bump*poly up to order 2 via mp.diff on the product
        g = lambda t: mp.exp(-1 / ((t - s1) * (s2 - t))) * J(t)[0]
        return [g(sv), mp.diff(g, sv), mp.diff(g, sv, 2)]
    return f
alphaF = bumped(1 + S); aF = bumped(2 - 3 * S**2); PsiF = bumped(S - sp.Rational(1, 3)); zetaF = bumped(1 + 2 * S)
def make(alpha=None, a=None, Psi=None):
    def field(x):
        x = [mp.mpf(c) for c in x]; sv = sum(c * c for c in x); o = ic.OmJ(sv)
        al = alpha(sv)[0] if alpha else 0; aa = a(sv)[0] if a else 0
        ps = Psi(sv)[0] if Psi else 0
        h = [[al * (1 if i == j else 0) + aa * x[i] * x[j] for j in range(4)] for i in range(4)]
        g_c = [[[o[0]**2 * (1 if i == j else 0), h[i][j], 0] for j in range(4)] for i in range(4)]
        return g_c, [[[0, 0, 0]] * 4] * 4, [0, ps, 0], [0, 0, 0]
    return field
def DH_at(field, sv):
    r = mp.sqrt(sv)
    xs = [r * mp.mpf('0.5'), r * mp.mpf('-0.5'), r * mp.mpf('0.5'), r * mp.mpf('0.5')]
    D = ic.numeric_D(field, xs, {'g': [0, 1], 'pi': []})
    H, _ = ic.constraints_at(D, sv)
    return H[1]
out = {}
# (G) gauge
def ag(sv):
    z = zetaF(sv); o = ic.OmJ(sv)
    return [2 * o[0]**2 * z[0] + 4 * sv * o[0] * o[1] * z[0]]
def aag(sv):
    z = zetaF(sv); o = ic.OmJ(sv)
    return [4 * o[0]**2 * z[1]]
def psig(sv):
    z = zetaF(sv); f = ic.PhpJ(sv)
    return [2 * sv * f[0] * z[0]]
gres = []
for sv in ('0.7', '1.0', '1.3'):
    sv = mp.mpf(sv)
    val = DH_at(make(ag, aag, psig), sv)
    gres.append(mp.nstr(abs(val) / abs(ag(sv)[0]), 5))
out['gauge_residual_rel'] = gres
print('gauge', gres, flush=True)
mp.mp.dps = 25
Om = lambda sv: ic.OmJ(sv)[0]
lhs = mp.quad(lambda sv: sv**mp.mpf(1.5) * Om(sv)**5 * DH_at(make(alphaF, aF, PsiF), sv), [s1, 1, s2], maxdegree=5)
rhs = mp.quad(lambda sv: -mp.mpf(3) / 2 * Om(sv) * mp.sqrt(sv) * (3 * alphaF(sv)[0] + sv * aF(sv)[0]), [s1, 1, s2], maxdegree=7)
out['weighted_identity_lhs'] = mp.nstr(lhs, 15)
out['weighted_identity_rhs'] = mp.nstr(rhs, 15)
out['weighted_identity_rel_diff'] = mp.nstr(abs(lhs - rhs) / abs(rhs), 5)
print(out, flush=True)
pint = mp.quad(lambda sv: sv**mp.mpf(1.5) * Om(sv)**5 * DH_at(make(Psi=PsiF), sv), [s1, 1, s2], maxdegree=5)
pscale = mp.quad(lambda sv: abs(sv**mp.mpf(1.5) * Om(sv)**5 * DH_at(make(Psi=PsiF), sv)), [s1, 1, s2], maxdegree=4)
out['psi_only_weighted_integral_rel'] = mp.nstr(abs(pint) / pscale, 5)
print(json.dumps(out, indent=2))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'indep_l0.json'), 'w'), indent=2)
