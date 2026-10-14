"""Boost-profile-weighted l=0 integrand W(s) = s^{3/2} Omega^5 f0(s) for the
commuting state, reduced modulo total s-derivatives to a quadratic form."""
import json, pickle
import sympy as sp
from engine import O, F, register, bg_rules, dS
from sphere import sphere_average
d = pickle.load(open('../results/sources_cc.pkl', 'rb'))
f0 = sp.together(sphere_average(d['SH']))
s = sp.Symbol('s_')
G = register('G'); K = register('K')
rules, _ = bg_rules()
Wd = sp.expand(sp.together(s**sp.Rational(3, 2) * O[0]**5 * f0) )
Wd = sp.together(Wd)
print('f0 =', sp.factor(f0))
# split into GG, GK, KK parts
num, den = sp.fraction(sp.factor(Wd))
parts = {}
for lab, sub in [('GG', {k: 0 for k in K}), ('KK', {g: 0 for g in G})]:
    parts[lab] = sp.factor(Wd.subs(sub))
parts['GK'] = sp.factor(sp.expand(Wd - parts['GG'] - parts['KK']))
for k, v in parts.items():
    print(k, ':', v)
json.dump({'f0': str(sp.factor(f0)), **{k: str(v) for k, v in parts.items()}},
          open('../results/l0_weight.json', 'w'), indent=2)
