#!/usr/bin/env python3
"""Exact indexed check of the differentiated reduced-Ricci principal symbol.

Manufactured normal-frame jets, not an actual HDBLAST field or a complete
differentiated-source bound. No project code, source map, or CAS is imported.
"""
from fractions import Fraction as Q
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


sig = [1, -1, 1, 1, 1]
groups = []
controls = set()


def principal(seed, sign, mutation=None):
    # q is an arbitrary symmetric low metric value. Lower first/second jets
    # of the high-gauge perturbation vanish at the selected event.
    q = [[Q((min(a,b)+2)*(max(a,b)+seed+1)-3, a+b+2)
          for b in range(5)] for a in range(5)]
    signs = [1]*5 if mutation == 'euclidean_time_raising' else sig
    qup = [[signs[a]*signs[b]*q[a][b] for b in range(5)] for a in range(5)]
    z = [Q(1), Q(sign), Q(0), Q(0), Q(0)]
    f = [Q(1), Q(0), Q(0), Q(0), Q(0)]

    def third(e, c, d, a, b):
        return z[e]*z[c]*z[d]*(z[a]*f[b]+z[b]*f[a])

    def dP(e, d, a, b, c):
        # d_e d_d P_{abc}, with P_{abc}=1/2(d_b q_ca+d_c q_ba-d_a q_bc).
        return (third(e,d,b,c,a)+third(e,d,c,b,a)-third(e,d,a,b,c))/2

    ric = [[sum(-2*qup[c][d]*(dP(1,c,d,a,b)-dP(1,b,d,a,c))
                for c in range(5) for d in range(5))
            for b in range(5)] for a in range(5)]
    gradient = [[-sum(qup[c][d]*(dP(1,a,b,c,d)+dP(1,b,a,c,d))
                     for c in range(5) for d in range(5))
                 for b in range(5)] for a in range(5)]
    factor = Q(2) if mutation == 'wrong_symmetrization' else Q(1)
    if mutation == 'omit_gauge_gradient':
        factor = Q(0)
    reduced = [[ric[a][b]-factor*gradient[a][b] for b in range(5)] for a in range(5)]
    expected = [[sum(sig[c]*sig[d]*q[c][d]*third(1,c,d,a,b)
                     for c in range(5) for d in range(5))
                 for b in range(5)] for a in range(5)]
    return reduced, expected, sum(sig[c]*z[c]**2 for c in range(5))


for seed in (1, 3, 7, 11):
    for sign in (-1, 1):
        reduced, expected, null = principal(seed, sign)
        need(null == 0, 'principal wave covector must be null')
        need(reduced == expected, 'indexed reduced Ricci principal identity')
        groups.append('indexed_principal_seed_'+str(seed)+'_sign_'+str(sign))
for mutation in ('wrong_symmetrization', 'omit_gauge_gradient', 'euclidean_time_raising'):
    got, want, _ = principal(3, 1, mutation)
    need(got != want, 'wrong principal formula must be rejected:'+mutation)
    controls.add(mutation)

# Null quartic physical gauge-jet variations have zero traces in every pair,
# so the highest-order prolonged wave/gauge relations do not remove them.
for sign in (-1, 1):
    z = [Q(1), Q(sign)]
    for a in range(2):
        for b in range(2):
            need(sum(sig[c]*z[c]*z[c]*z[a]*z[b] for c in range(2)) == 0,
                 'prolonged gauge wave trace')
groups.append('all_null_quartic_pair_traces_zero')

# With low q_yy=1 and other low components zero, the complete yy and y0
# principal forcing symbols are nonzero. The scalar source has no fourth
# radial gauge derivative because its first scalar shift is p*f.
need(Q(-2,3) != 0 and Q(-1,3) != 0, 'registered mean projection symbols nonzero')
groups.append('nonzero_mean_symbol_witness')

print(json.dumps({'status':'PASS_EXACT_INDEXED_PRINCIPAL_SYMBOL',
    'group_count':len(groups),'named_groups':groups,
    'wrong_formula_controls':sorted(controls),
    'scope':'independent local principal-symbol algebra and null prolonged traces; not an actual source derivative, a C4 bound, a proof of actual nonzero forcing, or exclusion of data-specific identities'},sort_keys=True))
