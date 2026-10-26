#!/usr/bin/env python3
"""Small exact algebra controls; no actual fields or project imports."""
from fractions import Fraction as Q
from itertools import product
import json


def need(ok, label):
    if not ok:
        raise ValueError(label)


def eps(a,b,c):
    if len({a,b,c}) < 3:
        return 0
    return 1 if (a,b,c) in ((0,1,2),(1,2,0),(2,0,1)) else -1


def delta(a,b):
    return int(a == b)


# Product metric dz^2+dS4: explicit 5D Schouten reconstruction, all625 slots.
sign = (1,-1,1,1,1)
metric = lambda a,b: sign[a]*delta(a,b)
schouten = lambda a,b: Q(-1,2)*delta(a,0)*delta(b,0) + (Q(1,2)*metric(a,b) if a and b else 0)
for a,b,c,d in product(range(5),repeat=4):
    riemann = metric(a,c)*metric(b,d)-metric(a,d)*metric(b,c) if all(i > 0 for i in (a,b,c,d)) else 0
    wedge = metric(a,c)*schouten(b,d)-metric(a,d)*schouten(b,c)-metric(b,c)*schouten(a,d)+metric(b,d)*schouten(a,c)
    need(riemann == wedge, 'PRODUCT_WEYL_ZERO')

# Scalar curls of trace/Hessian reduce to kappa_bc v_d-kappa_bd v_c.
# The Hessian reduction uses the exact round-sphere commutator identity
# D_c Hess_bd Y-D_d Hess_bc Y=kappa_bc D_dY-kappa_bd D_cY.
for component in range(3):
    v = [int(i == component) for i in range(3)]
    curl = lambda b,c,d: delta(b,c)*v[d]-delta(b,d)*v[c]
    for a,b in product(range(3),repeat=2):
        symdual = sum(Q(1,4)*(eps(a,c,d)*curl(b,c,d)+eps(b,c,d)*curl(a,c,d)) for c,d in product(range(3),repeat=2))
        need(symdual == 0, 'SCALAR_SYMDUAL_ZERO')

# Check the sign/factor for every independent component of D_c T_ab,
# using C_b0cd=D_c T_bd-D_d T_bc and normalized symmetrization.
for i,j,k in ((i,j,k) for i in range(3) for j in range(i,3) for k in range(3)):
    deriv = lambda a,b,c: int(tuple(sorted((a,b))) == (i,j) and c == k)
    curvature = lambda b,c,d: deriv(b,d,c)-deriv(b,c,d)
    for a,b in product(range(3),repeat=2):
        H = sum(Q(1,4)*(eps(a,c,d)*curvature(b,c,d)+eps(b,c,d)*curvature(a,c,d)) for c,d in product(range(3),repeat=2))
        symcurl = sum(Q(1,2)*(eps(a,c,d)*deriv(b,d,c)+eps(b,c,d)*deriv(a,d,c)) for c,d in product(range(3),repeat=2))
        need(H == symcurl, 'TT_CODAZZI_SYMCURL_FACTOR_SIGN')

print(json.dumps({'status':'PASS_EXACT_SCREEN_ANGULAR_ALGEBRA_CONTROLS',
                  'checks':['625 product-Weyl components','scalar trace/Hessian reduced-curl symmetric dual','18 independent angular tensor derivative components'],
                  'limits':'Conditional tensor identities only; the sphere Hessian commutator is stated analytically. No registered basis curl matrix, observer anchoring, quadratic map, or actual curvature interval computed.'},sort_keys=True))
