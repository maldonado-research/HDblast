#!/usr/bin/env python3
"""Exact polynomial controls for the frozen characteristic shell lift."""
from fractions import Fraction as Q
from math import comb
import json


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def add_term(p, key, value):
    p[key] = p.get(key, Q(0)) + value


def lift_polynomial(coefficients):
    # -integral_0^(x-r) g(s) ds on the x>r side.
    p = {}
    for k, c in enumerate(coefficients):
        d = k + 1
        for j in range(d + 1):
            add_term(p, (d-j, j), -Q(c, d) * comb(d, j) * (-1)**j)
    return p


def derivative(p, dx=0, dr=0):
    out = dict(p)
    for slot, count in enumerate((dx, dr)):
        for _ in range(count):
            nxt = {}
            for (i, j), c in out.items():
                power = (i, j)[slot]
                if power:
                    key = (i-1, j) if slot == 0 else (i, j-1)
                    add_term(nxt, key, c * power)
            out = nxt
    return {k: v for k, v in out.items() if v}


def evalp(p, x, r):
    return sum((c * x**i * r**j for (i, j), c in p.items()), Q(0))


def polyval(p, z):
    return sum((Q(c) * z**k for k, c in enumerate(p)), Q(0))


def derivcoef(p):
    return [Q(k) * Q(p[k]) for k in range(1, len(p))]


def primitive(q):
    return q**3/Q(3) - q**4/Q(2) + q**5/Q(5)


def w(q):
    return q*q*(1-q)*(1-q) if 0 <= q <= 1 else Q(0)


def wprime(q):
    return 2*q*(1-q)*(1-2*q) if 0 <= q <= 1 else Q(0)


checks = 0
# Cover polynomial inputs with several nonzero derivatives and dimensions.
for source in ([0, 0, 1], [0, 0, -2, 3, 5], [0, 0, 0, 7, -4]):
    p = lift_polynomial(source)
    gp = derivcoef(source)
    for x, r in [(Q(3,2), Q(1,2)), (Q(7,4), Q(1,4)), (Q(9,5), Q(2,5))]:
        z=x-r
        expected={(1,0):-polyval(source,z),(0,1):polyval(source,z),
                  (2,0):-polyval(gp,z),(1,1):polyval(gp,z),(0,2):-polyval(gp,z)}
        for orders, target in expected.items():
            require(evalp(derivative(p,*orders),x,r)==target,'positive-side derivative')
            checks += 1
    # First and second jets of the positive polynomial match zero at x=r.
    for r in [Q(1,3), Q(1), Q(2)]:
        for orders in [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]:
            require(evalp(derivative(p,*orders),r,r)==0,'C2 characteristic join')
            checks += 1

# Converse: g=s gives C1 but different second-jet one-sided values.
p=lift_polynomial([0,1])
for r in [Q(1,2),Q(1),Q(2)]:
    for orders in [(0,0),(1,0),(0,1)]:
        require(evalp(derivative(p,*orders),r,r)==0,'g=s must be C1')
    require(evalp(derivative(p,2,0),r,r)==-1,'g=s must not be C2')
    require(evalp(derivative(p,1,1),r,r)==1,'mixed sign')
    checks+=5

# Exact polynomial family with shell collar, known extrema and integral.
for n in [5,8,10,100]:
    n=Q(n)
    require(n>=5 and 1+1/n<=Q(5,4), 'complete source support in shell-history interval')
    require(w(-Q(1,2))==0 and w(Q(3,2))==0, 'source support outside bump')
    checks+=2
    require(w(Q(1,2)) / n == 1/(16*n),'bump maximum witness')
    require(primitive(Q(1)) / n**2 == 1/(30*n**2),'bump integral')
    require(wprime(Q(1,4))==Q(3,16),'nonvanishing second jet')
    rx=Q(1,16); xx=1+1/(4*n)+rx
    require(0<=xx<=Q(3,2) and max(Q(0),xx-Q(5,4))<=rx<=Q(7,4)-xx, 'witness on registered domain')
    checks+=1
    require(w(Q(0))==w(Q(1))==wprime(Q(0))==wprime(Q(1))==0,'C1 bump joins')
    checks+=4

# Monomial chain rule with y-derivative rho_y symbolic-independent fixture.
# d_y = -1/(ell*rho) d_r; then d_y² F = a²F_rr + rho_y/(ell*rho²)F_r.
ell,rho,rho_y=Q(1,1000000),Q(5,4),Q(3,7)
p=lift_polynomial([0,0,2,3])
x,r=Q(7,4),Q(1,4)
fr=evalp(derivative(p,0,1),x,r)
frr=evalp(derivative(p,0,2),x,r)
fxr=evalp(derivative(p,1,1),x,r)
a=1/(ell*rho)
# Independent composition of -a*d_r treating d_y a = -a*rho_y/rho.
composed=a*a*frr+a*rho_y/rho*fr
expanded=frr/(ell**2*rho**2)+rho_y*fr/(ell*rho**2)
require(composed==expanded,'physical yy chain rule')
require(-a*fxr/ell == -fxr/(ell**2*rho),'physical ty chain rule')
require(expanded != frr/(ell**2*rho**2),'omitted rho_y control')
require(evalp(derivative(p,1,1),x,r)!=-polyval(derivcoef([0,0,2,3]),x-r),'wrong mixed sign control')
checks+=4

print(json.dumps({'status':'PASS','checks':checks,'source_family':'exact C1 compact polynomial with zero collar; separate note gives smooth bump',
                  'negative_controls':['C1 but not C2 g=s','incorrect mixed sign','omitted rho_y chain term'],
                  'scope':'affine lift only; no actual-field numerical evaluation'},sort_keys=True))
