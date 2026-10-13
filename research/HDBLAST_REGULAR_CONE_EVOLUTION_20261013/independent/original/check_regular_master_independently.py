#!/usr/bin/env python3
"""Exact algebra for the independent cone-master derivation.

No asserts: checks remain effective under python -O. This checks algebra,
not the cited local hyperbolic or Cauchy-Kowalevskaya theorems.
"""
import argparse
import json
import sys
import sympy as s

parser = argparse.ArgumentParser()
parser.add_argument('--mutation', choices=['master-potential', 'scalar-source', 'null-jet'])
args = parser.parse_args()
checks = []

def verify(name, expression, expect_zero=True):
    value = s.factor(s.cancel(s.together(expression)))
    passed = (value == 0) if expect_zero else (value != 0)
    checks.append({'name': name, 'passed': bool(passed), 'residual': str(value)})

r = s.symbols('r', positive=True)
Om = s.Function('Omega')(r)
V = s.Function('V')(r)  # V(r) = Phi_s(r^2)
F = s.Function('F')(r)
B = s.Function('Box_gamma_F')(r)
rho = r*Om
p = 2*r*V/Om
P = r**4*F
Dy = lambda z: s.diff(z, r)/Om
alpha = 1/(rho**2*p**2)
q = 2/(3*rho**2)
w = 1/(rho**4*p**2)
master = -Dy(alpha*Dy(P)) + q*P - w*r**4*(B+4*F)
L = (s.diff(Om,r)/Om + 2*s.diff(V,r)/V)/(2*r)
potential_coefficient = 9 if args.mutation == 'master-potential' else 8
regular = s.diff(F,r,2)+4*s.diff(F,r)/r+B/r**2-2*r*L*s.diff(F,r) \
    -(potential_coefficient*L+s.Rational(8,3)*r**2*V**2)*F
verify('exact radial master to regular ambient operator', master+regular/(4*Om**2*V**2))
verify('wrong potential coefficient rejected', master+(regular-L*F)/(4*Om**2*V**2), False)
verify('wrong drift coefficient rejected', master+(regular+r*L*s.diff(F,r))/(4*Om**2*V**2), False)
verify('scalar reconstruction', -3*Dy(P)/(rho**2*p)+3*(4*F+r*s.diff(F,r))/(2*Om**2*V))
Psi = P/rho**2
verify('metric eta coefficient', 2*Psi*Om**2-2*r**2*F)
verify('metric X_A X_B coefficient', -6*Psi*Om**2/r**2+6*F)

# Independent ADM constraints and scalar residual, with lambda representing
# the commuting operator D=Box_gamma+4. Polynomial identity in lambda means
# no harmonic division, no projection away from exceptional subspaces.
y = s.symbols('y')
a, pp, H, U1 = [s.Function(n)(y) for n in ['a','p','H','U1']]
U2, lam, psi, chi = s.symbols('U2 lambda psi chi')
bg = {s.diff(a,y):-H-pp**2/3, s.diff(pp,y):U1-4*a*pp,
      s.diff(H,y):-2*a*H, s.diff(U1,y):U2*pp}
psip = -2*a*psi-pp*chi/3
chip = (U1/pp-4*a)*chi+(3*H*lam/pp-2*pp)*psi

def total_derivative(expr):
    return s.diff(expr,y).subs(bg, simultaneous=True) + s.diff(expr,psi)*psip + s.diff(expr,chi)*chip

chipp = total_derivative(chip)
N = -2*psi
d = psip-a*N
hamiltonian = 12*a*d+3*H*lam*psi-pp*(chip-pp*N)+U1*chi
verify('unreduced Hamiltonian constraint', hamiltonian)
verify('unreduced momentum coefficient', -3*d-pp*chi)
verify('unreduced trace-free lapse coefficient', N+2*psi)
scalar_source = 3 if args.mutation == 'scalar-source' else 4
scalar = chipp+4*a*chip+H*(lam-4)*chi-U2*chi+6*pp*psip+scalar_source*U1*psi
verify('full linearized scalar residual', scalar)
verify('wrong scalar source rejected', scalar-U1*psi, False)

# At q=(T,R)=(t0,t0), compact cutoff Cauchy data are spatially constant
# on an annulus: F1=a0, F1_T=0 and F2=0, F2_T=b0. The master supplies TT.
t0, a0, b0, L0, C, ellam, f0 = s.symbols('t0 a0 b0 L0 C lambda_l f0', nonzero=True)

def cone_scalar_jet(value, ft, fr, ftt, ftr, frr):
    field = -C*(4*value+t0*(ft+fr))
    euler_derivative = 4 if args.mutation == 'null-jet' else 5
    deriv = -C*s.Rational(1,2)*(euler_derivative*(ft+fr)+t0*(ftt+2*ftr+frr))
    return field, deriv
x1,x1v = cone_scalar_jet(a0,0,0,(-ellam/t0**2-8*L0)*a0,0,0)
x2,x2v = cone_scalar_jet(0,b0,0,-2*L0*t0*b0,0,0)
wronskian = s.expand(x1*x2v-x2*x1v)
verify('actual Cauchy-jet null Wronskian', wronskian-(10+ellam/2)*C**2*a0*b0)
verify('nonzero Wronskian', wronskian, False)
verify('background cancellation in Wronskian', s.diff(wronskian,L0))
verify('l=2 null Wronskian', wronskian.subs(ellam,8)-14*C**2*a0*b0)
z1 = 4*f0*(2*t0*x1v-x1)
z2 = 4*f0*(2*t0*x2v-x2)
verify('gauge invariant relational scalar first jet', z1-4*f0*C*a0*(ellam+4+8*t0**2*L0))
verify('gauge invariant relational scalar second jet', z2-4*f0*C*t0*b0*(-4+2*t0**2*L0))
verify('nonzero relational scalar first jet', z1, False)
verify('nonzero relational scalar second jet', z2, False)

# Independent conjugation check: A_s/A=-L and F=A^(-1/2) psi.
# In five dimensions Box_eta k(s)=4s k_ss+10 k_s.
scoord, LL, LLs, ff = s.symbols('s L L_s f')
k_s_over_k = LL/2
k_ss_over_k = LLs/2+LL**2/4
conjugated_potential = 4*scoord*k_ss_over_k+10*k_s_over_k \
    -4*LL*scoord*k_s_over_k-8*LL-s.Rational(8,3)*scoord*ff**2
Vcan = 3*LL-2*scoord*LLs+scoord*LL**2+s.Rational(8,3)*scoord*ff**2
verify('canonical scalar conjugation potential', conjugated_potential+Vcan)
verify('wrong canonical scalar sign rejected', conjugated_potential-Vcan, False)


# Exact density scaling from the inherited radial current, with r^2=uv.
u,v = s.symbols('u v', positive=True)
om = s.symbols('Omega', positive=True)
rr = s.sqrt(u*v)
cosht = (u+v)/(2*rr)
# rho^4 cosh^3(t) * (v/rho) * dt/dv, where dt/dv=1/v at fixed r.
scale = (rr*om)**4*cosht**3*(v/(rr*om))/v
verify('null scalar current density coefficient', s.limit(scale,u,0).subs(om,1)-v**3/8)

result = {'schema_version':1, 'sympy':s.__version__, 'optimized':not __debug__,
          'mutation':args.mutation, 'passed':all(c['passed'] for c in checks),
          'check_count':len(checks), 'checks':checks}
print(json.dumps(result, indent=2))
sys.exit(0 if result['passed'] else 1)
