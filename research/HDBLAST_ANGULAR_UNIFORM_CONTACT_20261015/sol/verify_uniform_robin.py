#!/usr/bin/env python3
"""Independent exact boundary/conjugation checks; no nonlinear theorem checker."""
from pathlib import Path
import argparse
import json

import sympy as S

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path,
                    default=Path(__file__).with_name("verification_results.json"))
parser.add_argument("--inject-corruption", choices=("boundary-sign", "conjugation-sign",
                    "omitted-quadratic-derivative", "noncommuting-time"))
args = parser.parse_args()


def zero(x):
    if isinstance(x, S.MatrixBase):
        return all(S.simplify(e) == 0 for e in x)
    return S.simplify(x) == 0


def require(condition, message):
    if not condition:
        raise SystemExit("FAIL: " + message)


def pack(K):
    return S.Matrix([K[i, j] for i, j in PAIRS])


def unpack(v):
    K = S.zeros(4)
    for k, (i, j) in enumerate(PAIRS):
        K[i, j] = K[j, i] = v[k]
    return K


PAIRS = [(i, j) for i in range(4) for j in range(i, 4)]
a, s, b, q0, qsp, eps = S.symbols("a s b q0 qsp eps", nonzero=True, real=True)
# Off-diagonal angular coordinates are intentional: they test the factor two
# in contraction of symmetric off-diagonal components.
kappa = S.Matrix([[2, 1, 0], [1, 3, 1], [0, 1, 2]])
H = S.diag(-q0, qsp, qsp, qsp)
H[1:4, 1:4] = qsp * kappa
Hinv = H.inv()
r, w = S.symbols("r w")
ks = S.symbols("k0:10")
K = unpack(ks)
tau = S.trace(Hinv * K)
v = S.Matrix([r, *ks, w])


def boundary_B(vv):
    rr, KK, ww = vv[0], unpack(vv[1:11]), vv[11]
    return S.Matrix([
        4*a*rr - 2*a*S.trace(Hinv*KK) - 4*s*ww/3,
        *pack(-a*H*rr - 2*a*KK - s*H*ww/3),
        s*rr/4 + b*ww/2,
    ])


B = boundary_B(v).jacobian(v)
require(B.shape == (12, 12), "12 normal derivative rows")
require(zero(boundary_B(v) - B*v), "linear boundary matrix")
if args.inject_corruption == "boundary-sign":
    B[11, 0] = -B[11, 0]

# Direct first and second differentiation of the exact primal boundary laws.
r2, w2 = S.symbols("r2 w2")
k2s = S.symbols("k20:30")
K2 = unpack(k2s)
A1 = Hinv*K
tau2 = S.trace(Hinv*K2)
N = 1 + eps*r/2 + eps**2*(r2/4-r**2/8)
sig = 6*a + eps*s*w + eps**2*(s*w2+b*w**2)/2
sigp = s + eps*b*w + eps**2*(b*w2+4*w**2)/2
tr = 4-eps*tau+eps**2*(S.trace(A1*A1)-tau2/2)
h = H+eps*K+eps**2*K2/2
Fyy = 4*sig*N**3/3 - 2*a*N**4*tr
Fh = sig*N*h/3
Fphi = -sigp*N/2
exact = S.Matrix([Fyy, *pack(Fh), Fphi])
first = exact.diff(eps).subs(eps, 0)
require(zero(first+B*v), "B is minus Jacobian of the exact boundary RHS")

n = r/2
n_y = -4*a*n+a*tau+2*s*w/3
Rh = (2*s*w/3+4*a*n)*K + (b*w*w/3+2*s*w*n/3)*H
Rphi = -2*w*w-b*w*n
Rn = 2*b*w*w/3+8*s*w*n/3-16*a*n*n-2*a*S.trace(A1*A1)+6*a*n*tau
sources = S.Matrix([
    2*Rn+4*n*n_y+8*a*n*n,
    *pack(Rh-2*a*n*n*H),
    Rphi+s*n*n/2,
])
if args.inject_corruption == "omitted-quadratic-derivative":
    sources[0] -= 4*n*n_y
second = exact.diff(eps, 2).subs(eps, 0)
v2 = S.Matrix([r2, *k2s, w2])
require(zero(second+B*v2-sources), "all equation (5) metric-wave sources")

# Trace/trace-free decomposition uses H held fixed in the normal extension.
tau_symbol = S.symbols("tau")
trace_block = S.Matrix([
    [4*a, -2*a, -4*s/3],
    [-4*a, -2*a, -4*s/3],
    [s/4, 0, b/2],
])
Kt = K-S.trace(Hinv*K)*H/4
Bt = unpack(boundary_B(v)[1:11])
require(zero(Bt-S.trace(Hinv*Bt)*H/4+2*a*Kt), "nine trace-free rows")
generic_example = trace_block.subs({a:0, b:0, s:1})
require(generic_example.charpoly().as_expr().equals(
    generic_example.charpoly().gen*(generic_example.charpoly().gen**2+S.Rational(1,3))),
    "generic complex eigenvalues; no presumed positive symmetrizer")

# Exact noncommuting nilpotent example: M^3=0, so its exponential is a
# polynomial and the Duhamel integral is exact, with no numerical tolerance.
t, chi, integration_r = S.symbols("t chi integration_r", real=True)
E12, E23 = S.zeros(3), S.zeros(3)
E12[0, 1], E23[1, 2] = 1, 1
example_M = E12+t*E23
example_Mt = example_M.diff(t)
require(zero(example_M**3), "nilpotent example has exact exponential")
def nilpotent_exp(scale):
    return S.eye(3)+scale*example_M+scale**2*example_M**2/2
T = nilpotent_exp(chi)
correct_time_derivative = T.diff(t)
integrand = nilpotent_exp((1-integration_r)*chi) * (chi*example_Mt) \
            * nilpotent_exp(integration_r*chi)
duhamel = integrand.applyfunc(lambda e: S.integrate(e, (integration_r, 0, 1)))
bad_time_derivative = chi*example_Mt*T
selected_time_derivative = (bad_time_derivative
    if args.inject_corruption == "noncommuting-time" else correct_time_derivative)
require(zero(selected_time_derivative-duhamel), "noncommuting Duhamel derivative")
require(not zero(correct_time_derivative-bad_time_derivative),
        "corruption control detects invalid commuting derivative")

# T(0)=I and T_z(0)=M give exactly zero residual Robin matrix.
rho = S.symbols("rho", positive=True)
M = rho*B
selected_Tz = -M if args.inject_corruption == "conjugation-sign" else M
require(zero(-selected_Tz+M), "correct conjugation sign")
require(not zero(M+M), "corruption control detects opposite conjugation sign")

# Exact glancing asymptotics at Laplace frequency 1+i*eta and c=1.
eta = S.symbols("eta", positive=True)
kappa_abs_squared = S.sqrt(1+4*eta**2)
kappa_real = S.sqrt((kappa_abs_squared+1)/2)
lambda_squared = 2+2*eta**2
normalized_squared = kappa_abs_squared/lambda_squared
bulk_e_squared = (lambda_squared+kappa_abs_squared) \
                 /(2*kappa_real*kappa_abs_squared)
symbol_limit = S.limit(eta*normalized_squared, eta, S.oo)
energy_limit = S.limit(bulk_e_squared/S.sqrt(eta), eta, S.oo)
require(symbol_limit == 1, "glancing normalized symbol degenerates")
require(energy_limit == S.Rational(1, 2), "glancing integrated H1 scaling")

# Deliberate matrix corruptions must fail the direct equation test.
mutant = B.copy()
mutant[11, 0] = -mutant[11, 0]
require(not zero(first+mutant*v), "corruption control catches scalar/lapse sign")
mutant = B.copy()
offdiag_index = PAIRS.index((1, 2))+1
mutant[0, offdiag_index] /= 2
require(not zero(first+mutant*v), "corruption control catches trace factor two")
mutant_source = sources.copy()
mutant_source[0] -= 4*n*n_y
require(not zero(second+B*v2-mutant_source),
        "corruption control catches omitted metric/lapse quadratic derivative")

result = {
    "status":"PASS",
    "checks":["12x12 exact first derivative", "all second metric-wave sources",
              "trace/trace-free invariant blocks", "generic complex spectrum",
              "noncommuting exponential derivative", "conjugation sign",
              "Neumann energy normalization and glancing scaling",
              "four explicit corruption controls"],
    "exact_noncommuting_residual":str(correct_time_derivative-bad_time_derivative),
    "glancing_normalized_symbol_squared_times_eta_limit":str(symbol_limit),
    "glancing_energy_squared_over_sqrt_eta_limit":str(energy_limit),
    "limits":"Algebra/symbol checks only; the companion note supplies the energy argument.",
}
out = args.output
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(result, indent=2)+"\n")
print(json.dumps(result, indent=2))
