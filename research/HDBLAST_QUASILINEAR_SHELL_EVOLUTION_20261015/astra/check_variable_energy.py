#!/usr/bin/env python3
"""Finite algebra checks only; not a proof of PDE existence or tame estimates."""
import argparse
import json
import sys
from pathlib import Path
import sympy as s

parser = argparse.ArgumentParser()
parser.add_argument('--mutation', choices=['none', 'reverse_flux', 'drop_shift', 'mix_DN'], default='none')
parser.add_argument('--output', type=Path)
args = parser.parse_args()
t, z, x = s.symbols('t z x', real=True)
u = s.Function('u')(t, z, x)
a, bz, bx, czz, czx, cxx = [s.Function(n)(t, z, x) for n in ('a', 'bz', 'bx', 'czz', 'czx', 'cxx')]
ut, uz, ux = [s.diff(u, q) for q in (t, z, x)]
p = a*s.diff(u,t,2)-2*bz*s.diff(u,t,z)-2*bx*s.diff(u,t,x)-czz*s.diff(u,z,2)-2*czx*s.diff(u,z,x)-cxx*s.diff(u,x,2)
e = (a*ut**2+czz*uz**2+2*czx*uz*ux+cxx*ux**2+u**2)/2
qz = bz*ut**2+(czz*uz+czx*ux)*ut
qx = bx*ut**2+(czx*uz+cxx*ux)*ut
bulk = (s.diff(a,t)/2-s.diff(bz,z)-s.diff(bx,x))*ut**2
bulk += (s.diff(czz,t)*uz**2+2*s.diff(czx,t)*uz*ux+s.diff(cxx,t)*ux**2)/2
bulk -= (s.diff(czz,z)+s.diff(czx,x))*uz*ut+(s.diff(czx,z)+s.diff(cxx,x))*ux*ut
bulk += u*ut
flux_sign = -1 if args.mutation == 'reverse_flux' else 1
if args.mutation == 'drop_shift':
    qz -= bz*ut**2
    qx -= bx*ut**2
local_residual = s.expand(s.diff(e,t)-p*ut-bulk-flux_sign*(s.diff(qz,z)+s.diff(qx,x)))
local_identity = local_residual == 0

# Two orthogonal components suffice to exhibit the 4+12 split cancellation.
dz, nt, theta, gamma, shift, cross, nx = s.symbols('dz nt theta gamma shift cross nx', real=True)
normal_derivative = s.Matrix([dz, 0])
time_derivative = s.Matrix([0, nt])
principal_fiber = s.Matrix([[gamma, theta], [theta, gamma]]) if args.mutation == 'mix_DN' else gamma*s.eye(2)
homogeneous_boundary_flux = -(normal_derivative.T*principal_fiber*time_derivative)[0]
mixed_flux_zero = s.expand(homogeneous_boundary_flux) == 0

# Nonzero normal shift and normal/angular cross terms are real obstructions
# to this elementary zero-flux argument, even with homogeneous D/N rows.
shift_obstruction = -shift*nt**2
cross_obstruction = -cross*nx*nt
matrix_obstruction = -(normal_derivative.T*s.Matrix([[gamma,theta],[theta,gamma]])*time_derivative)[0]
controls_detect_obstruction = all(s.expand(v) != 0 for v in (shift_obstruction, cross_obstruction, matrix_obstruction))

# Fundamental theorem of calculus jet bookkeeping: the j-th lower jet
# H^(m-1-j) differentiates to the (j+1)-st high jet with exactly that exponent.
jet_bookkeeping = all((m-(j+1)) == (m-1-j) for m in (12,14) for j in range(m))
# Every quadratic principal commutator total is <=15, hence one factor
# has <=7 derivatives; X_10 leaves >=3 spatial derivatives for L-infinity.
tame_count = all(min(q,m+1-q) <= 7 and 10-min(q,m+1-q) >= 3
                 for m in (12,14) for q in range(1,m))
checks = {'local_energy_identity': local_identity,
          'homogeneous_split_flux_zero': mixed_flux_zero,
          'excluded_terms_are_nonzero_controls': controls_detect_obstruction,
          'FTC_jet_index_identity': jet_bookkeeping,
          'finite_tame_derivative_count': tame_count,
          'variation_and_smoothing_orders_10_through_15': all(
              min(q,r+1-q) <= (r+1)//2 and
              ((r+1)//2+3)-min(q,r+1-q) >= 3 and
              ((r+1)//2+3) <= r-1
              for r in range(10,16) for q in range(1,r))}
result = {'scope': 'finite symbolic identities; no PDE existence certification',
          'mutation': args.mutation, 'checks': checks, 'passed': all(checks.values())}
if args.output:
    args.output.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
sys.exit(0 if result['passed'] else 1)
