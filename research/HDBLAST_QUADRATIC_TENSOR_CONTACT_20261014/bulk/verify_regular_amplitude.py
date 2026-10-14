#!/usr/bin/env python3
"""Check the regular l=2 amplitude substitution in the derived tensor equation."""
import argparse
import json
from pathlib import Path
import sympy as s

ap=argparse.ArgumentParser()
ap.add_argument('--output',type=Path)
args=ap.parse_args()
t,r,z,A=s.symbols('t r z A',real=True)
J=s.Function('J')(t,z)
U=s.Function('U')(t,z)
H=r**2*J.subs(z,r**2)
F=r**2*U.subs(z,r**2)
actual=(-s.diff(H,t,2)+s.diff(H,r,2)+3*s.diff(H,r)/r
        +6*A*(t*s.diff(H,t)+r*s.diff(H,r))-8*H/r**2)/r**2
expected=(-s.diff(J,t,2)+4*z*s.diff(J,z,2)+16*s.diff(J,z)
          +6*A*(t*s.diff(J,t)+2*z*s.diff(J,z)+2*J)).subs(z,r**2)
scalar_actual=(t*s.diff(F,t)+r*s.diff(F,r)+4*F)/r**2
scalar_expected=(t*s.diff(U,t)+2*z*s.diff(U,z)+6*U).subs(z,r**2)
checks={}
for name,value in [('regular_tensor_operator',actual-expected),
                   ('regular_scalar_reconstruction',scalar_actual-scalar_expected)]:
    value=s.simplify(value)
    if value != 0:
        raise RuntimeError(f'{name}: {value}')
    checks[name]='passed'
out={'status':'passed','checks':checks,'scope':'Differential identities only; no origin or boundary existence theorem.'}
if args.output:
    args.output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
