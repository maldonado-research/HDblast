#!/usr/bin/env python3
"""Check the chosen unit-observer current by differentiating its components."""
import argparse
import json
from pathlib import Path
import sympy as s


def main():
    AT,AZ,BT,BZ,t,z,tt,tz,zz,U,Up,e2B=s.symbols('AT AZ BT BZ phiT phiZ phiTT phiTZ phiZZ U Uphi exp2B')
    core=(t*t+z*z)/2+e2B*U
    # Divide by the common factor exp(3(A-A0)-B)/kappa5^2.
    Et=(3*AT-BT)*core+t*tt+z*tz+e2B*(2*BT*U+Up*t)
    Fz=-(3*AZ-BZ)*t*z-tz*z-t*zz
    kg=zz-3*AT*t+3*AZ*z-e2B*Up
    raw=(Et+Fz).subs(tt,kg)
    G=(3*AT-BT)*z*z/2-(3*AT+BT)*t*t/2+(3*AT+BT)*e2B*U+BZ*t*z
    if s.simplify(raw-G)!=0:raise RuntimeError('unit-observer current identity failed')
    wrongs={'omit_metric_work':raw,
            'omit_BZ_term':raw-(G-BZ*t*z),
            'reverse_wall_flux_sign':(Et-Fz).subs(tt,kg)-G,
            'omit_bulk_potential_terms':raw-(G-(3*AT+BT)*e2B*U)}
    for name,value in wrongs.items():
        if s.simplify(value)==0:raise RuntimeError('nondiscriminating wrong formula: '+name)
    report=dict(scope='Finite algebra for the stated scalar stress current only.',identity_zero=True,
                incorrect_formula_residuals={k:str(s.factor(v)) for k,v in wrongs.items()})
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')


if __name__=='__main__':main()
