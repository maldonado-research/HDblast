#!/usr/bin/env python3
"""Independent finite algebra checks of the displayed action/junction budget.
These identities do not establish a numerical solution or a global reservoir law.
"""
import argparse
import json
from pathlib import Path
import sympy as s


def main():
    positive=[];negative=[]
    def identity(name,expr):
        value=s.factor(expr)
        if value!=0:raise RuntimeError(name+': '+str(value))
        positive.append(name)
    def control(name,expr):
        value=s.factor(expr)
        if value==0:raise RuntimeError('nondiscriminating control: '+name)
        negative.append(dict(name=name,residual=str(value)))
    k2,sigma,sp,rho,p,J,v,H,rd=s.symbols('kappa2 sigma sigma_phi rho p J v H rho_dot', nonzero=True)
    ks=(sigma+k2*rho)/6;k0=(sigma-k2*(2*rho+3*p))/6;w=-(sp+k2*J)/2
    identity('metric_time_row',-3*ks+k2*(sigma/k2+rho)/2)
    identity('metric_spatial_row',-k0-2*ks-k2*(-sigma/k2+p)/2)
    identity('two_copy_scalar_boundary_variation',-2*w/k2-sp/k2-J)
    ward=rd+3*H*(rho+p)-J*v
    total=sp*v/k2+rd+3*H*(rho+p)
    identity('total_shell_flux_equals_matter_Ward',total+2*w*v/k2-ward)
    ksd=(sp*v+k2*rd)/6
    identity('Codazzi_equals_Ward',ksd+H*(ks-k0)+w*v/3-k2*ward/6)
    momentum=-3*ksd-3*H*(ks-k0)-w*v
    identity('normal_momentum_constraint_equals_Ward',momentum+k2*ward/2)
    control('omit_second_interior',total+w*v/k2-ward)
    control('reverse_scalar_junction_only',total-2*w*v/k2-ward)
    control('omit_tension_from_total_budget',rd+3*H*(rho+p)+2*w*v/k2-ward)
    # One particle mode: p redshifts, its mass changes, and decay is time dilated.
    n,q,m,g,phi_rate,Gamma,number_in=s.symbols('n q m g phi_rate Gamma number_in',nonzero=True)
    omega=s.sqrt(q*q+m*m)
    energy=n*omega;pressure=n*q*q/(3*omega);source=n*g/omega
    nd=-3*H*n-Gamma*m*n/omega+number_in
    # g = (1/2) d(m^2)/dphi, so m_dot=g phi_rate/m.
    ed=s.diff(energy,n)*nd+s.diff(energy,q)*(-H*q)+s.diff(energy,m)*(g*phi_rate/m)
    decay=Gamma*m*n;production=number_in*omega
    identity('particle_kinetic_energy_ledger',ed+3*H*(energy+pressure)-source*phi_rate+decay-production)
    control('production_without_paid_work',ed+3*H*(energy+pressure)-source*phi_rate+decay)
    control('massive_gas_as_radiation_pressure',ed+4*H*energy-source*phi_rate+decay-production)
    R=s.symbols('R');Rdot=-4*H*R+decay
    identity('decay_internal_transfer_cancels',ed+Rdot+3*H*(energy+pressure)+4*H*R-source*phi_rate-production)
    # Coordinate scalar energy current for xi=partial_T. Keep metric work.
    AT,AZ,BT,phT,phZ,phTT,phZZ,phTZ,U,Up,e2B=s.symbols('AT AZ BT phT phZ phTT phZZ phTZ U Up exp2B')
    # Common factor exp(3A)/kappa2 is omitted.
    E=phT**2/2+phZ**2/2+e2B*U
    dE=3*AT*E+phT*phTT+phZ*phTZ+e2B*(2*BT*U+Up*phT)
    dF=-3*AZ*phT*phZ-phTZ*phZ-phT*phZZ
    KG=-3*AT*phT+3*AZ*phZ-e2B*Up
    work=-s.Rational(3,2)*AT*phT**2+s.Rational(3,2)*AT*phZ**2+(3*AT+2*BT)*e2B*U
    identity('bulk_coordinate_energy_has_metric_work',(dE+dF).subs(phTT,phZZ+KG)-work)
    control('drop_bulk_metric_work',(dE+dF).subs(phTT,phZZ+KG))
    report=dict(scope='Finite symbolic identities only; no independent bulk flux measurement.',
                identities=positive,negative_controls=negative)
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(len(positive),'identities;',len(negative),'nonzero wrong-formula residuals')


if __name__=='__main__':
    main()
