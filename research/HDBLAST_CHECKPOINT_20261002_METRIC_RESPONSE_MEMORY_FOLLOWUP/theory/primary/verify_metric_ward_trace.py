#!/usr/bin/env python3
"""Exact closed metric Ward and independently differentiated finite anomaly."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import traceback
import sympy as S
import metric_contact_algebra as M


def audit():
    data=M.inventory(); closed=M.closed_quantities(data)
    L,V,H=M.L,M.V,M.H
    q=S.symbols('q0:5',real=True)
    pi2=2*S.pi**2
    def D(expr):return S.expand(L**2*S.diff(expr,L)-L*V*(1-V**2)*S.diff(expr,V)+sum(H[i+1]*S.diff(expr,H[i]) for i in range(6))+sum(q[i+1]*S.diff(expr,q[i]) for i in range(4)))
    checks=[]
    def zero(name,expr):
        residual=S.cancel(expr)
        if residual != 0:raise RuntimeError(name+': '+str(S.factor(residual)))
        checks.append(name)
    Q0=closed['physical_baseline_Q_times_2pi2']/pi2
    R0=closed['scaled_baseline_rho_times_2pi2']/pi2
    P0=closed['scaled_baseline_p_times_2pi2']/pi2
    g=data['g'];g1=data['g_prime'];g2=data['g_second']
    J5=V**3/6;J7=(V**3/3-V**5/5)/4;J9=(V**3/3-2*V**5/5+V**7/7)/8
    Rm=(3*L**2*q[0]-L*q[1])/2+L**2*Q0*g/2+(3*L**2*g-L*g1)*V**3/(96*S.pi**2)
    Pm=(q[2]-3*L*q[1]-3*L**2*q[0])/6-L**2*Q0*g/6+(70*L**2*g*J9-(30*L**2*g+10*L*g1)*J7+(g2-L*g1-9*L**2*g)*J5)/(48*S.pi**2)
    R=Rm-4*H[0]*R0+closed['Er_times_2pi2']/pi2
    P=Pm-4*H[0]*P0+closed['Ep_times_2pi2']/pi2
    Q=q[0]-2*H[0]*L**2*Q0-closed['Cq_times_2pi2']/pi2
    zero('complete_metric_fixed_K_Ward',D(R)-L*R+3*L*P+3*H[1]*(R0+P0))
    # Independent finite-K trace anomaly, differentiated from general physical FRW jets.
    r=S.Integer(2)
    delta_H=H[1]/L-H[0]
    delta_d=H[2]/L**2-2*H[1]/L
    delta_e=H[3]/L**3-4*H[2]/L**2+2*H[1]/L
    delta_f=H[4]/L**4-7*H[3]/L**3+10*H[2]/L**2-2*H[1]/L
    cn={5:S.Integer(0),7:5*r**2/4,9:385*r**3/64,11:-693*r**4/64,13:1155*r**5/256}
    dcn={5:r*(-18*delta_d-9*delta_e-delta_f)/16,
         7:r**2*(160*delta_H+24*delta_d-5*delta_e-delta_f)/32,
         9:7*r**3*(220*delta_H+48*delta_d+4*delta_e)/64,
         11:-231*r**4*(12*delta_H+delta_d)/64,
         13:1155*r**5*delta_H/64}
    J={n:S.expand(r**S.Rational(3-n,2)*sum((-1)**j*S.binomial((n-5)//2,j)*V**(2*j+3)/S.Integer(2*j+3) for j in range((n-5)//2+1))) for n in cn}
    dV=-H[0]*V*(1-V**2)
    delta_AK=S.expand(sum(dcn[n]*J[n]+cn[n]*S.diff(J[n],V)*dV for n in cn)/pi2)
    direct_trace=-R+3*P
    closed_trace=(D(D(Q))-2*L*D(Q)-6*L**2*Q)/2+L**2*H[1]*D(Q0)-L**2*H[0]*(D(D(Q0))+2*L*D(Q0))+L**4*delta_AK
    zero('direct_stress_vs_independent_finite_K_metric_trace',direct_trace-closed_trace)
    # Curvature forcing identity, and actual-root extra term, are exact action identities.
    z,x=S.symbols('z x',positive=True)
    delta_R=6/L**2*(H[2]+2*L*H[1]-4*L**2*H[0])
    zero('curvature_forcing_at_selected_point',g+L**2*delta_R/6)
    delta_x=S.symbols('delta_x',real=True)
    generic=2*L**2*x*H[0]+L**2*delta_x-2*L*H[1]-H[2]
    zero('generic_actual_root_curvature_decomposition',generic-(L**2*delta_x+2*L**2*(x-2)*H[0]-L**2*delta_R/6))
    zero('delta_cosmic_Hdot',D(delta_H)/L-delta_d)
    zero('delta_cosmic_Hddot',D(delta_d)/L-delta_e)
    zero('delta_cosmic_Hthird',D(delta_e)/L-delta_f)
    return {'passed':True,'identity_count':len(checks),'identities':checks,'scope':'Exact symbolic algebra only; no profile values, physical responses, quadratures or modes evaluated.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():raise RuntimeError('Refusing to overwrite development evidence')
    provenance={'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'contact_algebra_sha256':hashlib.sha256(Path(M.__file__).read_bytes()).hexdigest(),'sympy_version':S.__version__,'python_optimization':sys.flags.optimize}
    try:result=audit();result['provenance']=provenance
    except BaseException as exc:
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps({'passed':False,'exception':repr(exc),'traceback':traceback.format_exc(),'provenance':provenance},indent=2)+'\n');raise
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('Pure metric Ward/trace PASS:',result['identity_count'],'identities.')
if __name__=='__main__':main()
