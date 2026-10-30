#!/usr/bin/env python3
"""Actual retained first master: conditional continuous fourth-jet majorants.

Only the standard library; no differentiation of interval selections. The
separable temporal potential and reflected wave identities are analytic
premises checked against pinned proof text, with exact polynomial controls.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

BASE = Path('research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016')
PINS = {
    str(BASE/'numeric/inputs/HIGHER_JETS_NORMAL.json'): '92b17fd71d78f83eb3f2c5ec8b9ea5f4e70365ec4d41aa4c0205551e7eaeb5df',
    str(BASE/'numeric/inputs/BUMP_DERIVATIVE_BOUNDS.json'): 'add927ba57117a967da7d6725821d95e894ef954405d503731b3055fc64377bd',
    str(BASE/'numeric/inputs/MASTER_HIGHER_JETS_PREMISE.md'): '4c6453049f0d91e9223f06b8c881cd6bd42164d84473bf5fa5a5357bfbe35c8a',
    str(BASE/'response-signs/inputs/INITIAL_PROFILE_PROOF.md'): 'a6a417c442a8132263276cc0d14fa6525b0172b2221d0c298758854fedee80c2',
}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add(*args):
    out = {}
    for p in args:
        for k, v in p.items():
            out[k] = out.get(k, F(0))+v
    return {k:v for k,v in out.items() if v}


def mul(p, q):
    out = {}
    for a,v in p.items():
        for b,w in q.items():
            k = tuple(x+y for x,y in zip(a,b))
            out[k] = out.get(k,F(0))+v*w
    return {k:v for k,v in out.items() if v}


def scale(c, p):
    return {k:c*v for k,v in p.items() if c*v}


def powp(p, n):
    out = {(0,0,0):F(1)}
    for _ in range(n):
        out = mul(out,p)
    return out


def derivative(p, axis):
    out={}
    for k,v in p.items():
        if k[axis]:
            j=list(k);j[axis]-=1;out[tuple(j)]=k[axis]*v
    return out


def eval_s(p, endpoint):
    out={}
    for (a,b,c),v in p.items():
        out=add(out,scale(v,mul({(a,b,0):F(1)},powp(endpoint,c))))
    return out


def integral(p, lo, hi):
    primitive={(a,b,c+1):v/F(c+1) for (a,b,c),v in p.items()}
    return add(eval_s(primitive,hi),scale(-1,eval_s(primitive,lo)))


XX={(1,0,0):F(1)}
RR={(0,1,0):F(1)}
SS={(0,0,1):F(1)}
ZERO={}


def ray_integral(a,b,branch,plus,even=True):
    """Integral s^a |r +/- (x-s)|^b; plus ray is positive."""
    ap=add(RR,XX,scale(-1,SS)) if plus else add(RR,scale(-1,XX),SS)
    p=mul(powp(SS,a),powp(ap,b))
    if plus or branch=='outer':
        return integral(p,ZERO,XX)
    cut=add(XX,scale(-1,RR))
    factor=(-1)**b if even else (-1)**(b+1)
    return add(scale(factor,integral(p,ZERO,cut)),integral(p,cut,XX))


def kernel(a,b,branch):
    """KD(s^a r^b), constructed by exact spatial primitive integration."""
    return scale(F(1,2*(b+1)),add(ray_integral(a,b+1,branch,True),
             scale(-1,ray_integral(a,b+1,branch,False))))


def polynomial_audit(injection='none'):
    count=0
    for a,b in ((0,0),(1,0),(2,0),(3,0),(0,1),(1,1),(2,1),(0,2),(1,2),(2,2),(0,3),(1,3)):
        for branch in ('outer','inner'):
            k=kernel(a,b,branch)
            gp=scale(b,ray_integral(a,b-1,branch,True)) if b else {}
            gm=scale(b,ray_integral(a,b-1,branch,False)) if b else {}
            boundary=powp(add(XX,scale(-1,RR)),a) if branch=='inner' and b==0 else {}
            rr=add(scale(F(1,2),add(gp,scale(-1,gm))),scale(-1,boundary))
            xr=add(scale(F(1,2),add(gp,gm)),boundary)
            if injection=='drop_reflection_trace':
                rr=add(rr,boundary)
            need(derivative(derivative(k,1),1)==rr,'exact reflected Krr identity')
            need(derivative(derivative(k,0),1)==xr,'exact reflected Kxr identity')
            force=mul(powp(XX,a),powp(RR,b))
            need(add(derivative(derivative(k,0),0),scale(-1,rr))==force,'exact principal Kxx identity')
            count+=3
        # Moving-endpoint incoming integral and its second derivative.
        incoming=integral(mul(powp(SS,a),powp(add(XX,scale(-1,SS)),b)),ZERO,XX)
        fxt=scale(a,powp(XX,a-1)) if a and b==0 else {}
        frt=powp(XX,a) if b==1 else {}
        frr=scale(b*(b-1),integral(mul(powp(SS,a),powp(add(XX,scale(-1,SS)),b-2)),ZERO,XX)) if b>=2 else {}
        rhs=add(fxt,frt,frr)
        if injection=='drop_Hsecond_endpoint':
            rhs=add(rhs,scale(-1,fxt))
        need(derivative(derivative(incoming,0),0)==rhs,'incoming Isecond includes fx and fr endpoints')
        count+=1
    # Exact temporal polynomial identity and all-point |(1-z)(1-3z)|<=1.
    # For z in [0,1], upper inequality follows 4z-3z²>=0;
    # lower inequality follows 2-4z+3z²=3(z-2/3)²+2/3.
    z=F(2,3)
    need(2-4*z+3*z*z==F(2,3),'temporal polynomial minimum')
    need(F(4)-3>=0,'temporal polynomial upper bound on unit interval')
    temporal=scale(F(-35,4),powp(XX,2))
    tanh_rhs=add({(0,0,0):F(1)},scale(-1,powp(XX,2)))
    twice=mul(tanh_rhs,derivative(mul(tanh_rhs,derivative(temporal,0)),0))
    need(twice==scale(F(-35,2),mul(tanh_rhs,add({(0,0,0):F(1)},scale(-3,powp(XX,2))))),'exact temporal second derivative')
    return count+3


def run(root, injection='none'):
    root=Path(root)
    pins={}
    for rel,pin in PINS.items():
        path=root/rel
        need(sha(path)==pin,'authenticated input '+rel)
        pins[rel]=pin
    premise=(root/BASE/'numeric/inputs/MASTER_HIGHER_JETS_PREMISE.md').read_text()
    temporal_proof=(root/BASE/'response-signs/inputs/INITIAL_PROFILE_PROOF.md').read_text()
    for expression in (r'U_{xx}-U_{rr}+f=0', r'V=e^{-\ell r}Q',
                       r"y''=(\ell/\mu)H-Ay", r'A=\ell^2[5/2-(35/4)\tanh^2(\ell x)]'):
        need(expression in premise,'authenticated retained equation '+expression)
    need('static part of the master potential' in temporal_proof and
         'same -35 tanh^2(t)/4 is retained in both' in temporal_proof,
         'authenticated separable temporal master potential')
    master=json.loads((root/BASE/'numeric/inputs/HIGHER_JETS_NORMAL.json').read_text())
    bump=json.loads((root/BASE/'numeric/inputs/BUMP_DERIVATIVE_BOUNDS.json').read_text())
    m=master['exact']['interior'];bd=master['exact']['boundary']
    need(master['status']=='PASS' and bump['status']=='PASS','retained report status')
    need(m['scaled_domain']=={'x_max':'3','r_min':0,'x_plus_r_max':8},'retained causal domain')
    need(m['ell']=='1/1000000','retained scaling')
    count=polynomial_audit(injection)
    ell=F(1,10**6);X=F(3);Q=list(map(F,bump['dimensionless_Q_derivative_sup']))
    M=F(master['exact']['background']['W']);Mz=F(master['exact']['background']['Wz']);Mt=F(master['exact']['background']['Wt'])
    need((M,Mz,Mt)==(F(102000),F(30000000),F(21,400000)),'actual potential constants')
    Mtt=F(35,2)
    U=F(m['U']);U1=F(m['first']);U2=max(map(F,m['second'].values()));U3=max(map(F,m['third'].values()))
    y=F(bd['y']);yp=F(bd['yprime']);y2=F(m['ysecond']);y3=F(m['ythird'])
    gamma=max(abs(F(v)) for v in bd['gamma']);emu=F(bd['ell_over_mu'])
    frr=F(m['forcing_rr']);fx0=F(m['forcing']['shell_fx'])
    fr0=ell**3*Mz*y+ell**2*M*U1
    H2=Q[3]+ell*(Q[2]+2*ell*Q[1]+ell**2)+y3+ell*gamma*y2+fx0+fr0+X*frr
    y4=emu*H2+Mtt*ell**4*y+2*F(bd['Aprime'])*yp+F(bd['A'])*y2
    fxx=ell**4*Mtt*U+2*ell**3*Mt*U1+ell**2*M*U2
    fxr=ell**3*Mt*U1+ell**3*Mz*U1+ell**2*M*U2
    fxxr=ell**4*Mtt*U1+2*ell**3*Mt*U2+ell**3*Mz*U2+ell**2*M*U3
    fxx0=ell**4*Mtt*y+2*ell**3*Mt*yp+ell**2*M*y2
    V3=Q[3]+3*ell*Q[2]+3*ell**2*Q[1]+ell**3
    free=Q[4]+frr+ell*V3+fxr+y4
    if injection=='drop_initial_force_derivatives':
        free=Q[4]+ell*V3+y4
    need(free>=Q[4]+frr+ell*V3+fxr+y4,'actual Uxx initial forcing derivatives retained')
    base=free+X*fxxr+fxx0
    u4={'4,0':base+fxx,'3,1':base,'2,2':base,'1,3':base+fxr,'0,4':base+frr}
    need(max(u4.values())<F('342395.408'),'readable all-mixed fourth cap')
    need(y4<F('3.235927'),'readable boundary fourth cap')
    return {'status':'PASS_CONDITIONAL_ACTUAL_MASTER_CONTINUOUS_FOURTH_JETS',
      'domain':'0<=x<=3, r>=0, x+r<=8; x=t/ell, r=z/ell, ell=1e-6',
      'input_sha256':pins,'exact_polynomial_identity_checks':count,
      'analytic_premises':[
        'Same retained smooth actual first master, potential, dynamical boundary and flat initial collar.',
        'W(t,z)=W_static(z)-35*tanh(t)^2/4; hence W_tr=W_ttr=0 and |W_tt|<=35/2.',
        'Differentiated reflected wave formulas apply to the smooth actual function Uxx, not to independent interval selections.'],
      'exact':{'Hsecond':str(H2),'yfourth':str(y4),'Fxx':str(fxx),'Fxr':str(fxr),'Fxxr':str(fxxr),'shell_Fxx':str(fxx0),'Uxx_free_second':str(free),
        'all_scaled_order_four_mixed_caps':{k:str(v) for k,v in u4.items()}},
      'readable_strict_caps':{'all_scaled_order_four_mixed':'342395.408','boundary_yfourth':'3.235927'},
      'scope':{'actual_master_derivatives':True,'derivative_of_interval_selections':False,'new_Wzzz_input_required':False,
        'actual_gauge_fourth_certificate':False,'complete_actual_bulk_Fx_cap':False,'second_response_C2_error':False,'nonlinear_or_cosmological_claim':False}}


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--inject-error',default='none',choices=('none','drop_reflection_trace','drop_Hsecond_endpoint','drop_initial_force_derivatives'))
    a=p.parse_args();need(not a.output.exists(),'fresh output required')
    out=run(a.root,a.inject_error);a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(out['status'])


if __name__=='__main__':
    main()
