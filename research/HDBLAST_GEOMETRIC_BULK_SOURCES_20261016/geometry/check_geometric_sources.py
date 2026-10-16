#!/usr/bin/env python3
"""Exact generic-jet witnesses and manufactured geometric source checks."""
import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT.parent/'projection'))
import angular_generic as R
from geometric_sources import evaluate_sources, frame_two_jet

MUTATIONS={'omit_derivative_slot':'COVARIANT_DERIVATIVE_SLOT',
           'ricci_square_half':'CONFORMAL_RICCI_HESSIAN',
           'matter_half':'PARENT_MATTER_HESSIAN',
           'lower_time_sign':'LOWER_TIME_GAUGE_SIGN',
           'reduction_half':'REDUCTION_SYMMETRIZATION',
           'forcing_sign':'UNIT_PRINCIPAL_METRIC_SIGN'}

def need(ok,label):
    if not ok:raise ValueError(label)

def same(a,b):return not R.reduce(R.add(a,R.scale(-1,b)))
def add(*p):return R.add(*p)
def mul(a,b):return R.mul(a,b)
def scale(a,p):return R.scale(a,p)
def expected(coeff,basis):
    return {key:add(*(scale(c,base[key]) for c,base in zip(coeff,basis))) for key in basis[0]}

def fixture():
    m={k+s:F((i*7+j*3)%17-8,11) for i,k in enumerate('ABCDELM')
       for j,s in enumerate(('','y','t','yy','yt','tt'))}
    s={k:F(i+1,7) for i,k in enumerate(('s','sy','st','syy','syt','stt'))}
    b=dict(rho=F(3,2),c=F(5,4),theta=F(3,5),a=F(2,7),ay=F(-4,9),p=F(3,8),
           py=F(-2,5),Vp=F(2,3),Vpp=F(4,3),Vppp=F(7,3),ell=F(1,1000000))
    return m,s,b

def closed_scalar_hessian(s,b):
    z,zy,zt,zyy,zyt,ztt=(s[k] for k in ('s','sy','st','syy','syt','stt'))
    rho,c,a,th=(b[k] for k in ('rho','c','a','theta')); h=R.hess(R.Y);d=R.grad(R.Y)
    out={(i,j):{} for i,j in product(range(6),repeat=2)}
    out[0,0]=scale(zyy,R.Y);out[0,1]=out[1,0]=scale((zyt-a*zt)/rho,R.Y)
    out[1,1]=scale(ztt/rho**2-a*zy,R.Y)
    for i in range(4):
        out[0,i+2]=out[i+2,0]=scale((zy-a*z)/(rho*c),d[i,])
        out[1,i+2]=out[i+2,1]=scale((zt-th*z)/(rho*rho*c),d[i,])
        for j in range(4):
            out[i+2,j+2]=add(scale(a*zy-th*zt/rho**2,mul(R.Y,R.P[i,j])),scale(z/(rho*rho*c*c),h[i,j]))
    return out

def run(mutation=None):
    for name,pin in json.loads((ROOT/'INPUT_PINS.json').read_text()).items():
        raw=(ROOT/name).read_bytes()
        need(len(raw)==pin['bytes'] and hashlib.sha256(raw).hexdigest()==pin['sha256'],'INPUT_AUTHENTICATION_'+name)
    m,s,b=fixture();I=range(6);sig=(1,-1,1,1,1,1);checks=0
    def eq(a,b0,label):
        nonlocal checks
        need(same(a,b0),label);checks+=1
    # The scalar Hessian has an independently written six-sector formula.
    cm={key+suffix:(-1 if key=='D' else 1)*s['s'+suffix] if key in ('A','D','L') else F(0)
        for key in 'ABCDELM' for suffix in ('','y','t','yy','yt','tt')}
    cj=frame_two_jet(cm,s,b,mutation=mutation)
    h=closed_scalar_hessian(s,b)
    for k in h:eq(cj['ddeta'][k].p,h[k],'SCALAR_HESSIAN_FRAME_CONNECTION')
    for d,a,b0 in product(I,repeat=3):
        eq(cj['dq'][d,a,b0].p,mul(cj['deta'][d,].p,cj['g0'][a,b0].p),'CONFORMAL_FIRST_COVARIANT_JET')
    for e,d,a,b0 in product(I,repeat=4):
        eq(cj['ddq'][e,d,a,b0].p,mul(h[e,d],cj['g0'][a,b0].p),'COVARIANT_DERIVATIVE_SLOT')
    # For q=f g in dimension5, independent conformal transformation law gives
    # Ric''=3 f Hess(f)+(9/2)df df+g(f Delta f-(1/2)|df|²).
    cr=evaluate_sources(cm,s,b,mutation=mutation)
    f=scale(s['s'],R.Y)
    lap=add(*(scale(sig[i],h[i,i]) for i in I))
    norm=add(*(scale(sig[i],mul(cj['deta'][i,].p,cj['deta'][i,].p)) for i in I))
    for a,b0 in product(I,repeat=2):
        rhs=add(scale(3,mul(f,h[a,b0])),scale(F(9,2),mul(cj['deta'][a,].p,cj['deta'][b0,].p)),
                mul(cj['g0'][a,b0].p,add(mul(f,lap),scale(F(-1,2),norm))))
        eq(cr['QRic'][a,b0],rhs,'CONFORMAL_RICCI_HESSIAN')
    r=evaluate_sources(m,s,b,mutation=mutation)
    spec=importlib.util.spec_from_file_location('frozen_source_adapters',ROOT/'inputs'/'source_adapters.py')
    parent=importlib.util.module_from_spec(spec);spec.loader.exec_module(parent)
    pq=parent.gauge_source(m,b);pm=parent.einstein_matter_source(m,s,b)
    gauge={}
    for i,key in ((0,'y'),(1,'0')):gauge[i,]=expected(pq[key],R.SCALARS)[()]
    for (i,),p in expected(pq['a'],R.GRADS).items():gauge[i+2,]=p
    if mutation=='lower_time_sign':r['QC'][1,]=scale(-1,r['QC'][1,])
    for key in gauge:eq(r['QC'][key],gauge[key],'LOWER_TIME_GAUGE_SIGN' if key==(1,) else 'PARENT_GAUGE_HESSIAN')
    eq(r['Qphi'],expected(parent.scalar_bulk_source(m,s,b),R.SCALARS)[()],'PARENT_SCALAR_HESSIAN')
    for pair,key in (((0,0),'yy'),((0,1),'y0'),((1,1),'00')):
        eq(r['QE_matter'][pair],expected(pm[key],R.SCALARS)[()],'PARENT_MATTER_HESSIAN')
    for fixed,key in ((0,'ya'),(1,'0a')):
        for (i,),p in expected(pm[key],R.GRADS).items():eq(r['QE_matter'][fixed,i+2],p,'PARENT_MATTER_HESSIAN')
    for (i,j),p in expected(pm['ab'],R.TENSORS).items():eq(r['QE_matter'][i+2,j+2],p,'PARENT_MATTER_HESSIAN')
    for name in ('QRic','QE','QEred','QEred_general','QE_matter'):
        for a,b0 in product(I,repeat=2):eq(r[name][a,b0],r[name][b0,a],'SOURCE_SYMMETRY_'+name)
    for a,b0 in product(I,repeat=2):
        gradient=scale(F(1,2),add(r['dQC'][a,b0],r['dQC'][b0,a]))
        reduced=r['QEred'][a,b0]
        if mutation=='reduction_half':reduced=add(r['QE'][a,b0],scale(F(-1,2),gradient))
        eq(reduced,add(r['QE'][a,b0],scale(-1,gradient)),'REDUCTION_SYMMETRIZATION')
        forcing=r['metric_forcing'][a,b0]
        if mutation=='forcing_sign':forcing=scale(-1,forcing)
        eq(forcing,scale(-2*b['ell']**2*b['rho']**2,reduced),'UNIT_PRINCIPAL_METRIC_SIGN')
    need(any(r['C1'].values()),'GENERIC_WITNESS_MUST_BE_OFF_HARMONIC')
    need(any(not same(r['QC'][k],r['QC_background'][k]) for k in r['QC']),'GAUGE_LOWERING_DISTINCTION')
    need(any(not same(r['QEred'][k],r['QEred_general'][k]) for k in r['QEred']),'REDUCTION_OFF_HARMONIC_DISTINCTION')
    eq(r['scalar_forcing'],scale(b['ell']**2*b['rho']**2,r['Qphi']),'UNIT_PRINCIPAL_SCALAR_SIGN')
    if mutation is None:
        # Public integer inputs must remain in the exact rational ring.
        # rho=2 forces an inverse that native Python 1/2 would turn into float.
        im={key:value.numerator for key,value in m.items()}
        iss={key:value.numerator for key,value in s.items()}
        ib=dict(rho=2,c=1,theta=0,a=1,ay=-1,p=2,py=3,Vp=1,Vpp=2,Vppp=3,ell=1)
        ir=evaluate_sources(im,iss,ib)
        fr=evaluate_sources({k:F(v) for k,v in im.items()},
                            {k:F(v) for k,v in iss.items()},
                            {k:F(v) for k,v in ib.items()})
        def exact_tree(left,right):
            if isinstance(left,dict):
                return isinstance(right,dict) and left.keys()==right.keys() and all(exact_tree(v,right[k]) for k,v in left.items())
            return isinstance(left,(int,F)) and isinstance(right,(int,F)) and left==right
        need(exact_tree(ir,fr),'INTEGER_VERSUS_FRACTION_SOURCE_EXACTNESS')
        try:frame_two_jet(im,iss,{**ib,'rho':2.0})
        except ValueError as exc:need(str(exc)=='FLOAT_INPUT_REQUIRES_EXACT_COEFFICIENT_RING','FLOAT_REJECTION_REASON')
        else:raise ValueError('FLOAT_INPUT_NOT_REJECTED')
    encoded={name:{','.join(map(str,k)):{','.join(map(str,e)):str(v) for e,v in sorted(p.items())}
                        for k,p in sorted(r[name].items())} for name in ('QRic','QE','QEred','QEred_general','QC','QC_background')}
    digest=hashlib.sha256(json.dumps(encoded,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return {'status':'pass','exact_polynomial_comparisons':checks,'generic_source_digest_sha256':digest,
            'manufactured':'five-dimensional conformal metric Hessian, all sphere points',
            'witness_scope':'exact generic rational jets; universal independent-jet coefficient identities are checked in projection',
            'harmonic_specialization_explicit':True,'full_lowered_and_background_lowered_QC_distinguished':True,
            'integer_fraction_full_source_equality':True,'raw_float_inputs_rejected':True}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);ap.add_argument('--mutation',choices=MUTATIONS)
    ap.add_argument('--run-controls',action='store_true');args=ap.parse_args()
    if args.output and args.output.exists():raise ValueError('OUTPUT_EXISTS')
    report=run(args.mutation)
    if args.run_controls:
        controls=[]
        for mutation,reason in MUTATIONS.items():
            cmd=[sys.executable,'-I','-S','-B']+(['-O'] if not __debug__ else [])+[str(Path(__file__).resolve()),'--mutation',mutation]
            cp=subprocess.run(cmd,text=True,capture_output=True)
            need(cp.returncode!=0 and reason in cp.stderr,'MUTATION_REASON_'+mutation)
            controls.append({'mutation':mutation,'required_diagnostic':reason,'rejected':True})
        report['controls']=controls
    raw=json.dumps(report,sort_keys=True,indent=2)+'\n'
    if args.output:args.output.write_text(raw)
    else:print(raw,end='')

if __name__=='__main__':main()
