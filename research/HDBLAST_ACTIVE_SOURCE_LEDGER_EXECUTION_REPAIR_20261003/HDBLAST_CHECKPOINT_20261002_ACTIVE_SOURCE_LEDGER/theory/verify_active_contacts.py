#!/usr/bin/env python3
"""Pure finite-K contact proof; zero source samples and no saved array inputs."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import sympy as S


def require(ok,message):
    if not ok: raise RuntimeError(message)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference-module',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();require(not a.output.exists(),'Fresh receipt required')
    spec=importlib.util.spec_from_file_location('frozen_metric_contact_pure_reference',a.reference_module)
    M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M)
    data=M.inventory();closed=M.closed_quantities(data)
    L,V,H=M.L,M.V,M.H
    A=S.symbols('A',real=True)
    def D(e):return S.expand(L**2*S.diff(e,L)-L*V*(1-V**2)*S.diff(e,V)-L*V**3*S.diff(e,A)+sum(H[i+1]*S.diff(e,H[i]) for i in range(6)))
    checks=[];mutations=[]
    def zero(n,e):
        require(S.cancel(e)==0,n);checks.append(n)
    def reject(n,e):
        require(S.cancel(e)!=0,'Mutation survived: '+n);mutations.append(n)
    # Recheck all fixed-reference inventory Ward grades before using contacts.
    for i,(r,pres) in enumerate(data['grades']):
        ward=M.sub(M.plus(M.deriv(r),M.scale(M.mul((L,H[1]),pres),3)),M.mul((L,H[1]),r))
        zero('subtraction grade '+str(2*i)+' background Ward',ward[0])
        zero('subtraction grade '+str(2*i)+' metric-direction Ward',ward[1])
    pi2=2*S.pi**2
    Q0=closed['physical_baseline_Q_times_2pi2']/pi2
    R0=closed['scaled_baseline_rho_times_2pi2']/pi2
    P0=closed['scaled_baseline_p_times_2pi2']/pi2
    g,g1,g2=data['g'],data['g_prime'],data['g_second']
    J5=V**3/6;J7=(V**3/3-V**5/5)/4;J9=(V**3/3-2*V**5/5+V**7/7)/8
    qsub=g*A/(8*S.pi**2)
    qsub1=D(qsub);qsub2=D(qsub1)
    rho_local=-4*H[0]*R0+closed['Er_times_2pi2']/pi2
    p_local=-4*H[0]*P0+closed['Ep_times_2pi2']/pi2
    density_contact=(3*L**2*g-L*g1)*V**3/(96*S.pi**2)
    pressure_contact=(70*L**2*g*J9-(30*L**2*g+10*L*g1)*J7+(g2-L*g1-9*L**2*g)*J5)/(48*S.pi**2)
    # K²=2L²v²/(1-v²), so M=integral k²/(2pi²)/(2k) dk=K²/(8pi²).
    moment=L**2*V**2/(4*S.pi**2*(1-V**2))
    CR=(3*L**2*qsub-L*qsub1)/2+L**2*Q0*g/2+density_contact+rho_local
    CP=(qsub2-3*L*qsub1-3*L**2*qsub)/6-L**2*Q0*g/6+pressure_contact+p_local-moment*g/3
    FC=L*(CR-3*CP)-3*H[1]*(R0+P0)
    zero('fixed comoving K moment derivative',D(moment))
    zero('fixed comoving K renormalized baseline Ward',D(R0)-L*(R0-3*P0))
    zero('analytic finite-K local contact Ward with source work',D(CR)-FC+moment*L*g)
    zero('analytic finite-K contact integrand primitive identity',FC-(D(CR)+moment*L*g))
    Md=S.symbols('M_discrete',positive=True)
    CPmix=CP+(moment-Md)*g/3
    FCmix=L*(CR-3*CPmix)-3*H[1]*(R0+P0)
    zero('explicit mixed-target pressure shift',CPmix-CP-(moment-Md)*g/3)
    zero('explicit mixed-target integrand shift',FCmix-FC-(Md-moment)*L*g)
    zero('mixed-target contact Ward compensates discrete source moment',D(CR)-FCmix+Md*L*g)
    reject('omit pressure mode force compensation',-moment*L*g)
    reject('omit baseline work in finite-K contact integrand',3*H[1]*(R0+P0))
    reject('wrong source compensation sign',2*moment*L*g)
    reject('use frozen A instead of derivative A prime',L*V**3*g/(8*S.pi**2))
    receipt={'status':'PASS_PURE_CONTACT_ALGEBRA','identity_count':len(checks),'identities':checks,
        'mutation_count':len(mutations),'mutations_rejected':mutations,'physical_evaluations':0,
        'saved_scientific_values_loaded':0,'python_optimization':sys.flags.optimize,'sympy_version':S.__version__,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'reference_module_sha256':hashlib.sha256(a.reference_module.read_bytes()).hexdigest(),
        'scope':'Rechecked arbitrary-jet W0/W2/W4 fixed-reference Ward grades and independently assembled exact finite-K contact correction. No physical source values, finite integrals, modes or stored arrays evaluated.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','identity_count','mutation_count','physical_evaluations']}))
if __name__=='__main__':main()
