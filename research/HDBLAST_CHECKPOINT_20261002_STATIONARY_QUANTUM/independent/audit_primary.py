#!/usr/bin/env python3
"""Post-run aggregator using pre-frozen independent formulas and gates.

This audit selects no new parameters and runs no new radial solve. It checks
all producer rows, without reading or trusting the producer validator status.
"""
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import mpmath as mp

ROOT=Path('/workspace/HDblast/research/HDBLAST_CHECKPOINT_20261002_STATIONARY_QUANTUM')
OUT=Path(__file__).resolve().parent
PRIMARY=Path('/workspace/hdblast-research-work/stationary-primary-001/results.json')
INDEPENDENT=OUT/'ROOTS.json'
SOURCE=ROOT/'independent/independent_stationary.py'
spec=importlib.util.spec_from_file_location('frozen_independent',SOURCE)
ind=importlib.util.module_from_spec(spec)
spec.loader.exec_module(ind)
producer=json.loads(PRIMARY.read_text())
review=json.loads(INDEPENDENT.read_text())
fixed=json.loads((ROOT/'FIXED_MODEL.json').read_text())
ref=float(fixed['r_ref']);eta0=float(fixed['eta_ref'])
gates=[]
def check(name,actual,expected,absolute,relative=0):
    error=abs(actual-expected)
    threshold=absolute+relative*abs(expected)
    gates.append(dict(name=name,actual=actual,expected=expected,error=error,
                      threshold=threshold,passed=bool(error<=threshold)))

def rho_x(z,x,r):
    with mp.workdps(50):
        z,x,r=map(mp.mpf,map(str,(z,x,r)))
        u=x/z;v=r/z;nu=mp.sqrt(mp.mpf(9)/4-u)
        L=mp.digamma(mp.mpf(3)/2+nu)+mp.digamma(mp.mpf(3)/2-nu)-mp.log(v)
        if abs(nu)<mp.mpf('1e-20'):
            Pu=-mp.polygamma(2,mp.mpf(3)/2)
        else:
            Pu=(mp.polygamma(1,mp.mpf(3)/2-nu)-mp.polygamma(1,mp.mpf(3)/2+nu))/(2*nu)
        Tu=2*(u-1)*L+u*(u-2)*Pu-3*u+mp.mpf(10)/3+2*v
        return float(mp.re(z*Tu/(64*mp.pi**2)))

# Exact case coverage independently reconstructed from frozen protocol.
protocol=json.loads((ROOT/'EXPERIMENT.json').read_text())
expected=set()
for setting in protocol['resolutions']:
    gs=protocol['gamma_values'] if setting['coverage']=='all' else [0,protocol['gamma_values'][-1]]
    expected.update((setting['name'],b,g) for b in protocol['b_values'] for g in gs)
observed=[(r['resolution'],r['b'],r['gamma']) for r in producer['runs']]
gates.append(dict(name='primary_exact_case_coverage',passed=len(observed)==len(expected) and set(observed)==expected))
gates.append(dict(name='primary_no_recorded_failures',passed=not producer['failures']))
for row in producer['runs']:
    b=row['b'];z=row['H2'];factor=1+b*(row['eta']-eta0)/2
    x=ref*factor**2;xphi=b*ref*factor
    rho,Q=ind.sources(z,x,ref)
    derivatives=ind.four_derivatives(z,x,ref,b)
    label=f"{row['resolution']}/b{b}/gamma{row['gamma']}"
    check(label+'/rho',row['source']['rho'],rho,1e-22,1e-10)
    check(label+'/Q',row['source']['Q'],Q,1e-22,1e-10)
    check(label+'/current',row['source']['j'],xphi*Q/2,1e-22,1e-10)
    derivative_map={'rho_z':'rho_z','rho_eta':'rho_phi','j_z':'J_z','j_eta':'J_phi','Q_x':'Q_x','Q_z':'Q_z'}
    for primary_key,independent_key in derivative_map.items():
        check(label+'/'+primary_key,row['source'][primary_key],derivatives[independent_key],1e-20,1e-8)
    check(label+'/rho_x',row['source']['rho_x'],rho_x(z,x,ref),1e-20,1e-8)
    check(label+'/common_action_pairing',derivatives['common_action_pairing'],0.,1e-20)

# All selected endpoints at finer independently solved radial resolution.
primary_index={(r['resolution'],r['b'],r['gamma']):r for r in producer['runs']}
selected=[r for r in review['runs'] if r['setting']['rtol']==3e-14]
gates.append(dict(name='independent_six_endpoint_coverage',passed=len(selected)==6))
for row in selected:
    other=primary_index[('refined',row['b'],row['gamma'])]
    label=f"independent_endpoint/b{row['b']}/gamma{row['gamma']}"
    check(label+'/phi',row['phi'],other['phi'],2e-9)
    check(label+'/H2',row['H2'],other['H2'],0.,2e-8)
    check(label+'/ell',row['ell'],other['coordinates'][0],2e-7)
    check(label+'/y_b',row['y_b'],other['coordinates'][1],2e-7)

# Independently reconstruct the correct shell equations at all three wrong roots.
controls=[]
for row in producer['negative_controls']:
    eta=row['eta'];z=row['H2'];w=row['w'];b=row['b'];gamma=row['gamma']
    factor=1+b*(eta-eta0)/2;x=ref*factor*factor;xphi=b*ref*factor
    rho,Q=ind.sources(z,x,ref)
    W=1/3+eta*eta+eta**3/3;p=2*eta+eta*eta
    tau=.001*(1+.5975949350280132*(1+eta))+gamma*rho
    C=z+(w-p)*(w+p)/12-W*tau/9-tau*tau/36
    B=w+ind.tension_phi(eta)/2+gamma*xphi*Q/4
    detected=abs(C)>10*2e-16 or abs(B)>10*2e-13
    controls.append(dict(mutation=row['mutation'],C=C,B=B,detected=detected))
    gates.append(dict(name='correct_residual_at_wrong_root/'+row['mutation'],passed=detected))
    check('wrong_control_reconstruction/'+row['mutation']+'/C',C,row['correct_equation_residual'][0],3e-18)
    check('wrong_control_reconstruction/'+row['mutation']+'/B',B,row['correct_equation_residual'][1],3e-17)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
result=dict(status='PASS' if all(g['passed'] for g in gates) else 'FAIL',
            checked_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
            type='Post-run independent audit with frozen source formulas and predeclared source/endpoint gates; no new radial/source grid.',
            hashes={str(p):sha(p) for p in (PRIMARY,INDEPENDENT,SOURCE,Path(__file__),ROOT/'FULL_REGISTRATION.json')},
            primary_rows_checked=len(producer['runs']),independent_endpoints_checked=len(selected),
            source_derivative_rows_checked=len(producer['runs']),gates=gates,wrong_model_controls=controls)
(OUT/'PRIMARY_CROSS_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
print(result['status'],len(gates),'gates;',sum(not g['passed'] for g in gates),'failures')
for g in gates:
    if not g['passed']:print(g)
if result['status']!='PASS':raise RuntimeError('Independent raw-row primary audit failed')
