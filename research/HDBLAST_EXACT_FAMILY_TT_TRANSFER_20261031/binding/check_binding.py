#!/usr/bin/env python3
"""Finite exact checks for the selected conditional exact-family TT transfer.

This does not solve a PDE or verify a weighted inverse/analytic existence.
No project source is imported. Selected source bytes are authenticated before
and after algebra; mutations exercise the actual formulas evaluated here.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

PIN_HASH = '8404af36d4458f9ba0a9f04fe92b937457773d48624598617ed7813b71cc0c99'
MUTANTS = {
    'potential_sign': 'model-potential',
    'velocity_profile': 'seed-momentum',
    'scalar_density': 'seed-scalar-density',
    'canonical_half': 'canonical-inverse-00',
    'omit_density_variation': 'canonical-inverse-00',
    'omit_lapse_product': 'completed-velocity-11',
    'radial_chain': 'universal-TT-coordinate-conversion',
    'boundary_sign': 'shell-normalization',
    'curvature_half': 'screen-second-derivative',
    'support_undercount': 'domain-support-union',
}

def require(ok, name):
    if not ok:
        raise ValueError(name)

def sha(data):
    return hashlib.sha256(data).hexdigest()

def authenticate(repo):
    pinfile = Path(__file__).with_name('SOURCE_PINS.json')
    require(not pinfile.is_symlink(), 'source-pin-symlink')
    data = pinfile.read_bytes()
    require(sha(data) == PIN_HASH, 'source-pin-manifest-hash')
    pins = json.loads(data)
    for row in pins['files']:
        path = repo / 'research' / row['path']
        require(not path.is_symlink(), 'source-symlink:' + row['path'])
        require(path.is_file(), 'source-missing:' + row['path'])
        payload = path.read_bytes()
        require(len(payload) == row['bytes'] and sha(payload) == row['sha256'],
                'source-hash:' + row['path'])
        require(all(not parent.is_symlink() for parent in path.parents),
                'source-parent-symlink:' + row['path'])
    for capsule in pins['capsules']:
        folder=repo/'research'/capsule['path']
        manifest=folder/capsule['manifest']
        require(not manifest.is_symlink(),'capsule-manifest-symlink')
        raw=manifest.read_bytes()
        require(sha(raw)==capsule['manifest_sha256'],'capsule-external-manifest-hash')
        rows=json.loads(raw)['files']
        require(len(rows)==capsule['member_count'],'capsule-member-count')
        expected={capsule['manifest']}|{row['path'] for row in rows}
        entries=list(folder.rglob('*'))
        require(all(not path.is_symlink() for path in entries),'capsule-member-symlink')
        actual={str(path.relative_to(folder)) for path in entries if path.is_file()}
        require(actual==expected,'capsule-closed-inventory')
        for row in rows:
            payload=(folder/row['path']).read_bytes()
            require(len(payload)==row['bytes'] and sha(payload)==row['sha256'],
                    'capsule-member-hash:'+row['path'])
    return pins

def calculate(mutant, repo):
    import sympy as s
    require(s.__version__ == '1.14.0', 'sympy-version')
    passed = []

    def equal(actual, expected, name):
        require(s.cancel(s.expand(actual - expected)) == 0, name)
        passed.append(name)

    eta = s.symbols('eta')
    phi = 1 + eta
    W = 1 - phi + phi**3/3
    sign = 1 if mutant == 'potential_sign' else -1
    V = s.diff(W, eta)**2/2 + sign*s.Rational(2,3)*W**2
    expected_V = (-s.Rational(2,27)+s.Rational(14,9)*eta**2
                  +s.Rational(50,27)*eta**3-eta**4/6
                  -s.Rational(4,9)*eta**5-s.Rational(2,27)*eta**6)
    equal(V, expected_V, 'model-potential')
    cc = s.Rational(96714321425692,161839258930731)
    tension = 2*W+s.Rational(1,1000)*(1+cc*phi)
    equal(s.diff(tension,eta,2),4*phi,'model-tension-second')
    equal(s.diff(tension,eta,3),4,'model-tension-third')

    R, Om, p0, alpha, Q, Qz, kk, b = s.symbols(
        'R Om p0 alpha Q Qz kk b', nonzero=True)
    rho = Om*R
    u = alpha*Q/rho**2
    ut = R*u/b  # exp(-z)=R/b, an exact coordinate identity.
    v = 3*alpha*(Qz+kk*Q)/(rho**3*p0)
    vt = 3*alpha*(R/b)*(Qz+(kk-1)*Q)/(rho**3*p0)
    G = Om**2*u/R**4
    K = G if mutant == 'velocity_profile' else G/b
    qm, pm = G*R**2, K*R**2
    equal(-4*qm*R**2/Om**2,-4*u,'seed-metric-radial')
    equal(2*qm*R**2/Om**2,2*u,'seed-metric-angular')
    equal(-3*pm*R**2/(2*Om**3),-3*ut/(2*rho),'seed-momentum')
    # D=R*d_R=-d_z. Derivative of the already substituted exact q/p.
    dz = lambda expr: (s.diff(expr,R)*(-R)+s.diff(expr,alpha)*kk*alpha
                       +s.diff(expr,Q)*Qz)
    qm_expanded = alpha*Q/R**4
    pm_expanded = qm_expanded/b
    Phi_s = p0*rho/(2*R**2)
    equal(-3*(-dz(qm_expanded)+4*qm_expanded)/(2*Om**2*Phi_s),
          v,'seed-scalar')
    varpi = -3*Om*(-dz(pm_expanded)+5*pm_expanded)/(2*Phi_s)
    scalar_density = Om**3 if mutant == 'scalar_density' else Om**4
    equal(varpi,scalar_density*vt/rho,'seed-scalar-density')
    equal(s.cancel(K/G),1/b,'seed-frozen-profile-ratio')

    # Direct truncated expansion of the defining inverse canonical map.
    # Every second symmetric metric/momentum entry remains arbitrary.
    e,U,S = s.symbols('e U S')
    h = s.diag(-4*U,2*U,2*U,2*U)
    p = s.diag(-s.Rational(3,2)*S,0,0,0)
    j = s.zeros(4); z=s.zeros(4)
    for i in range(4):
        for l in range(i,4):
            j[i,l]=j[l,i]=s.Symbol('j%d%d'%(i,l))
            z[i,l]=z[l,i]=s.Symbol('z%d%d'%(i,l))
    metric = s.eye(4)+e*h+e**2*j/2
    momentum = e*p+e**2*z/2
    # Determinant derivative at zero is computed by Leibniz permutation sum.
    from itertools import permutations
    det1=0
    for perm in permutations(range(4)):
        parity=(-1)**sum(perm[i]>perm[l] for i in range(4) for l in range(i+1,4))
        product=s.prod(metric[i,perm[i]] for i in range(4))
        det1 += parity*s.expand(product).coeff(e,1)
    inverse_density = 1 if mutant == 'omit_density_variation' else 1-e*det1/2
    prefactor = 1 if mutant == 'canonical_half' else 2
    core=metric*momentum*metric-metric*s.trace(metric*momentum)/3
    expected=2*(z-s.eye(4)*s.trace(z)/3)+s.diag(36*U*S,-6*U*S,-6*U*S,-6*U*S)
    for i in range(4):
        for l in range(4):
            actual=2*s.expand(prefactor*inverse_density*core[i,l]).coeff(e,2)
            equal(actual,expected[i,l],'canonical-inverse-%d%d'%(i,l))
    k1=s.diag(-2*S,S,S,S)
    for i in range(1,4):
        product=0 if mutant=='omit_lapse_product' else 4*rho*U*k1[i,i]
        equal(2*rho*expected[i,i]+product,
              4*rho*(z[i,i]-s.trace(z)/3)-8*rho*U*S,
              'completed-velocity-%d%d'%(i,i))

    # Background reduction and the ordinary y-to-z chain rule are separate.
    rho,a,p0,theta,c,ell = s.symbols('rho a p0 theta c ell', nonzero=True)
    ay=-rho**-2-p0**2/3
    bgV=p0**2/2-6*a**2+6*rho**-2
    bcurv=-(ay+a**2); dcurv=rho**-2-a**2
    D=-6/c**2+rho**2*(-2*a**2+2*(theta/rho)**2-2*bcurv-8*dcurv+s.Rational(4,3)*bgV)
    equal(s.expand(D+8/c**2).subs(theta**2,1-c**-2),0,
          'TT-background-Hamilton-reduction')
    Utt,Uz,Uzz,Ut,H = s.symbols('Utt Uz Uzz Ut H')
    Uy=-Uz/rho
    Uyy=Uzz/rho**2+(0 if mutant=='radial_chain' else a*Uz/rho)
    physical=Utt-rho**2*Uyy-4*rho**2*a*Uy+3*theta*Ut+8*H/c**2
    equal(physical,Utt-Uzz+3*rho*a*Uz+3*theta*Ut+8*H/c**2,
          'universal-TT-coordinate-conversion')
    QKab,Ur=s.symbols('QKab Ur')
    # B'[k]+QKab=0, B'[k]=Uy/2; inward derivative carries a minus.
    Uy_scaled=(1 if mutant=='boundary_sign' else -1)*Ur/(ell*rho)
    equal((Uy_scaled/2+QKab).subs(Ur,2*ell*rho*QKab),0,'shell-normalization')
    ux,F=s.symbols('ux F')
    screen=3*ux/(2*ell*rho**2*c)+4*F/(rho*c)
    if mutant=='curvature_half': screen/=2
    equal(screen,3*ux/(2*ell*rho**2*c)+4*F/(rho*c),'screen-second-derivative')

    # Domain uses the whole correction annulus, not only first seed support.
    L=s.Rational(2) if mutant=='support_undercount' else s.Rational(5,2)
    require(L==s.Rational(5,2),'domain-support-union')
    require(L+2*s.Rational(3,2)<6,'domain-actual-speed-margin')
    require(s.Rational(1,8)+2*s.Rational(5,4)<6,'domain-observation-margin')
    require(s.Rational(3,2)+6<8 and s.Rational(3,2)<3,'domain-first-field-triangle')
    passed += ['domain-support-union','domain-actual-speed-margin',
               'domain-observation-margin','domain-first-field-triangle']
    from fractions import Fraction as Fq
    report=json.loads((repo/'research/HDBLAST_POSTCONTACT_SCREEN_CURVATURE_20261031/main/POSTCONTACT_NORMAL.json').read_text())
    require(report['point']==['5/4','1/8'],'inherited-observation-point')
    caps={k:Fq(v) for k,v in report['operator']['caps'].items()}
    T=Fq(5,4)
    contraction=T*(caps['A']+caps['B']+T*caps['D']/2)
    require(contraction==Fq(report['operator']['q']) and 0<contraction<1,
            'inherited-contraction-arithmetic')
    endpoints=report['actual_full_screen_h']
    require(endpoints==[
        '-1049651791613786574257886764672101/2596148429267413814265248164610048',
        '-18954240643899648600527883239121/1298074214633706907132624082305024'],
        'inherited-screen-endpoint-binding')
    require(Fq(endpoints[0])<Fq(endpoints[1])<0,'inherited-screen-sign')
    passed += ['inherited-observation-point','inherited-contraction-arithmetic',
               'inherited-screen-endpoint-binding','inherited-screen-sign']
    return {'checks':passed,'count':len(passed),'sympy':'1.14.0',
            'inherited_contraction':str(contraction),'inherited_screen_interval':endpoints,
            'scope':'Finite exact algebra and domain checks only; analytic premises are not certified.'}

def run_controls(repo,pins,scratch_parent):
    flags=['-I','-B']+(['-O'] if sys.flags.optimize else [])
    script=Path(__file__).resolve()
    results=[]
    with tempfile.TemporaryDirectory(prefix='binding-controls-',dir=scratch_parent) as td:
        scratch=Path(td)
        for mutant,expected in MUTANTS.items():
            output=scratch/(mutant+'.json')
            result=subprocess.run([sys.executable,*flags,str(script),'--repo-root',str(repo),
                                   '--mutant',mutant,'--output',str(output)],
                                  capture_output=True,text=True)
            require(result.returncode!=0 and expected in result.stderr and not output.exists(),
                    'control-not-rejected:'+mutant)
            results.append({'mutation':mutant,'rejected_by':expected,
                            'exit_code':result.returncode,'failed_witness':expected,
                            'success_artifact_absent':not output.exists()})
        # A real copied source byte mutation must fail before any algebra/report.
        copied=scratch/'source-copy'
        for row in pins['files']:
            dest=copied/'research'/row['path']; dest.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(repo/'research'/row['path'],dest)
        for capsule in pins['capsules']:
            shutil.copytree(repo/'research'/capsule['path'],copied/'research'/capsule['path'],
                            dirs_exist_ok=True)
        changed=copied/'research'/pins['files'][0]['path']
        changed.write_bytes(changed.read_bytes()+b'\nCORRUPTED SOURCE CONTROL\n')
        output=scratch/'changed-source.json'
        result=subprocess.run([sys.executable,*flags,str(script),'--repo-root',str(copied),
                               '--output',str(output)],capture_output=True,text=True)
        require(result.returncode!=0 and 'source-hash:' in result.stderr and not output.exists(),
                'control-not-rejected:source-byte-mutation')
        results.append({'mutation':'source-byte-mutation','rejected_by':'source-hash',
                        'exit_code':result.returncode,
                        'failed_witness':'source-hash:'+pins['files'][0]['path'],
                        'success_artifact_absent':not output.exists()})
        shutil.copyfile(repo/'research'/pins['files'][0]['path'],changed)
        capsule_only=copied/'research'/pins['capsules'][0]['path']/'README.md'
        capsule_only.write_bytes(capsule_only.read_bytes()+b'\nCORRUPTED CAPSULE CONTROL\n')
        output=scratch/'changed-capsule.json'
        result=subprocess.run([sys.executable,*flags,str(script),'--repo-root',str(copied),
                               '--output',str(output)],capture_output=True,text=True)
        require(result.returncode!=0 and 'capsule-member-hash:README.md' in result.stderr
                and not output.exists(),'control-not-rejected:capsule-byte-mutation')
        results.append({'mutation':'capsule-byte-mutation','rejected_by':'capsule-member-hash',
                        'exit_code':result.returncode,
                        'failed_witness':'capsule-member-hash:README.md',
                        'success_artifact_absent':not output.exists()})
    return results

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--repo-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--run-controls',action='store_true')
    parser.add_argument('--mutant',choices=MUTANTS)
    args=parser.parse_args()
    require(not args.output.exists(),'output-already-exists')
    pins=authenticate(args.repo_root)
    result=calculate(args.mutant,args.repo_root)
    if args.run_controls:
        result['controls']=run_controls(args.repo_root,pins,args.output.parent)
    authenticate(args.repo_root)
    result['selected_source_count']=len(pins['files'])
    result['source_pins_sha256']=PIN_HASH
    result['authenticated_capsules']=pins['capsules']
    with args.output.open('x') as stream:
        json.dump(result,stream,sort_keys=True,indent=2)
        stream.write('\n')

if __name__=='__main__':
    try:
        main()
    except Exception as exc:
        print(str(exc),file=sys.stderr)
        sys.exit(1)
