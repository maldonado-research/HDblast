#!/usr/bin/env python3
"""Independent exact algebra, not PDE evolution or a physical sign claim."""
import argparse
from functools import lru_cache
import hashlib
import json
from math import factorial
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
PINS_SHA256 = '0da81d7eb58b5ef408891447fbddf1843d3bf2abf9038ed2af97a5bc26c85c49'
CHECKS = []
FIXTURES = {}


def eq(label, left, right=0):
    residue = s.factor(s.expand(left-right))
    if residue != 0:
        raise ValueError('IDENTITY_FAILED:'+label+':'+str(residue))
    CHECKS.append(label)


def sphere_average(expr, coords):
    """Uniform S3 moments, derived by normalized Gaussian radial moments."""
    out = 0
    for powers, coefficient in s.Poly(s.expand(expr), *coords).terms():
        if any(p % 2 for p in powers):
            continue
        halves = [p//2 for p in powers]
        moment = s.Integer(1)
        for h in halves:
            moment *= s.Rational(factorial(2*h), 2**h*factorial(h))
        for j in range(sum(halves)):
            moment /= 4+2*j
        out += coefficient*moment
    return s.factor(out)


def gaussian_contraction(mutation):
    # Arbitrary symmetric K and its orthonormal-frame time derivative J.
    symbols = s.symbols('k00 k01 k02 k11 k12 k22 j00 j01 j02 j11 j12 j22')
    def matrix(v):
        return s.Matrix([[v[0],v[1],v[2]],[v[1],v[3],v[4]],[v[2],v[4],v[5]]])
    K, J = matrix(symbols[:6]), matrix(symbols[6:])
    ric3 = matrix(s.symbols('r00 r01 r02 r11 r12 r22'))
    scalar3 = s.trace(ric3)
    metric = [-1,1,1,1]
    delta = lambda a,b: s.Integer(a == b)
    def r3(a,b,c,d):
        return (delta(a,c)*ric3[b,d]-delta(a,d)*ric3[b,c]
                -delta(b,c)*ric3[a,d]+delta(b,d)*ric3[a,c]
                -scalar3*(delta(a,c)*delta(b,d)-delta(a,d)*delta(b,c))/2)
    @lru_cache(None)
    def r4(a,b,c,d):
        if a == b or c == d:
            return s.Integer(0)
        if a and b and c and d:
            i,j,k,l = a-1,b-1,c-1,d-1
            return r3(i,j,k,l)+K[i,k]*K[j,l]-K[i,l]*K[j,k]
        if (a == 0) != (b == 0) and (c == 0) != (d == 0):
            sign = (1 if a == 0 else -1)*(1 if c == 0 else -1)
            i = (b if a == 0 else a)-1
            j = (d if c == 0 else c)-1
            return -sign*(J+K*K)[i,j]
        return s.Integer(0)  # Codazzi slots do not enter these contractions.
    ric4 = s.Matrix(4,4,lambda b,d: sum(metric[a]*r4(a,b,a,d) for a in range(4)))
    scalar4 = s.expand(sum(metric[a]*ric4[a,a] for a in range(4)))
    theta, theta_dot = s.trace(K), s.trace(J)
    shear2 = s.trace(K*K)-theta**2/3
    eq('Gaussian Ric00 direct contraction',ric4[0,0],-theta_dot-theta**2/3-shear2)
    eq('Gaussian scalar direct contraction',scalar4,scalar3+2*theta_dot+4*theta**2/3+shear2)
    eq('Gaussian G00 direct contraction',ric4[0,0]+scalar4/2,scalar3/2+theta**2/3-shear2/2)
    # Independent normalized-integral derivative, with variance symbolic.
    H, Hp, var, sh, r = s.symbols('H Hp variance shear r3')
    mean_theta_dot = 3*Hp + (var if mutation == 'commutator' else -var)
    mean_theta2 = 9*H**2+var
    mean_r4 = r+2*mean_theta_dot+4*mean_theta2/3+sh
    q = 2*var/3-sh
    eq('Averaged scalar commutator',mean_r4,r+6*Hp+12*H**2-q)


def sphere_and_curvature(mutation):
    n=s.Matrix(s.symbols('n0:4'))
    Q=s.diag(1,-1,0,0)
    Y=(n.T*Q*n)[0]
    proj=s.eye(4)-n*n.T
    grad=2*proj*Q*n
    hess=2*proj*Q*proj-2*Y*proj
    tf=hess+s.Rational(8,3)*Y*proj
    avg=lambda p: sphere_average(p,list(n))
    y2=avg(Y**2)
    grad2=avg((grad.T*grad)[0])
    hess2=avg(s.trace(hess*hess))
    tf2=avg(s.trace(tf*tf))
    eq('Ambient S3 Y square',y2,s.Rational(1,6))
    eq('Ambient S3 gradient square',grad2,s.Rational(4,3))
    eq('Ambient S3 Hessian square',hess2,8)
    eq('Ambient S3 tracefree Hessian square',tf2,s.Rational(32 if mutation == 'shear' else 40,9))
    E,A,T,L=s.symbols('E A T L',nonzero=True)
    theta1=(3*T+8*L)*Y
    var2=2*avg(theta1**2)
    sh2=2*L**2*tf2
    eq('First-only expansion variance',var2,(3*T+8*L)**2/3)
    eq('First-only shear derivative',sh2,s.Rational(80,9)*L**2)
    eq('First-only Buchert Q',2*var2/3-sh2,s.Rational(2,9)*((3*T+8*L)**2-40*L**2))
    # DC=0 for arbitrary first trace; hence no h2 can contribute.
    trace=s.symbols('trace')
    dV=trace/2
    dI=trace/A**2
    dmean=dI-6*dV/A**2
    dref=-4*dV/A**2
    eq('Round normalized curvature is critical',dmean,dref)
    # Two curves with the same first tangent and distinct second tangents.
    # h_exp=e^(2e E Y)h0; h_lin=(1+2e E Y)h0.
    for name, volume_factor, conformal_factor in [('exponential',9,1),('linear',3,-1)]:
        volume2=volume_factor*E**2*y2
        integral2=E**2*(6*conformal_factor*y2+4*grad2)/A**2
        mean2=integral2-6*volume2/A**2
        ref2=-4*volume2/A**2
        coeff=s.Rational(8 if mutation == 'curvature' else 10,3)
        eq('Curvature defect '+name+' second jet',mean2-ref2,coeff*E**2/A**2)
    # Exact finite counterexample to <R3>=6/Rv²: h=(1+aY)^4 h0.
    # The positive conformal factor obeys 1+aY>=9/10. The Yamabe
    # transformation gives N=<6f²+8|Df|²>, D=<f^6>. C>0 iff N³>216D.
    amplitude=s.Rational(1,10)
    conformal=1+amplitude*Y
    numerator=avg(6*conformal**2+8*amplitude**2*(grad.T*grad)[0])
    volume=avg(conformal**6)
    gap=s.factor(numerator**3-216*volume)
    if not gap>0:
        raise ValueError('IDENTITY_FAILED:Exact finite non-round curvature defect')
    CHECKS.append('Exact finite non-round curvature defect is positive')
    FIXTURES['conformal_phi4_a_1_over_10']={
        'I_over_V0_times_A2':str(numerator),'V_over_V0':str(volume),
        'positive_cube_gap':str(gap),'positive_factor_lower':'9/10'}


def clock_and_trace(mutation):
    u,up,H,Hp,P,b,bp,d,dp,J=s.symbols('u up H Hp P b bp d dp J')
    rho,A=s.symbols('rho A',nonzero=True)
    clock=P+rho*b
    graph_isotropic=u+H*rho*b-H*clock
    graph_hessian=d
    transport=-d-J
    if mutation == 'transport':
        transport=0
    eq('First clock cancels time gauge',graph_isotropic,u-H*P)
    eq('Angular transport cancels spatial gauge',graph_hessian+transport,-J)
    # Exact scalar trace coefficients in the current shell boundary laws.
    p,a,sigma2,Vp=s.symbols('p a sigma2 Vp',nonzero=True)
    sigmap=-2*p
    sigma=6*a
    Fp=-sigmap*sigma2/2+2*Vp+2*sigma*sigmap/3
    B=Vp/p-4*a+sigma2/2
    eq('Current nonzero scalar mean coefficient',Fp,0 if mutation == 'monopole' else 2*p*B)
    # Lorentzian conformal variation, with Laplacian eigenvalue -8.
    upp=s.symbols('upp')
    box=-upp-3*H*up-8*u/A**2
    dr=-2*u*12/rho**2-6*box
    eq('Actual first intrinsic trace ODE',dr,6*(upp+3*H*up+(8/A**2-4/rho**2)*u))
    # The averaged second trace equation, including measure-averaged phi2.
    v,vp,Fpp,M,v2,v2p,v2pp,Q2,E=s.symbols('v vp Fpp M v2 v2p v2pp Q2 E')
    r42=Fp*M+Fpp*v**2/6-vp**2/3+8*v**2/(3*A**2)
    r32=-4*v2/A**2+10*E**2/(3*A**2)
    rhs=Fp*M/2+Fpp*v**2/12-vp**2/6+4*v**2/(3*A**2)-5*E**2/(3*A**2)+Q2/2
    eq('Second averaged scalar-volume balance',(r42-r32+Q2)/2,rhs+2*v2/A**2)
    eq('Second HD derivative normalization',6*(v2pp/3+4*H*v2p/3),2*v2pp+8*H*v2p)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--repo-root',type=Path,default=Path.cwd())
    destination=parser.add_mutually_exclusive_group()
    destination.add_argument('--output',type=Path)
    destination.add_argument('--output-dir',type=Path)
    parser.add_argument('--mutation',default='none')
    parser.add_argument('--run-controls',action='store_true')
    args=parser.parse_args()
    pins_path=HERE/'SOURCE_PINS.json'
    if hashlib.sha256(pins_path.read_bytes()).hexdigest()!=PINS_SHA256:
        raise ValueError('SOURCE_PINS_HASH_FAILED')
    pins=json.loads(pins_path.read_text())
    repo_root=args.repo_root.resolve(strict=True)
    for entry in pins['inputs']:
        base=repo_root if entry['base']=='repo' else HERE
        path=base/entry['path']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=entry['sha256']:
            raise ValueError('SOURCE_HASH_FAILED:'+entry['base']+':'+entry['path'])
    commit=pins['source_commit']
    # Authenticate the complete declared source set before loading SymPy.
    global s
    import sympy as s
    gaussian_contraction(args.mutation)
    sphere_and_curvature(args.mutation)
    clock_and_trace(args.mutation)
    controls=[]
    strict_control_rows=[]
    if args.run_controls:
        witnesses={
            'commutator':'ValueError: IDENTITY_FAILED:Averaged scalar commutator:4*variance',
            'shear':'ValueError: IDENTITY_FAILED:Ambient S3 tracefree Hessian square:8/9',
            'curvature':'ValueError: IDENTITY_FAILED:Curvature defect exponential second jet:2*E**2/(3*A**2)',
            'transport':'ValueError: IDENTITY_FAILED:Angular transport cancels spatial gauge:J + d',
            'monopole':'ValueError: IDENTITY_FAILED:Current nonzero scalar mean coefficient:2*Vp - 8*a*p + p*sigma2',
        }
        with tempfile.TemporaryDirectory(prefix='shell-expansion-controls-') as temp_dir:
            for name,expected in witnesses.items():
                child_output=Path(temp_dir)/(name+'.json')
                absent_before=not child_output.exists() and not child_output.is_symlink()
                if not absent_before:
                    raise ValueError('CONTROL_OUTPUT_NOT_FRESH:'+name)
                cmd=[sys.executable,'-I','-B']
                if sys.flags.optimize:
                    cmd+=['-O']
                cmd += [str(Path(__file__).resolve()),'--repo-root',str(repo_root),
                        '--mutation',name,'--output',str(child_output)]
                proc=subprocess.run(cmd,capture_output=True,text=True)
                lines=[line.strip() for line in proc.stderr.splitlines() if line.strip()]
                observed=lines[-1] if lines else ''
                absent_after=not child_output.exists() and not child_output.is_symlink()
                stdout_absent=proc.stdout.strip()==''
                if proc.returncode==0:
                    raise ValueError('CONTROL_ZERO_EXIT:'+name)
                if observed!=expected:
                    raise ValueError('CONTROL_WRONG_FINAL_WITNESS:'+name+':'+observed)
                if not stdout_absent or not absent_after:
                    raise ValueError('CONTROL_SUCCESS_OUTPUT_PRESENT:'+name)
                strict_control_rows.append({
                    'case':name,'expected_final_witness':expected,
                    'observed_final_witness':observed,'subprocess_exit_code':proc.returncode,
                    'nonzero_exit':True,'success_stdout_absent':stdout_absent,
                    'success_json_absent':stdout_absent,
                    'output_artifact_absent_before':absent_before,
                    'output_artifact_absent_after':absent_after})
                controls.append(name)
    report={'status':'FINITE_EXACT_IDENTITIES_PASS','checks':CHECKS,'exact_fixtures':FIXTURES,'rejected_formula_controls':controls,
            'strict_formula_control_rows':strict_control_rows,
            'source_commit':commit,'source_pins_sha256':hashlib.sha256((HERE/'SOURCE_PINS.json').read_bytes()).hexdigest(),
            'scope':'Exact geometry/clock algebra; no PDE solve or evaluated physical sign.'}
    raw=json.dumps(report,indent=2,sort_keys=True)+'\n'
    output=args.output
    if args.output_dir:
        args.output_dir.mkdir(parents=True,exist_ok=True)
        output=args.output_dir/'SHELL_EXPANSION.json'
    if output:
        output.write_text(raw)
    else:
        print(raw,end='')


if __name__=='__main__':
    main()
