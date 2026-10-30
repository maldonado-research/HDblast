#!/usr/bin/env python3
"""Exact first master/dynamic-shell/physical-junction implementation binding.

Uses no scientific dependencies and imports no historical project modules.
Source AST expressions are used only after complete manifest authentication.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F
import subprocess
import sys
import tempfile

CAPSULE_PINS = {
    'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016':
        '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2',
    'HDBLAST_INTERIOR_MASTER_JETS_20261015':
        'b043a13a228d11296fbbde292aa15f7a024deab9ff915bdce452c52c73745bbb',
    'HDBLAST_REGULAR_SHELL_SPECTRA_20261013':
        'd1d47c66cef1ed79979c22ae7145dc66d35fa35be69dba900f9912ac01f31e2f',
}

def need(ok, message):
    if not ok:
        raise ValueError(message)

class P:
    """Collected sparse Laurent polynomial, rational coefficients."""
    def __init__(self, value=0):
        if isinstance(value, P):
            self.terms = dict(value.terms)
        elif isinstance(value, dict):
            self.terms = {k:F(v) for k,v in value.items() if v}
        else:
            self.terms = {():F(value)} if value else {}
    @staticmethod
    def symbol(name):
        return P({((name, 1),):F(1)})
    @staticmethod
    def cv(value):
        return value if isinstance(value,P) else P(value)
    def __add__(self, other):
        out = dict(self.terms)
        for k,v in P.cv(other).terms.items():
            out[k] = out.get(k,F(0))+v
        return P(out)
    __radd__ = __add__
    def __neg__(self):
        return P({k:-v for k,v in self.terms.items()})
    def __sub__(self, other):
        return self + -P.cv(other)
    def __rsub__(self, other):
        return P.cv(other) + -self
    def __mul__(self, other):
        out = {}
        for k,v in self.terms.items():
            for j,w in P.cv(other).terms.items():
                d = dict(k)
                for n,e in j:
                    d[n] = d.get(n,0)+e
                key = tuple(sorted((n,e) for n,e in d.items() if e))
                out[key] = out.get(key,F(0))+v*w
        return P(out)
    __rmul__ = __mul__
    def __pow__(self, n):
        if n < 0:
            need(len(self.terms)==1, 'only monomial denominators are allowed')
            (key,v), = self.terms.items()
            return P({tuple((name,e*n) for name,e in key):v**n})
        out=P(1)
        for _ in range(n):
            out=out*self
        return out
    def __truediv__(self, other):
        return self*P.cv(other)**-1
    def __rtruediv__(self, other):
        return P.cv(other)*self**-1
    def __eq__(self, other):
        return self.terms==P.cv(other).terms
    def diff(self, name):
        out={}
        for key,v in self.terms.items():
            d=dict(key); e=d.get(name,0)
            if e:
                d[name]=e-1
                key2=tuple(sorted((n,p) for n,p in d.items() if p))
                out[key2]=v*e
        return P(out)
    def substitute(self, rules):
        out=P()
        for key,v in self.terms.items():
            q=P(v)
            for n,e in key:
                q=q*rules.get(n,P.symbol(n))**e
            out=out+q
        return out
    def derivative(self, rules):
        names={n for key in self.terms for n,e in key}
        return sum((self.diff(n)*rules.get(n,P()) for n in names),P())
    def dr(self):
        return self.derivative(DR)
    def dx(self):
        return self.derivative(DX)
    def encode(self):
        return [[[list(item) for item in key],str(v)] for key,v in sorted(self.terms.items())]

NAMES = ('rho','p','a','h','V2','s','alpha','vartheta','H','c','ell',
         'W','Wr','Wrr','Wx','Wxx','f','fr','fx','frr','frx','fxx',
         'b','br','bx','brr','brx','bxx','d','dr','dx','drr','drx','dxx')
G={n:P.symbol(n) for n in NAMES}
rho,p,a,h,V2,s,alpha,vartheta,H,c,ell=[G[n] for n in NAMES[:11]]
DR={
    'rho':-ell*rho**2*a,
    'p':-ell*rho*p*(h-4*a),
    'a':ell*(rho**-1+rho*p**2/3),
    'h':-ell*rho*(V2-h**2+4*a*h),
    'alpha':ell*rho*(F(5,2)*a-h)*alpha,
    'W':G['Wr'], 'Wr':G['Wrr'],
}
DX={'H':ell*(1-H**2), 'c':ell*H*c,
    'vartheta':-F(3,2)*ell*H*vartheta,
    'W':G['Wx'], 'Wx':G['Wxx']}
for name in ('f','b','d'):
    DR[name]=G[name+'r'];DR[name+'r']=G[name+'rr'];DR[name+'x']=G[name+'rx']
    DX[name]=G[name+'x'];DX[name+'x']=G[name+'xx'];DX[name+'r']=G[name+'rx']

def dy(value, reverse=False):
    return (-1 if not reverse else 1)*P.cv(value).dr()/(ell*rho)

def authenticate(repo):
    receipt={}
    for capsule,pin in CAPSULE_PINS.items():
        root=repo/'research'/capsule
        need(not root.is_symlink(),'source capsule may not be a symlink')
        raw=(root/'MANIFEST.json').read_bytes()
        need(hashlib.sha256(raw).hexdigest()==pin,'source manifest: '+capsule)
        members=json.loads(raw)['files']
        seen=set()
        for row in members:
            rel=Path(row['path'])
            need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,
                 'safe unique manifest path')
            seen.add(str(rel));content=(root/rel).read_bytes()
            need(len(content)==row['bytes'] and hashlib.sha256(content).hexdigest()==row['sha256'],
                 'source member: '+capsule+'/'+str(rel))
        entries=list(root.rglob('*'))
        need(not any(path.is_symlink() for path in entries),'source symlink rejected: '+capsule)
        actual={str(path.relative_to(root)) for path in entries if path.is_file()}
        need(actual==seen|{'MANIFEST.json'},'closed source membership: '+capsule)
        receipt[capsule]={'manifest_sha256':pin,'authenticated_members':len(members)}
    return receipt

def function_ast(path,name):
    tree=ast.parse(path.read_text())
    nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name]
    need(len(nodes)==1,'unique source function '+name)
    return nodes[0]

def run(repo, mutation='none'):
    before=authenticate(repo)
    checks=[];ast_receipts={}
    def same(label,left,right):
        residual=P.cv(left)-P.cv(right)
        need(not residual.terms,'IDENTITY FAILED: '+label+'; residual='+json.dumps(residual.encode()))
        checks.append(label)
    geo=repo/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016'
    master=repo/'research/HDBLAST_INTERIOR_MASTER_JETS_20261015'
    spectra=repo/'research/HDBLAST_REGULAR_SHELL_SPECTRA_20261013'
    frame_path=geo/'numeric/actual_first_jets.py'
    frame_node=function_ast(frame_path,'frame_amplitudes')
    env={'__builtins__':{},'F':F,'ELL':ell}
    exec(compile(ast.Module(body=[frame_node],type_ignores=[]),str(frame_path),'exec'),env)
    ast_receipts['frame_amplitudes']=hashlib.sha256(ast.dump(frame_node,include_attributes=False).encode()).hexdigest()
    bg_path=geo/'numeric/interval_core.py'
    bg_node=function_ast(bg_path,'background')
    assigns=[];alpha_rhs=None
    for node in ast.walk(bg_node):
        if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
            if node.targets[0].id in ('k','C','B'):
                assigns.append(node)
            if node.targets[0].id=='rhs' and isinstance(node.value,ast.Dict):
                for key,value in zip(node.value.keys,node.value.values):
                    if isinstance(key,ast.Constant) and key.value=='alpha':
                        alpha_rhs=value
    need(len(assigns)==3 and alpha_rhs is not None,'selected source coefficient closure')
    coeff_env={'__builtins__':{},'F':F,'rho':rho,'p':p,'a':a,'h':h,'alpha':alpha}
    exec(compile(ast.Module(body=assigns,type_ignores=[]),str(bg_path),'exec'),coeff_env)
    source_R,source_B,source_k=[coeff_env[n] for n in ('C','B','k')]
    source_alpha_z=eval(compile(ast.Expression(alpha_rhs),str(bg_path),'eval'),coeff_env)
    ast_receipts['background_selected_coefficients']=hashlib.sha256(
        ast.dump(ast.Module(body=assigns,type_ignores=[]),include_attributes=False).encode()).hexdigest()
    ast_receipts['background_alpha_z']=hashlib.sha256(ast.dump(alpha_rhs,include_attributes=False).encode()).hexdigest()
    same('actual R coefficient',source_R,alpha/rho**2)
    same('actual scalar B coefficient',source_B,3*alpha/(rho**3*p))
    same('actual k coefficient',source_k,rho*(F(5,2)*a-h))
    same('actual alpha inward derivative',source_alpha_z,alpha.dr()/ell)
    # Read retained analytic laws explicitly. They are analytic inputs, not an executed trace solver.
    trace_note=(master/'jets/HIGHER_JETS.md').read_text()
    need('H=U_r(x,0)-ell*gamma*y' in trace_note,
         'retained dynamic trace source')
    need(r'5/2-(35/4)\tanh^2(\ell x)' in trace_note,
         'retained dynamic trace time coefficient')
    need("y''=(\\ell/\\mu)H-Ay" in trace_note,'retained second trace equation')
    physical=(spectra/'operator/PHYSICAL_OPERATOR_AUDIT_20261013.md').read_text()
    need("\\mathcal J_\\phi=X'-p\\widetilde N+\\sigma''X/2+pB\\xi" in physical,
         'retained physical scalar boundary source')
    need("(\\Psi'-a\\widetilde N+pX/3)Y\\delta^\\mu{}_{\\nu}" in physical,
         'retained physical metric boundary source')
    # Source pin authenticates the full arbitrary-function physical boundary formulas.
    law_note=(master/'inputs/reference/PHYSICAL_SHELL_RESPONSE.md').read_text()
    need(r'\mu=(\rho_bB_b)^{-1}>0' in law_note and r'B_bP_y=\rho_b^{-2}(\Box_\gamma+4)P_b' in law_note,
         'retained physical-to-dynamic boundary law')
    # Current actual reconstruction with completely free longitudinal master 2-jet.
    W=G['W'];Wr=G['Wr'];Wrr=G['Wrr'];Wxx=G['Wxx']
    amps,scalar=env['frame_amplitudes'](W,[P(),P(),P()],rho,p,a,source_R,source_B,source_k,c,H,vartheta)
    u=source_R*vartheta*W
    Pmaster=rho**2*u
    chi=scalar
    if mutation=='reverse_inward_y':
        derivative=lambda q:dy(q,True)
    else:
        derivative=dy
    same('longitudinal metric lapse AST',amps['A'],-4*u)
    same('longitudinal metric time AST',amps['D'],-2*u)
    same('longitudinal metric sphere AST',amps['L'],2*u)
    same('longitudinal master momentum reconstruction',chi,
         (-3 if mutation!='wrong_scalar_reconstruction' else 3)*derivative(Pmaster)/(rho**2*p))
    # Unit normal is N^-1=1+2 epsilon*uY at this longitudinal first variation.
    # K1=1/2*m_y−1/2*h_yy*a*gamma; subtract sigma1*chi*gamma/6+sigma0*m/6.
    sigma0=6*a; sigma1=-2*p
    metric_from_affine=derivative(u)+2*a*u-sigma1*chi/(6 if mutation!='double_tension' else 3)
    same('one-copy physical metric shell residual',metric_from_affine,0)
    # The s symbol is sigma2/2; its actual value is established below separately.
    scalar_direct=derivative(chi)+(2 if mutation!='omit_unit_normal' else 0)*p*u+s*chi
    gamma=rho*(h-F(5,2)*a)
    if mutation=='wrong_gamma':
        gamma=-gamma
    invmu=rho*(h-4*a+s)
    if mutation=='wrong_mu':
        invmu=rho*(h-4*a-s)
    shell_identity_source=json.loads((master/'inputs/reference/BACKGROUND_DERIVATIVE_BOUNDS.json').read_bytes())
    definitions=shell_identity_source['exact_identities']
    need(definitions['B_b']=='h_b-4*a_b+W_second(eta_b)' and
         definitions['gamma']=='rho_b*(h_b-5*a_b/2)' and
         definitions['inverse_mu']=='rho_b*B_b','actual trace coefficient source definitions')
    source_shell_env={'__builtins__':{},'rho_b':rho,'h_b':h,'a_b':a,
                      'eta_b':P.symbol('eta'),'W_second':lambda ignored:s}
    source_Bcal=eval(compile(ast.parse(definitions['B_b'],mode='eval'),'<pinned B_b>','eval'),source_shell_env)
    source_shell_env['B_b']=source_Bcal
    source_gamma=eval(compile(ast.parse(definitions['gamma'],mode='eval'),'<pinned gamma>','eval'),source_shell_env)
    source_invmu=eval(compile(ast.parse(definitions['inverse_mu'],mode='eval'),'<pinned inverse_mu>','eval'),source_shell_env)
    same('actual registered gamma definition',gamma,source_gamma)
    same('actual registered inverse-mu definition',invmu,source_invmu)
    Htrace=Wr-ell*gamma*W
    if mutation=='wrong_H_orientation':
        Htrace=Wr+ell*gamma*W
    A0=F(5,2)-F(35,4)*H**2
    if mutation=='omit_boundary_constant':
        A0=-F(35,4)*H**2
    Dtrace=Wxx-ell*invmu*Htrace+ell**2*A0*W
    if mutation=='lose_ell':
        Dtrace=Wxx-invmu*Htrace+ell**2*A0*W
    Q=2*h**2-10*a*h+F(35,4)*a**2-V2-p**2/6
    potential_path=master/'potential/check_actual_potential.py'
    potential_node=function_ast(potential_path,'run')
    q_assigns=[node for node in potential_node.body if isinstance(node,ast.Assign)
               and len(node.targets)==1 and isinstance(node.targets[0],ast.Name)
               and node.targets[0].id=='Q']
    need(len(q_assigns)==1,'unique actual potential Q AST')
    q_env={'__builtins__':{},'F':F,'h':h,'a':a,'u2':V2,'p':p}
    exec(compile(ast.Module(body=q_assigns,type_ignores=[]),str(potential_path),'exec'),q_env)
    same('actual registered master potential Q AST',Q,q_env['Q'])
    ast_receipts['actual_potential_Q']=hashlib.sha256(
        ast.dump(ast.Module(body=q_assigns,type_ignores=[]),include_attributes=False).encode()).hexdigest()
    if mutation=='wrong_potential':
        Q=Q+2*V2
    Rmaster=Wxx-Wrr+ell**2*(rho**2*Q-F(35,4)*H**2)*W
    prefactor=3*alpha*vartheta/(ell**2*rho**4*p)
    same('physical scalar residual = prefactor*(master residual−dynamic trace residual)',
         scalar_direct,prefactor*(Rmaster-Dtrace))
    same('dynamic trace retains inertial time second derivative',Dtrace.diff('Wxx'),1)
    same('dynamic/master time accelerations cancel in physical scalar normal condition',
         (Rmaster-Dtrace).diff('Wxx'),0)
    same('dynamic/master Hubble terms cancel in physical scalar normal condition',
         (Rmaster-Dtrace).diff('H'),0)
    same('mu/gamma exact shell correlation',invmu-gamma,rho*(s-F(3,2)*a))
    # Bind sigma2 to the actual registered polynomial, not a freely selected interval.
    eta=P.symbol('eta');delta=F(1,1000);model_c=F(96714321425692,161839258930731)
    Wsuper=F(1,3)+eta**2+eta**3/3
    sigma=2*Wsuper+delta*(1+model_c+model_c*eta)
    model_path=master/'potential/check_actual_potential.py'
    model_run=function_ast(model_path,'run')
    model_assigns=[node for node in model_run.body if isinstance(node,ast.Assign)
                   and len(node.targets)==1 and isinstance(node.targets[0],ast.Name)
                   and node.targets[0].id in ('W','sigma')]
    need(len(model_assigns)==2,'unique registered model tension AST closure')
    def padd(left,right):
        return [(left[j] if j<len(left) else F(0))+
                (right[j] if j<len(right) else F(0)) for j in range(max(len(left),len(right)))]
    model_env={'__builtins__':{},'F':F,'delta':delta,'c':model_c,
               'add':padd,'scale':lambda seq,factor:[v*factor for v in seq]}
    exec(compile(ast.Module(body=model_assigns,type_ignores=[]),str(model_path),'exec'),model_env)
    source_sigma=sum((coeff*eta**j for j,coeff in enumerate(model_env['sigma'])),P())
    same('actual registered tension polynomial AST',source_sigma,sigma)
    ast_receipts['actual_model_tension']=hashlib.sha256(
        ast.dump(ast.Module(body=model_assigns,type_ignores=[]),include_attributes=False).encode()).hexdigest()
    same('actual affine-detuning tension sigma2',sigma.diff('eta').diff('eta'),4+4*eta)
    same('actual half sigma2',sigma.diff('eta').diff('eta')/2,2+2*eta)
    provenance=(master/'potential/ACTUAL_POTENTIAL_WZZ_DERIVATION.md').read_text()
    need(r'\sigma=2W_{\rm super}+\delta(1+c+c\eta)' in provenance,
         'registered actual tension source')
    # Gauge trace is tangent, but tangential b,d values remain arbitrary.
    gauge_amps,gauge_scalar=env['frame_amplitudes'](W,[G['f'],G['b'],G['d']],
        rho,p,a,source_R,source_B,source_k,c,H,vartheta)
    trace_zero={n:P() for n in ('f','fx','fxx','br','brx','dr','drx')}
    if mutation=='omit_gauge_scalar_normal':
        gauge_scalar=chi
    scalar_gauge=(derivative(gauge_scalar)-p*gauge_amps['A']/2+s*gauge_scalar).substitute(trace_zero)
    same('physical scalar shell residual survives tangential-gauge pullback',scalar_gauge,scalar_direct)
    for key in ('B','C'):
        same('mixed shell metric trace '+key,gauge_amps[key].substitute(trace_zero),0)
    # Parallel physical-frame first metric R1=1/2 m_y−a A/2*eta−sigma1 chi/6*eta.
    for key,sign in (('D',-1),('L',1)):
        residual=(derivative(gauge_amps[key])/2-
                  sign*a*gauge_amps['A']/2-sign*sigma1*gauge_scalar/6).substitute(trace_zero)
        same('physical metric shell '+key+' survives tangential-gauge pullback',residual,sign*metric_from_affine)
    for key in ('E','M'):
        same('physical offdiagonal/anisotropic shell '+key,
             derivative(gauge_amps[key]).substitute(trace_zero)/2,0)
    # Explicitly retain nonzero tangential labels and independent normal lapse gauge.
    same('free tangential b value retained',gauge_amps['L'].substitute(trace_zero).diff('b'),2*H)
    same('free tangential d value retained',gauge_amps['M'].substitute(trace_zero).diff('d'),2)
    same('free radial normal gauge retained',gauge_amps['A'].substitute(trace_zero).diff('fr'),-2/(ell*rho))
    # A formal residual witness satisfies the master equation, but not the dynamic trace.
    witness={'W':P(),'Wr':P(),'Wrr':P(1),'Wxx':P(1)}
    same('master-only counterexample has zero master residual',Rmaster.substitute(witness),0)
    same('master-only counterexample has nonzero dynamic trace',Dtrace.substitute(witness),1)
    same('master-only counterexample retains physical scalar residual',scalar_direct.substitute(witness),-prefactor)
    if mutation=='omit_dynamic_premise':
        same('premature physical scalar cancellation from master alone',scalar_direct.substitute(witness),0)
    after=authenticate(repo)
    need(before==after,'complete source readback authentication')
    return {
        'status':'PASS_FIRST_MASTER_DYNAMIC_SHELL_PHYSICAL_JUNCTION_BINDING',
        'schema':'hdblast-first-shell-master-binding-v1',
        'arithmetic':'Exact collected rational Laurent coefficients; no floating decisions or project imports',
        'identities':checks,
        'source_closures':before,
        'selected_ast_sha256':ast_receipts,
        'factorization':{
            'scalar':'P1=3*alpha*vartheta/(ell^2*rho_b^4*p_b)*(Rmaster-Dtrace)',
            'metric':'R1_ab=(u_y+2*a*u+p*chi/3)*gamma0_ab=0',
            'Dtrace':'W_xx-(ell/mu)*(W_r-ell*gamma*W)+ell^2*(5/2-35*Hubble^2/4)*W',
            'Rmaster':'W_xx-W_rr+ell^2*(rho^2*Q-35*Hubble^2/4)*W',
            'Q':'2*h_ratio^2-10*a*h_ratio+35*a^2/4-V_second-p^2/6',
            'mu_inverse':'rho_b*(h_ratio-4*a_b+sigma_second/2)',
            'gamma':'rho_b*(h_ratio-5*a_b/2)',
        },
        'conventions':{
            'one_copy_outward':'K_out-sigma*gamma0/6 and n_out(phi)+sigma_prime/2',
            'inward_radial':'partial_y=-(ell*rho)^-1 partial_r',
            'gauge_shell':'f=0, b_r=d_r=0 and tangential derivatives; b,d,fr remain free',
            'trace_equation':'Retained analytic law, authenticated; not a newly executed trace solver',
        },
        'scope':{
            'existing_analytic_boundary_proof_retained':True,
            'actual_first_AST_bound':True,
            'source_reconstruction_and_trace_law_coefficient_join':True,
            'new_dynamic_trace_solution_constructed':False,
            'first_master_and_gauge_smooth_solution_existence_reproved':False,
            'second_order_solution_or_nonlinear_family_established':False,
            'cosmological_origin_established':False,
        },
    }

MUTATIONS=(
    'reverse_inward_y','wrong_scalar_reconstruction','double_tension',
    'omit_unit_normal','wrong_gamma','wrong_mu','wrong_H_orientation',
    'omit_boundary_constant','lose_ell','wrong_potential','omit_gauge_scalar_normal',
    'omit_dynamic_premise',
)
EXPECTED_FAILURE={
    'reverse_inward_y':'longitudinal master momentum reconstruction',
    'wrong_scalar_reconstruction':'longitudinal master momentum reconstruction',
    'double_tension':'one-copy physical metric shell residual',
    'omit_unit_normal':'physical scalar residual = prefactor*(master residual−dynamic trace residual)',
    'wrong_gamma':'actual registered gamma definition',
    'wrong_mu':'actual registered inverse-mu definition',
    'wrong_H_orientation':'physical scalar residual = prefactor*(master residual−dynamic trace residual)',
    'omit_boundary_constant':'physical scalar residual = prefactor*(master residual−dynamic trace residual)',
    'lose_ell':'physical scalar residual = prefactor*(master residual−dynamic trace residual)',
    'wrong_potential':'physical scalar residual = prefactor*(master residual−dynamic trace residual)',
    'omit_gauge_scalar_normal':'physical scalar shell residual survives tangential-gauge pullback',
    'omit_dynamic_premise':'premature physical scalar cancellation from master alone',
}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--repo-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--controls',action='store_true')
    parser.add_argument('--inject-error',default='none',choices=('none',)+MUTATIONS)
    args=parser.parse_args()
    need(not args.output.exists(),'output must be fresh')
    try:
        report=run(args.repo_root,args.inject_error)
        controls=[]
        if args.controls:
            with tempfile.TemporaryDirectory(prefix='first-shell-controls-') as temp:
                for mutation in MUTATIONS:
                    target=Path(temp)/(mutation+'.json')
                    cmd=[sys.executable,'-I','-S','-B']
                    if sys.flags.optimize:
                        cmd.append('-O')
                    cmd += [str(Path(__file__).resolve()),'--repo-root',str(args.repo_root),
                            '--output',str(target),'--inject-error',mutation]
                    proc=subprocess.run(cmd,capture_output=True,text=True)
                    need(proc.returncode==1 and not target.exists() and
                         proc.stderr.startswith('IDENTITY FAILED: '+EXPECTED_FAILURE[mutation]+'; residual='),
                         'corruption did not fail mathematically: '+mutation+' '+proc.stderr)
                    controls.append({'mutation':mutation,'returncode':1,'output_absent':True,
                                     'expected_failed_identity':EXPECTED_FAILURE[mutation],
                                     'failure':proc.stderr.strip()})
        report['corrupted_subprocesses']=controls
        args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    except (ValueError, OSError) as exc:
        print(str(exc),file=sys.stderr)
        return 1
    return 0

if __name__=='__main__':
    raise SystemExit(main())
