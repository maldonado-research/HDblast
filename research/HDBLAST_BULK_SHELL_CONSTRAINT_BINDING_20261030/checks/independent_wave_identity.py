#!/usr/bin/env python3
"""Independent stdlib exact Laurent-polynomial check; imports no producer result."""
from __future__ import annotations
import argparse, ast, hashlib, json, subprocess, sys, tempfile
from fractions import Fraction as F
from pathlib import Path

SOURCE_BASE_COMMIT='d98c641aabd658ff1bf456f29ae8b69886d52e54'
MANIFEST_PIN='adb6e3a7e00a8e477af80a49bc6b7a0f4037fafd0dd42ab8175c8b55e3e4465c'
EXPECTED_PINS={
    'actual':'a68fa8025d2170efca6c3e85d443910a796d0c30a9b7dd17ed988cd9e00bbb00',
    'dc':'b53ccfd7c638c9b98df23803b50a6821ae46416adfb95f8a68a9e35b48b93f4f',
    'initial_checker':'6b5ec8e23a6dea34273b9ddac1fdd85d3a0461ddf66bea810f3042096dd83ee5',
    'original_operator':'163553f2bdb9469b3ae40bd8f1c9bcf01077b690abac9a2a7b2513619ac99011',
    'pullback':'d15f8c65ffffd0085d195064ccea3c61bb4e719695652d0dd73472c8a06c224d',
    'raw_system_note':'1fec12bc9f729ef7ca93fd4ec06fdce8fc1f61a2d295b6b04e3925ba609b479c',
    'system':'00f654647d3cd5578e3f9afcb00dc3ec81739036f40adbe43f53e84b87506e98',
}
CONTROLS={
    'wrong_scalar_shift':'scalar plus pullback',
    'missing_radial_rho':'pullback A',
    'raw_f_normalized':'pullback A',
    'wrong_clock_term':'pullback L',
    'wrong_radial_drift':'primary operator row f',
    'wrong_source_sign':'primary operator source',
    'wrong_time_sign':'time normalization',
    'missing_angular_c':'angular normalization',
    'omit_variable_rho':'ordinary physical derivative yy',
    'zero_initial_acceleration':'initial velocity B',
}


class P:
    def __init__(self, value=0):
        if isinstance(value, P): self.terms = dict(value.terms)
        elif isinstance(value, dict): self.terms = {m:F(v) for m,v in value.items() if v}
        else: self.terms = {():F(value)} if value else {}
    @staticmethod
    def var(s): return P({((s,1),):1})
    def __add__(self, other):
        other=P(other); out=dict(self.terms)
        for m,v in other.terms.items():
            out[m]=out.get(m,0)+v
            if not out[m]: del out[m]
        return P(out)
    __radd__=__add__
    def __neg__(self): return P({m:-v for m,v in self.terms.items()})
    def __sub__(self,other): return self+-P(other)
    def __rsub__(self,other): return P(other)+-self
    def __mul__(self,other):
        out={}
        for ma,va in self.terms.items():
            for mb,vb in P(other).terms.items():
                m=dict(ma)
                for s,n in mb: m[s]=m.get(s,0)+n
                m=tuple(sorted((s,n) for s,n in m.items() if n))
                out[m]=out.get(m,0)+va*vb
        return P(out)
    __rmul__=__mul__
    def __pow__(self,n):
        if n<0:
            if len(self.terms)!=1: raise ValueError('inverse must be monomial')
            ((m,v),)=self.terms.items()
            return P({tuple((s,n*k) for s,k in m):v**n})
        result=P(1)
        for _ in range(n): result=result*self
        return result
    def __truediv__(self,other): return self*P(other)**-1
    def __rtruediv__(self,other): return P(other)*self**-1
    def derivative(self,axis):
        out=P()
        for mon,coef in self.terms.items():
            for s,n in mon:
                ds=dvar(s,axis)
                if ds.terms:
                    rest=dict(mon); rest[s]-=1
                    rest=tuple(sorted((k,v) for k,v in rest.items() if v))
                    out+=P({rest:coef*n})*ds
        return out
    def dx(self): return self.derivative('x')
    def dr(self): return self.derivative('r')
    def dy(self): return -self.dr()/(ell*rho)
    def dt(self): return self.dx()/ell
    def serial(self):
        return [{'coefficient':str(v),'monomial':dict(m)} for m,v in sorted(self.terms.items())]


ell,rho,a,ay,H,q=(P.var(s) for s in ('ell','rho','a','ay','H','q'))
c=1/q


def dvar(s,axis):
    if s=='rho': return P() if axis=='x' else -ell*a*rho**2
    if s=='a': return P() if axis=='x' else -ell*rho*ay
    if s=='H': return ell*q**2 if axis=='x' else P()
    if s=='q': return -ell*H*q if axis=='x' else P()
    if s in ('ell','ay','k','p'): return P()
    name,_,suffix=s.partition(':')
    if name not in ('u','f','b','d'): return P()
    return P.var(name+':'+''.join(sorted(suffix+axis)))


def need_zero(p,label):
    if p.terms: raise RuntimeError('IDENTITY_FAILED:'+label+':'+json.dumps(p.serial(),sort_keys=True,separators=(',',':')))


def reduce_cosh(p):
    # q^2=sech^2(t)=1-H^2. Canonical q exponent 0 or 1.
    # Clear Laurent q denominators by a nonzero monomial before reducing.
    minimum=min((dict(mon).get('q',0) for mon in p.terms),default=0)
    if minimum<0: p=p*q**(-minimum)
    out=P()
    for mon,coef in p.terms.items():
        m=dict(mon); n=m.pop('q',0)
        if n<0: raise RuntimeError('unexpected negative q in identity remainder')
        out+=P({tuple(sorted(m.items())):coef})*q**(n%2)*(1-H**2)**(n//2)
    return out


def checked(p,label): need_zero(reduce_cosh(p),label)


def dc_rows(amplitudes):
    # Independent retained connection contraction; M multiplies untraced HessY.
    A,B,C,D,E,L,M=(amplitudes[k] for k in ('A','B','C','D','E','L','M'))
    J=L-F(8,3)*M; beta=F(-10,3)
    normal=(A.dy()+D.dy()-3*J.dy())/2-B.dt()/rho+a*(4*A+D-3*J)-3*H*B/rho-8*q*C/rho
    time=B.dy()-(A.dt()+D.dt()+3*J.dt())/(2*rho)+5*a*B-3*H*(D+J)/rho-8*q*E/rho
    angular=C.dy()-E.dt()/rho+5*a*C-4*H*E/rho+q*(-A/2+D/2-J/2+beta*M)/rho
    return normal,time,angular


def raw_residuals(f,b,d):
    # Registered raw P,Q,R, transcribed from the retained regularity premise.
    return (
        f.dx().dx()-f.dr().dr()+3*ell*H*f.dx()+2*ell*rho**2*a*b.dx()+3*ell*rho*a*f.dr()
        +ell**2*((8*q**2+4*rho**2*ay+8*rho**2*a**2)*f+6*rho**2*a*H*b-16*rho**2*a*d),
        b.dx().dx()-b.dr().dr()+2*ell*a*f.dx()+3*ell*H*b.dx()+5*ell*rho*a*b.dr()
        -ell**2*((3+3*H**2-8*q**2)*b-16*H*d),
        d.dx().dx()-d.dr().dr()+5*ell*H*d.dx()+5*ell*rho*a*d.dr()
        -ell**2*(2*a*q**2*f+2*H*q**2*b-4*q**2*d),
    )


def pin(path): return hashlib.sha256(path.read_bytes()).hexdigest()


class JetAdapter:
    def __init__(self,p): self.p=p
    def val(self): return self.p
    def derivative(self,i,j):
        p=self.p
        for _ in range(i): p=p.dx()
        for _ in range(j): p=p.dr()
        return p


def initial(p,acceleration):
    out=P()
    for mon,v in p.terms.items():
        term=P(v)
        for s,n in mon:
            name,_,suffix=s.partition(':')
            if name in ('f','b','d'):
                replacement=acceleration if name=='f' and suffix=='xx' else P()
                term*=replacement**n
            else: term*=P.var(s)**n
        out+=term
    return out


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--repo-root','--repo',dest='repo_root',required=True,type=Path); ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--run-controls',action='store_true');ap.add_argument('--mutation',choices=CONTROLS)
    args=ap.parse_args(); root=args.repo_root.resolve()
    if args.output.exists(): raise RuntimeError('REFUSE_OVERWRITE')
    files={
        'actual':root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py',
        'pullback':root/'research/HDBLAST_COMMON_GAUGE_FORCING_20261016/first-fields/FIRST_HARMONIC_FIELDS_UNPUBLISHED.md',
        'dc':root/'research/HDBLAST_HARMONIC_CONSTRAINT_VALUES_20261030/theory/CONSTRAINT_FORMULA_NOTE.md',
        'system':root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/inputs/GAUGE_WAVE_REGULARITY_PREMISE.md',
        'original_operator':root/'research/HDBLAST_FIRST_GAUGE_WAVE_20261016/operator/check_gauge_operator.py',
        'raw_system_note':root/'research/HDBLAST_FIRST_GAUGE_WAVE_20261016/coefficients/inputs/FIRST_GAUGE_VECTOR_SYSTEM_UNPUBLISHED.md',
        'initial_checker':root/'research/HDBLAST_REGISTERED_INITIAL_CONSTRAINTS_20261030/checks/check_registered_adm_binding.py',
    }
    before={k:pin(v) for k,v in files.items()}
    for key,path in files.items():
        if path.is_symlink() or before[key]!=EXPECTED_PINS[key]: raise RuntimeError('source authentication '+key)
    gauge=root/'research/HDBLAST_FIRST_GAUGE_WAVE_20261016'
    gauge_manifest=gauge/'MANIFEST.json'
    copied_manifest=root/'research/HDBLAST_COMMON_GAUGE_FORCING_20261016/bulk-producer/inputs/PUBLISHED_FIRST_GAUGE_MANIFEST.json'
    if pin(gauge_manifest)!=MANIFEST_PIN: raise RuntimeError('first-wave manifest authentication')
    if gauge_manifest.read_bytes()!=copied_manifest.read_bytes(): raise RuntimeError('copied first-wave manifest mismatch')
    manifest=json.loads(gauge_manifest.read_bytes())
    members={'MANIFEST.json'}
    for entry in manifest['files']:
        relative=Path(entry['path'])
        if relative.is_absolute() or '..' in relative.parts or relative.as_posix() in members: raise RuntimeError('unsafe capsule path')
        members.add(relative.as_posix())
        path=gauge/relative
        if path.is_symlink() or not path.is_file(): raise RuntimeError('capsule member type')
        value=path.read_bytes()
        if len(value)!=entry['bytes'] or hashlib.sha256(value).hexdigest()!=entry['sha256']:
            raise RuntimeError('first gauge capsule manifest mismatch '+entry['path'])
    if any(path.is_symlink() or (path.is_file() and path.relative_to(gauge).as_posix() not in members) for path in gauge.rglob('*')):
        raise RuntimeError('first gauge capsule closed tree')
    tree=ast.parse(files['actual'].read_text())
    function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='frame_amplitudes')
    conversion=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='physical_derivatives')
    isolated=ast.Module(body=[function,conversion],type_ignores=[])
    scope={'ELL':ell}; exec(compile(isolated,str(files['actual']),'exec'),scope)
    u,f,b,d=(P.var(s) for s in ('u','f','b','d')); p,k=P.var('p'),P.var('k')
    actual_derivatives=scope['physical_derivatives'](JetAdapter(f),JetAdapter(rho),JetAdapter(a))
    if args.mutation=='omit_variable_rho': actual_derivatives['yy']=f.dr().dr()/(ell*rho)**2
    derived_derivatives={'':f,'y':f.dy(),'t':f.dt(),'yy':f.dy().dy(),'yt':f.dy().dt(),'tt':f.dt().dt()}
    for suffix,value in actual_derivatives.items():
        need_zero(value-derived_derivatives[suffix],'ordinary physical derivative '+suffix)
    injected={'wrong_scalar_shift':'wrong_scalar_shift','missing_radial_rho':'drop_radial_rho'}.get(args.mutation,'none')
    raw_f=rho*f if args.mutation=='raw_f_normalized' else f
    actual,scalar=scope['frame_amplitudes'](u,(raw_f,b,d),rho,p,a,P(1),P(1),k,c,H,P(1),injected)
    if args.mutation=='wrong_clock_term': actual['L']=2*(u+a*f+b.dx()/ell)
    Rcoef,Bcoef,vartheta=(P.var(s) for s in ('Rcoef','Bcoef','vartheta'))
    full_amps,full_scalar=scope['frame_amplitudes'](u,(f,b,d),rho,p,a,Rcoef,Bcoef,k,c,H,vartheta)
    need_zero(full_scalar-(Bcoef*vartheta*(u.dr()/ell+k*u)+p*f),'unabridged scalar reconstruction')
    # Derived independently from coordinate Lie derivative and frame conversion.
    expected={'A':-4*u+2*f.dy(),'B':f.dt()/rho-rho*b.dy(),
              'C':f/(rho*c)+rho*c*d.dy(),'D':-2*(u+a*f+b.dt()),
              'E':c*d.dt()-b/c,'L':2*(u+a*f+H*b),'M':2*d}
    for name in expected: need_zero(actual[name]-expected[name],'pullback '+name)
    longitudinal_weight={'A':-4,'D':-2,'L':2}
    for name in expected:
        need_zero(full_amps[name]-expected[name]-longitudinal_weight.get(name,0)*(Rcoef*vartheta-1)*u,'unabridged reconstruction '+name)
    need_zero(scalar-(u.dr()/ell+k*u+p*f),'scalar plus pullback')
    normal,time,angular=dc_rows(actual)
    resf,resb,resd=raw_residuals(f,b,d)
    if args.mutation=='wrong_radial_drift': resf+=ell*rho*a*f.dr()
    source=-6*ell**2*rho**2*(u.dy()+4*a*u)
    if args.mutation=='wrong_source_sign': source=-source
    # Independently compare each of the 27 primary implementation entries.
    # Evaluate only list RHS ASTs, never the original checker or its reports.
    original_tree=ast.parse(files['original_operator'].read_text())
    build=next(n for n in original_tree.body if isinstance(n,ast.FunctionDef) and n.name=='build')
    op_env={'H':H,'a':a,'rho':rho,'ay':ay,'lam':P(8),'S':q**2,'u':u,'uz':u.dr()/ell}
    original={}
    for n in build.body:
        if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name) and n.targets[0].id in ('A','B','D','force'):
            original[n.targets[0].id]=eval(compile(ast.Expression(n.value),str(files['original_operator']),'eval'),{'__builtins__':{}},op_env)
    raw=(f,b,d); direct=[]
    for i in range(3):
        rhs=sum((ell*original['A'][i][j]*raw[j].dx()+ell*original['B'][i][j]*raw[j].dr()+ell**2*original['D'][i][j]*raw[j] for j in range(3)),P())
        direct.append(raw[i].dx().dx()-raw[i].dr().dr()-rhs)
    for name,ours,theirs in zip(('f','b','d'),(resf,resb,resd),direct):
        checked(ours-theirs,'primary operator row '+name)
    checked(ell**2*original['force'][0]-source,'primary operator source')
    checked(ell**2*rho**2*normal-(source-resf),'normal normalization')
    checked(ell**2*rho*time+(resb if args.mutation=='wrong_time_sign' else -resb),'time normalization')
    checked(ell**2*rho*(1 if args.mutation=='missing_angular_c' else q)*angular+resd,'angular normalization')
    # Specialize this same reconstruction/connection, using PDE-fixed f_xx.
    acceleration=P() if args.mutation=='zero_initial_acceleration' else source
    init={name:initial(value,acceleration) for name,value in actual.items()}
    init_x={name:initial(value.dx(),acceleration)/ell for name,value in actual.items()}
    wanted={'A':-4*u,'B':P(),'C':P(),'D':-2*u,'E':P(),'L':2*u,'M':P()}
    wanted_t={name:value.dt() for name,value in wanted.items()}
    wanted_t['B']=-6*rho*(u.dy()+4*a*u)
    for name in wanted:
        need_zero(init[name]-wanted[name],'initial '+name)
        need_zero(init_x[name]-wanted_t[name],'initial velocity '+name)
    # Read the exact completed initial-QC input table from its source AST.
    old_tree=ast.parse(files['initial_checker'].read_text())
    old_check=next(n for n in old_tree.body if isinstance(n,ast.FunctionDef) and n.name=='check')
    old_env={'L':P,'rho':rho,'Uy':u.dy(),'a':a,'U':u,'T':u.dt()}
    for n in old_check.body:
        is_m_or_bt=isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name) and n.targets[0].id in ('m','bt')
        is_update=isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Attribute) and isinstance(n.value.func.value,ast.Name) and n.value.func.value.id=='m' and n.value.func.attr=='update'
        if is_m_or_bt or is_update:
            exec(compile(ast.Module(body=[n],type_ignores=[]),str(files['initial_checker']),'exec'),{'__builtins__':{}},old_env)
    for name,value in actual.items():
        for suffix,generated in (('',initial(value,source)),('y',initial(value.dy(),source)),('t',initial(value.dt(),source))):
            need_zero(generated-old_env['m'][name+suffix],'completed initial QC table '+name+suffix)
    for name,row in zip(('normal','time','angular'),(normal,time,angular)):
        checked(initial(row,source),'initial DC '+name)
    # One fixture per dangerous normalization; report nonzero formal remainder.
    wrongs={
        'wrong_timelike_sign':ell**2*rho*time+resb,
        'missing_angular_c':ell**2*rho*angular+resd,
        'wrong_f_radial_drift':2*ell*rho*a*f.dr(),
        'zero_initial_acceleration':initial(actual['B'].dx(),P())/ell-wanted_t['B'],
        'wrong_scalar_shift':-2*p*f,
        'variable_rho_second_derivative_omitted':a*f.dr()/(ell*rho),
        'raw_f_replaced_by_rho_f':rho*f-f,
    }
    for label,value in wrongs.items():
        if not value.terms: raise RuntimeError('undetected control '+label)
    def verify_readback():
        after={k:pin(v) for k,v in files.items()}
        if before!=after or any(path.is_symlink() or not path.is_file() for path in files.values()):
            raise RuntimeError('source changed during independent review')
        if gauge_manifest.is_symlink() or pin(gauge_manifest)!=MANIFEST_PIN or gauge_manifest.read_bytes()!=copied_manifest.read_bytes():
            raise RuntimeError('gauge manifest readback mismatch')
        for entry in manifest['files']:
            path=gauge/entry['path'];value=path.read_bytes()
            if path.is_symlink() or len(value)!=entry['bytes'] or hashlib.sha256(value).hexdigest()!=entry['sha256']:
                raise RuntimeError('gauge capsule readback mismatch '+entry['path'])
        if any(path.is_symlink() or (path.is_file() and path.relative_to(gauge).as_posix() not in members) for path in gauge.rglob('*')):
            raise RuntimeError('gauge closed tree readback mismatch')
    verify_readback()
    report={'status':'PASS','arithmetic':'stdlib Fraction exact Laurent polynomials',
            'source_sha256':before,'source_readback_unchanged':True,
            'source_base_commit':SOURCE_BASE_COMMIT,'source_runtime_requires_git':False,'full_first_gauge_capsule_members_authenticated':len(manifest['files']),
            'full_first_gauge_capsule_manifest_sha256':pin(gauge_manifest),
            'identities':['ell^2 rho^2 DC_y = s_f - Res_f','ell^2 rho DC_0 = Res_b','ell^2 rho/c DC_a = -Res_d D_aY'],
            'checks':{'ordinary_physical_derivative_rows':6,'pullback_metric_rows':7,'scalar_pullback':1,'unabridged_reconstruction_rows':8,'primary_operator_rows':3,'primary_operator_source':1,'operator_normalizations':3,'initial_value_rows':7,'initial_velocity_rows':7,'initial_DC_rows':3,'completed_initial_QC_input_table_rows':21},
            'rejection_remainders':{k:v.serial() for k,v in wrongs.items()},
            'premises':['actual solution satisfies registered PDE as neighborhood identity','zero initial raw Cauchy data','retained analytic existence/regularity and physical domain','retained full frame DC formula and S3 identities'],
            'limitations':['no produced result imported','no analytic PDE solve','no shell/Codazzi binding','no nonlinear or cosmological conclusion']}
    if args.run_controls:
        failures=[]
        with tempfile.TemporaryDirectory(prefix='independent-first-wave-') as temporary:
            for mutation,label in CONTROLS.items():
                destination=Path(temporary)/(mutation+'.json')
                command=[sys.executable,'-I','-S','-B']+(['-O'] if sys.flags.optimize else [])+[str(Path(__file__).resolve()),'--repo-root',str(root),'--output',str(destination),'--mutation',mutation]
                result=subprocess.run(command,capture_output=True,text=True)
                marker='RuntimeError: IDENTITY_FAILED:'+label+':'
                if result.returncode==0 or destination.exists() or marker not in result.stderr:
                    raise RuntimeError('CONTROL_NOT_REJECTED:'+mutation+':'+result.stderr[-500:])
                remainder=json.loads(result.stderr.split(marker)[-1].splitlines()[0])
                if not remainder: raise RuntimeError('EMPTY_CONTROL_REMAINDER:'+mutation)
                failures.append({'mutation':mutation,'expected_identity':label,'returncode':result.returncode,'report_absent':True,'nonzero_exact_remainder':remainder})
        report['rejected_corrupted_executions']=failures
    report['checker_sha256']=pin(Path(__file__))
    verify_readback()
    out=args.output
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as fjson: json.dump(report,fjson,indent=2,sort_keys=True); fjson.write('\n')
    print('PASS independent exact first-wave normalization and initial specialization')


if __name__=='__main__': main()
