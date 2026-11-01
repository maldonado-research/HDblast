#!/usr/bin/env python3
"""Portable frozen actual-shell normal-clock intrinsic budget, stdlib only.

Requires a source repo root, authenticates closed capsule inventories before any
scientific import and after all calls, writes only a fresh external output dir.
The selected 64-history/8-target/one-point campaign is never adapted.
"""
import argparse,copy,hashlib,importlib.util,json,math,sys
from fractions import Fraction as F
from pathlib import Path
from types import SimpleNamespace
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
PINS_SHA='e07e57114f3f9cc2f134defeb7d72aa4434cfb62094d1b633fccbbddb1b532af'
FREEZE_SHA='2ba0aa65cf0942d29cd6e1bb3e3cbdfd3d4db2d07df61322aa537cc4408a4d9e'
PROTOCOL_SHA='a604c97a8f08bac1a038e09cf31990d4bfe5a01bdf91476b7c2a652dcec5cdc9'
def need(ok,msg):
    if not ok: raise ValueError(msg)

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def regular(p): return p.is_file() and not p.is_symlink() and all(not x.is_symlink() for x in p.absolute().parents)

def read(p): return json.loads(Path(p).read_bytes())

def encode(x):
    if hasattr(x,'lo') and hasattr(x,'hi'):return [str(x.lo),str(x.hi)]
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x

def authenticate(root,pins):
    need(root.is_dir() and not root.is_symlink(),'REPO_ROOT')
    result={}
    for label,pin in sorted(pins['capsules'].items()):
        cap=root/'research'/label
        need(cap.is_dir() and not cap.is_symlink(),'CAPSULE_ROOT')
        manifest=cap/'MANIFEST.json'
        need(regular(manifest) and sha(manifest)==pin,'MANIFEST:'+label)
        seen={'MANIFEST.json'}
        for row in read(manifest)['files']:
            rel=Path(row['path'])
            need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,'SAFE_MEMBER')
            q=cap/rel
            need(regular(q) and q.stat().st_size==row['bytes'] and sha(q)==row['sha256'],'MEMBER:'+label+'/'+str(rel))
            need(pins['all_manifest_member_pins'][label+'/'+str(rel)]=={'bytes':row['bytes'],'sha256':row['sha256']},'FLAT_MEMBER_PIN')
            seen.add(str(rel))
        actual=set()
        for q in cap.rglob('*'):
            need(not q.is_symlink(),'SYMLINK')
            if q.is_file():actual.add(str(q.relative_to(cap)))
            else:need(q.is_dir(),'SPECIAL_FILE')
        need(actual==seen,'CLOSED_INVENTORY')
        result[label]={'manifest_sha256':pin,'members_including_manifest':len(seen)}
    return result

def actual_master(jets,x,errors,data):
    """Actual longitudinal master u, using exact Wbar(x,0)=0 plus full errors.

    The comparison's characteristic odd extension gives zero for every purely
    temporal derivative at the physical shell. This cancellation precedes
    interval evaluation; no actual master value is replaced by its comparison.
    """
    ew=[F(v['exact_upper']) for v in errors['master_error_budget_by_scaled_total_order']['actual_minus_comparison']]
    W=jets.T({(i,0):jets.I(-ew[i],ew[i])/math.factorial(i) for i in range(4)})
    rho,eta,p,a,h,R,B,k,Vpp,decay=jets.core.background(jets.I(0),data,'none')
    c,H,vartheta=jets.temporal(x)
    U=jets.radial(R)*vartheta*W
    return {'W_scaled_time_derivatives':{str(i):W.derivative(i,0) for i in range(4)},
            'u_scaled_time_derivatives':{str(i):U.derivative(i,0) for i in range(4)},
            'u':U.val(),'u_t':U.derivative(1,0)/jets.ELL,
            'reconstruction_R':R.val(),'vartheta':vartheta.val()}

def evaluate(jets,errors,data,xb,kind):
    x=jets.I(*xb)
    m,s,b=jets.actual_jet(x,jets.I(0),mode='actual')
    master=actual_master(jets,x,errors,data)
    return {'kind':kind,'scaled_x':xb,'scaled_r':[F(0),F(0)],
            'actual_metric':m,'actual_scalar':s,'actual_background':b,'actual_master':master}

def outward(x,places=42):
    scale=10**places
    def dec(v):
        sign='-' if v<0 else '';v=abs(v)
        return sign+str(v//scale)+'.'+str(v%scale).zfill(places)
    vals=[(x.lo*scale).__floor__(),-(-x.hi*scale).__floor__()]
    out=[dec(v) for v in vals]
    need(F(out[0])<=x.lo<=x.hi<=F(out[1]),'OUTWARD_DECIMAL')
    return out

def weight(x,lo,hi,power):
    """Exact nonnegative integral on a possibly truncated source cell."""
    if x<=lo:return F(0)
    end=min(x,hi)
    return ((x-lo)**(power+1)-(x-end)**(power+1))/(power+1)

def forced_integral(I,history,xb,power):
    xl,xh=xb;value=I(0);terms=[]
    for row in history:
        lo,hi=row['scaled_x']
        if lo>=xh:break
        wl,wh=weight(xl,lo,hi,power),weight(xh,lo,hi,power)
        need(0<=wl<=wh,'NONNEGATIVE_MONOTONE_VOLTERRA_WEIGHTS')
        w=I(wl,wh);term=w*row['trace_forcing'];value+=term
        terms.append({'source_cell':[lo,hi],'kernel_power':power,'exact_weight':w,'term':term})
    return value,terms

def campaign(jets,N,ntarget,mutation=None):
    I=jets.I;data,errors,errorpins=jets.inputs();ell=jets.ELL;T=F(63,50)
    need(ell==F(1,10**6) and jets.core.GRID==1<<112,'ARITHMETIC_UNITS')
    need(errors['status']=='PASS_ACTUAL_MASTER_AND_GAUGE_CONTINUOUS_SCALED_C3_APPROXIMATION_ERRORS','ACTUAL_COMPLETE_SOURCE_ERRORS')
    shell=data['shell']['boxes']
    raw_shell={key:shell[key] for key in ['rho_b','p_b','a_b','eta_b','U_prime_b','sigma_second_b','B_b','mu']}
    p=I(*shell['p_b']['exact']);Bcal=I(*shell['B_b']['exact']);mu=I(*shell['mu']['exact'])
    rho=I(*shell['rho_b']['exact']);a=I(*shell['a_b']['exact']);Vp=I(*shell['U_prime_b']['exact']);sigma2=I(*shell['sigma_second_b']['exact'])
    Fprime=2*p*Bcal;direct=2*Vp-8*a*p+(2 if mutation=='gauss_coefficient' else 1)*sigma2*p
    need(p.hi<0 and Bcal.lo>0 and mu.lo>0 and Fprime.hi<0,'ACTUAL_TRACE_FORCE_COEFFICIENT_NONZERO')
    need(max(Fprime.lo,direct.lo)<=min(Fprime.hi,direct.hi),'TRACE_FORCE_IDENTITIES_COMPATIBLE')
    Fprime=I(max(Fprime.lo,direct.lo),min(Fprime.hi,direct.hi))
    history=[];acap=bcap=F(0);source_abs_integral=F(0)
    for i in range(N):
        xb=[T*i/N,T*(i+1)/N];row=evaluate(jets,errors,data,xb,'whole_actual_scalar_forcing_history_cell')
        b=row['actual_background'];rho=b['rho']
        if mutation=='error_gain':
            source=row['actual_scalar']['s'];mid=(source.lo+source.hi)/2;rad=(source.hi-source.lo)/4
            row['actual_scalar']['s']=I(mid-rad,mid+rad)
        forcing=ell**2*rho**2*Fprime*row['actual_scalar']['s']/6
        damp=3*ell*b['theta'];mass=ell**2*(8/b['c']**2-4)
        row.update({'trace_forcing':forcing,'trace_damping_coefficient':damp,'trace_mass_coefficient':mass})
        history.append(row);acap=max(acap,damp.mag());bcap=max(bcap,mass.mag())
        source_abs_integral+=(xb[1]-xb[0])*forcing.mag()
    need(history[0]['scaled_x'][0]==0 and history[-1]['scaled_x'][1]==T
         and all(history[i]['scaled_x'][1]==history[i+1]['scaled_x'][0] for i in range(N-1)),'WHOLE_HISTORY_COVER_FROM_ZERO')
    den=1-acap*T-bcap*T*T/2
    need(den>0,'VOLTERRA_GLOBAL_DENOMINATOR')
    Vcap=source_abs_integral/den;Ucap=T*Vcap;G=acap*Vcap+bcap*Ucap
    lo,hi=F(31,25),T
    target=[]
    specs=[('continuous_target_cell',[lo+(hi-lo)*i/ntarget,lo+(hi-lo)*(i+1)/ntarget]) for i in range(ntarget)]
    specs.append(('selected_target',[F(5,4),F(5,4)]))
    for kind,xb in specs:
        row=evaluate(jets,errors,data,xb,kind);rho=row['actual_background']['rho'];xmax=xb[1]
        vx,vterms=forced_integral(I,history,xb,0)
        uxforce,uterms=forced_integral(I,history,xb,2 if mutation=='kernel_order' else 1)
        pintforce,pterms=forced_integral(I,history,xb,2);pintforce/=(1 if mutation=='P_coefficient' else 2)
        ux=vx+I(-xmax*G,xmax*G)
        u=uxforce+I(-xmax*xmax*G/2,xmax*xmax*G/2)
        pint=pintforce+I(-xmax**3*G/6,xmax**3*G/6)
        P=rho*ell*pint
        b=row['actual_background'];A=rho*b['c'];H=b['theta']/rho;Hp=1/(A*A)
        up=ux/(ell if mutation=='clock_scale' else rho*ell);E=u-H*P;Ep=up-H*u-Hp*P;L=P/(A*A);theta=3*Ep+8*L
        C2=F(10,3)*E**2/(A*A);S2=F(80,9)*L**2;V2=theta**2/3
        Q2=F(2,9)*(theta**2-40*L**2)
        need(C2.lo>=0 and S2.lo>=0 and V2.lo>=0,'GEOMETRIC_NONNEGATIVITY')
        # This independent global-master enclosure must overlap every recovered
        # actual master interval; no intersection is used to manufacture signs.
        old=row['actual_master']
        need(max(u.lo,old['u'].lo)<=min(u.hi,old['u'].hi),'GLOBAL_MASTER_VALUE_COMPATIBILITY')
        oldux=ell*old['u_t']
        need(max(ux.lo,oldux.lo)<=min(ux.hi,oldux.hi),'GLOBAL_MASTER_DERIVATIVE_COMPATIBILITY')
        row.update({'forcing_ux_integral':vx,'forcing_u_integral':uxforce,'forcing_integrated_u':pintforce,
                    'forcing_ux_terms':vterms,'forcing_u_terms':uterms,'forcing_integrated_u_terms':pterms,
                    'ux_feedback_error_cap':xmax*G,'u_feedback_error_cap':xmax*xmax*G/2,
                    'integrated_u_feedback_error_cap':xmax**3*G/6,
                    'recovered_u':u,'recovered_u_x':ux,'history_integral_u_dx':pint,
                    'recovered_u_t':ux/ell,'u_prime':up,'P':P,'A':A,'H':H,'H_prime':Hp,
                    'E':E,'E_prime':Ep,'L':L,'theta_first_amplitude':theta,
                    'curvature_defect_second':C2,'mean_squared_shear_second':S2,
                    'expansion_variance_second':V2,'buchert_Q_second':Q2})
        target.append(row)
    need(target[0]['scaled_x'][0]==lo and target[ntarget-1]['scaled_x'][1]==hi
         and all(target[i]['scaled_x'][1]==target[i+1]['scaled_x'][0] for i in range(ntarget-1)),'CONTINUOUS_TARGET_COVER')
    keys=['recovered_u','recovered_u_x','P','E','E_prime','theta_first_amplitude','curvature_defect_second','mean_squared_shear_second','expansion_variance_second','buchert_Q_second']
    hulls={key:I(min(row[key].lo for row in target[:ntarget]),max(row[key].hi for row in target[:ntarget])) for key in keys}
    flags={'actual_master_strict_positive':hulls['recovered_u'].lo>0,'normal_clock_E_prime_strict_positive':hulls['E_prime'].lo>0,
           'curvature_defect_second_strict_positive':hulls['curvature_defect_second'].lo>0,
           'mean_squared_shear_second_strict_positive':hulls['mean_squared_shear_second'].lo>0,
           'expansion_variance_second_strict_positive':hulls['expansion_variance_second'].lo>0,
           'buchert_Q_second_strict_positive':hulls['buchert_Q_second'].lo>0}
    need((flags['curvature_defect_second_strict_positive'] and flags['expansion_variance_second_strict_positive'] and flags['buchert_Q_second_strict_positive']),'FROZEN_MAIN_SIGNS')
    return {'schema':'HDBLAST_ACTUAL_SHELL_TRACE_RECOVERY_NORMAL_CLOCK_V1',
            'status':'PASS_CONDITIONAL_ACTUAL_SHELL_TRACE_RECOVERED_NORMAL_CLOCK_GEOMETRY',
            'actual_error_input_pins':errorpins,'raw_actual_shell_coefficient_boxes':raw_shell,
            'F_prime':Fprime,'F_prime_direct':direct,'B_cal':Bcal,
            'history_rows':history,'target_rows':target,'selected_target':target[-1],
            'actual_oracle_calls':N+ntarget+1,'T':T,'damping_abs_cap':acap,'mass_abs_cap':bcap,
            'source_abs_integral':source_abs_integral,'global_denominator':den,
            'u_x_abs_cap':Vcap,'u_abs_cap':Ucap,'feedback_acceleration_error_cap':G,
            'continuous_interval_hulls':hulls,'strict_sign_flags':flags,
            'outward_decimal':{key:outward(value,42) for key,value in hulls.items()},
            'selected_target_outward_decimal':{key:outward(target[-1][key],42) for key in keys},
            'derivative_convention':'second=actual second parameter derivative; epsilon^2/2. u_x is scaled coordinate derivative; u_prime=u_x/(rho_b*ell) is proper tau derivative.',
            'scope':{'conditional_exact_family_and_actual_first_field_premises':True,
                     'linearized_intrinsic_shell_trace_identity_used':True,'actual_scalar_source_with_complete_inherited_errors':True,
                     'initial_u_and_u_x_exact_zero':True,'whole_history_from_zero_including_precontact_errors':True,
                     'finite_propagation_zero_before_one_used':False,'no_point_quadrature':True,
                     'second_metric_or_scalar_response_used':False,'completion_independent_given_same_first_shell_pullback':True,
                     'finite_amplitude_radius':False,'SI_force_or_detector_claim':False,'hot_matter_or_cosmological_origin_inference':False}}

def strip_origins(value):
    if isinstance(value,dict):return {k:strip_origins(v) for k,v in value.items() if k!='origin'}
    if isinstance(value,list):return [strip_origins(v) for v in value]
    return value

def scientific_sha(report,freeze):
    projection={key:report[key] for key in freeze['scientific_projection_keys']}
    return hashlib.sha256(json.dumps(encode(projection),sort_keys=True,separators=(',',':')).encode()).hexdigest()

def authenticate_essential():
    need(regular(HERE/'SOURCE_PINS.json') and sha(HERE/'SOURCE_PINS.json')==PINS_SHA,'SOURCE_PINS')
    need(regular(HERE/'FREEZE.json') and sha(HERE/'FREEZE.json')==FREEZE_SHA,'FREEZE')
    need(regular(HERE/'PROTOCOL_PINS.json') and sha(HERE/'PROTOCOL_PINS.json')==PROTOCOL_SHA,'PROTOCOL_PINS')
    pins,freeze,protocol=read(HERE/'SOURCE_PINS.json'),read(HERE/'FREEZE.json'),read(HERE/'PROTOCOL_PINS.json')
    for rel,pin in protocol.items():
        p=Path(rel)
        need(not p.is_absolute() and '..' not in p.parts and p.parts[0]=='protocol','SAFE_PROTOCOL_MEMBER')
        p=HERE/p
        need(regular(p) and p.stat().st_size==pin['bytes'] and sha(p)==pin['sha256'],'PROTOCOL_MEMBER:'+rel)
    actual={str(p.relative_to(HERE)) for p in (HERE/'protocol').rglob('*') if p.is_file()}
    need(actual==set(protocol),'PROTOCOL_CLOSED_INVENTORY')
    need(freeze['history_cells']==64 and freeze['target_cells']==8 and freeze['actual_oracle_calls']==73
         and freeze['history_interval']==['0','63/50'] and freeze['continuous_target_interval']==['31/25','63/50']
         and freeze['selected_target']==['5/4','0'] and freeze['source_pins_sha256']==PINS_SHA
         and freeze['no_adaptive_refinement'] is True,'FROZEN_CAMPAIGN')
    return pins,freeze

def install_raw_cache(jets,raw):
    """Review/control mode only: replay pinned raw rows, with no new oracle call."""
    rows=raw['history_rows']+raw['target_rows']
    cache={tuple(row['scaled_x']):row for row in rows}
    def interval(value):return jets.I(F(value[0]),F(value[1])) if isinstance(value,list) else F(value)
    def cached(jets_arg,errors,data,xb,kind):
        source=copy.deepcopy(cache[tuple(map(str,xb))])
        row={'kind':kind,'scaled_x':xb,'scaled_r':[F(0),F(0)]}
        for key in ['actual_metric','actual_scalar','actual_background']:
            row[key]={k:interval(v) for k,v in source[key].items()}
        master=source['actual_master']
        for key in ['W_scaled_time_derivatives','u_scaled_time_derivatives']:
            master[key]={k:interval(v) for k,v in master[key].items()}
        for key in ['u','u_t','reconstruction_R','vartheta']:master[key]=interval(master[key])
        row['actual_master']=master
        return row
    globals()['evaluate']=cached
    fake_data={'shell':{'boxes':raw['raw_actual_shell_coefficient_boxes']}}
    fake_errors={'status':'PASS_ACTUAL_MASTER_AND_GAUGE_CONTINUOUS_SCALED_C3_APPROXIMATION_ERRORS'}
    jets.inputs=lambda:(fake_data,fake_errors,raw['actual_error_input_pins'])

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',required=True,type=Path)
    ap.add_argument('--output-dir',required=True,type=Path)
    ap.add_argument('--raw-report',type=Path,help='Review/control replay only; no fresh oracle calls.')
    ap.add_argument('--inject-error',choices=['error_gain','kernel_order','clock_scale','P_coefficient','gauss_coefficient'])
    args=ap.parse_args()
    root=args.repo_root.resolve(strict=True);out=args.output_dir.absolute()
    need(not out.exists() and not out.is_symlink(),'FRESH_OUTPUT_DIRECTORY')
    need(all(not p.is_symlink() for p in out.parents),'OUTPUT_ANCESTORS')
    need(out!=HERE and HERE not in out.parents and out!=root and root not in out.parents,'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
    ownpin=sha(Path(__file__))
    pins,freeze=authenticate_essential()
    before=authenticate(root,pins)
    raw=None
    if args.raw_report:
        need(regular(args.raw_report),'RAW_REPORT_REGULAR_FILE')
        rawpin=sha(args.raw_report);raw=read(args.raw_report)
        need(scientific_sha(raw,freeze)==freeze['expected_scientific_witness_sha256'],'RAW_REPORT_SCIENTIFIC_WITNESS')
    path=root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py'
    spec=importlib.util.spec_from_file_location('authenticated_actual_intrinsic_shell_source',path)
    jets=importlib.util.module_from_spec(spec);spec.loader.exec_module(jets)
    if raw is not None:install_raw_cache(jets,raw)
    result=campaign(jets,64,8,args.inject_error)
    witness=scientific_sha(result,freeze)
    need(witness==freeze['expected_scientific_witness_sha256'],'FROZEN_SCIENTIFIC_WITNESS')
    after=authenticate(root,pins);pins_after,freeze_after=authenticate_essential()
    need(before==after and pins==pins_after and freeze==freeze_after and sha(Path(__file__))==ownpin,'ALL_SOURCES_POSTCHECK')
    if raw is not None:need(sha(args.raw_report)==rawpin,'RAW_REPORT_POSTCHECK')
    result=strip_origins(result)
    result.update({'source_authentication':before,'all_sources_authenticated_before_import_and_after':True,
                   'producer_sha256':ownpin,'source_pins_sha256':PINS_SHA,'freeze_sha256':FREEZE_SHA,
                   'protocol_pins_sha256':PROTOCOL_SHA,'frozen_campaign':freeze,
                   'scientific_witness_sha256':witness,'actual_new_oracle_calls':73 if raw is None else 0,
                   'report_kind':'FRESH_ACTUAL_ORACLE_CAMPAIGN' if raw is None else 'PINNED_RAW_REAGGREGATION_REVIEW_ONLY'})
    encoded=encode(result)
    serialized=json.dumps(encoded,sort_keys=True,indent=2)+'\n'
    need('/workspace' not in serialized and '/tmp/' not in serialized,'PATH_FREE_CANONICAL_REPORT')
    out.mkdir()
    target=out/'ACTUAL_SHELL_BUDGET.json';target.write_text(serialized)
    print(json.dumps({'status':encoded['status'],'actual_new_oracle_calls':encoded['actual_new_oracle_calls'],
                      'scientific_witness_sha256':witness,'report_sha256':sha(target),
                      'strict_sign_flags':encoded['strict_sign_flags']}))
    return 0

if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,OSError,KeyError) as exc:
        print('FAIL: '+str(exc),file=sys.stderr)
        raise SystemExit(1)
