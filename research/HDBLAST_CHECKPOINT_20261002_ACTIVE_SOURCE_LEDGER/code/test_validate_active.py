#!/usr/bin/env python3
"""Fabricated-record acceptance, negative science and fatal mutation checks only."""
from __future__ import annotations
import argparse
import copy
from decimal import Decimal, localcontext
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import subprocess

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('synthetic_authoritative_validator',HERE/'validate_active.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)


def text(x):
    x=Q(x)
    with localcontext() as c:
        c.prec=140
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def primary_row(setting,k,sab='-0.00001'):
    n=len(v.grid(setting));r={'K':k}
    for field in v.PRIMARY_PROFILES:r[field]=[text(t) for t in v.grid(setting)] if field=='times' else ['0']*n
    # A constant retained F gives the reported dyadic Simpson ledger. P is
    # chosen to respect F=L*(R-3P) with R=0; these remain fabricated results,
    # not solutions of a physical trajectory or a proof of any attribution.
    r['F']=r['stored_F']=[sab]*n
    r['P']=r['stored_P']=[text(v.number(sab)*t/3) for t in v.grid(setting)]
    for field in v.PRIMARY_ANCHORS:r[field]=['0']*3
    for field in v.PRIMARY_SCALARS+(v.FINE_SCALARS if setting=='fine' else ()):r[field]='0'
    pi=Q(14488038916154245685,4611686018427387904);m=Q(k*k)/(8*pi*pi)
    r['analytic_source_moment']=r['discrete_source_moment']=text(m)
    r['S_ab']=sab;r['S_direct_reset']=sab
    if setting=='fine':r['S_ab_double_global']=sab
    repair_primary(r,setting)
    return r


def repair_primary(r,setting):
    n=len(r['R']);mid=n//2;N=lambda k:v.number(r[k]);A=lambda k,j:v.number(r[k][j])
    # These algebraic fixture repairs manufacture internally coherent mocked
    # records. Tests independently assert their expected scientific semantics.
    vals={}
    vals['DeltaR']=A('stored_R',n-1)-A('stored_R',0)
    vals['predicted_delta_R']=A('R',n-1)-A('R',0)
    vals['canonical_constant_drift']=A('canonical_constant',n-1)-A('canonical_constant',0)
    vals['mixed_momentum_work']=(N('analytic_source_moment')-N('discrete_source_moment'))*N('source_work_integral')
    vals['I_A']=A('G',n-1)-A('G',0)-N('discrete_source_moment')*N('source_work_integral')+N('full_contact_ledger_integral')
    vals['contact_momentum_endpoint_difference']=(A('discrete_contact_R_anchors',2)-A('contact_R',n-1))-(A('discrete_contact_R_anchors',0)-A('contact_R',0))
    vals['E_momentum']=vals['contact_momentum_endpoint_difference']-vals['mixed_momentum_work']
    vals['I_ab']=vals['I_A']+vals['E_momentum']
    vals['D_S']=vals['DeltaR']-N('S_ab');vals['D_cont_A']=vals['DeltaR']-vals['I_A'];vals['E_Q_A']=vals['I_A']-N('S_ab')
    vals['D_cont']=vals['DeltaR']-vals['I_ab'];vals['E_Q']=vals['I_ab']-N('S_ab')
    vals['predicted_Ward_residual']=vals['predicted_delta_R']-vals['I_A']
    vals['primitive_contact_identity_gap']=vals['I_A']-(A('G',n-1)-A('G',0)+A('contact_R',n-1)-A('contact_R',0)+vals['mixed_momentum_work'])
    vals['E_operator_initial']=A('stored_R',0)-A('R',0)-(A('discrete_contact_R_anchors',0)-A('contact_R',0))
    vals['E_operator_end']=A('stored_R',n-1)-A('R',n-1)-N('E_flow')-(A('discrete_contact_R_anchors',2)-A('contact_R',n-1))
    vals['E_operator']=vals['E_operator_end']-vals['E_operator_initial']
    vals['source_primitive_residual']=vals['predicted_Ward_residual']-vals['canonical_constant_drift']+vals['mixed_momentum_work']
    vals['E_reconstruction']=vals['canonical_constant_drift']+vals['source_primitive_residual']
    vals['contact_endpoint_defect']=vals['E_operator']+vals['contact_momentum_endpoint_difference']
    vals['flow_operator_closure']=vals['D_cont_A']-N('E_flow')
    vals['matched_integral_correction']=vals['E_momentum']
    vals['S_global_minus_reset']=N('S_ab')-N('S_direct_reset')
    vals['midpoint_contact_plus_flow_density_gap']=A('stored_R',mid)-A('R',mid)
    vals['midpoint_E_momentum']=A('discrete_contact_R_anchors',1)-A('contact_R',mid)
    vals['midpoint_E_operator']=vals['midpoint_contact_plus_flow_density_gap']-N('midpoint_E_flow')-vals['midpoint_E_momentum']
    for field in ('decomposition_error','raw_analytic_decomposition_error','flow_contact_decomposition_error',
                  'complete_flow_momentum_operator_decomposition_error','matched_decomposition_error','midpoint_decomposition_error'):vals[field]=Q(0)
    vals['flow_exact_projection_gap']=N('E_flow_signed_projection')-N('E_flow')
    for f in ('R','P','F'):vals[f+'_profile_max_error']=max(abs(A(f,j)-A('stored_'+f,j)) for j in range(n))
    vals['baseline_contact_profile_max_error']=max(abs(A('baseline_contact',j)-A('stored_baseline_contact',j)) for j in range(n))
    if setting=='fine':vals['S_native_minus_double']=N('S_ab')-N('S_ab_double_global')
    for k,val in vals.items():r[k]=text(val)
    for j,ti in enumerate((0,mid,n-1)):
        r['discrete_operator_R_residual_anchors'][j]=r['E_operator_initial'] if j==0 else r['midpoint_E_operator'] if j==1 else r['E_operator_end']
        r['discrete_contact_R_minus_analytic_anchors'][j]=text(A('discrete_contact_R_anchors',j)-A('contact_R',ti))
        r['discrete_contact_P_minus_analytic_anchors'][j]=text(A('discrete_contact_P_anchors',j)-A('contact_P',ti))
        r['analytic_rho0_anchors'][j]=r['rho0'][ti];r['analytic_p0_anchors'][j]=r['p0'][ti]


def independent_row(setting,source,k,sab='-0.00001'):
    n=len(v.grid(setting));r={'source':source,'setting':setting,'K':k,'times':[text(t) for t in v.grid(setting)]}
    for field,subkeys in [('profiles',('R','P','F')),('contact_profiles',('R','P','F')),('stored_profiles',('R','P','F')),('baseline_profiles',('r0','p0'))]:
        r[field]={s:['0']*n for s in subkeys}
    r['profiles']['F']=r['stored_profiles']['F']=[sab]*n
    r['profiles']['P']=r['stored_profiles']['P']=[text(v.number(sab)*t/3) for t in v.grid(setting)]
    r['control_profiles']={c:copy.deepcopy(r['profiles']) for c in v.INDEPENDENT_CONTROLS}
    r['anchor_recomputed']={c:{f:r['profiles'][f][j] for f in ('R','P','F')} for c,j in (('initial',0),('midpoint',n//2),('endpoint',n-1))}
    pi=Q(14488038916154245685,4611686018427387904);m=text(Q(k*k)/(8*pi*pi));r['precisions']={}
    for dps in v.LEVELS:
        controls={c:{f:'0' for f in v.INDEPENDENT_CONTROL_SCALARS} for c in v.INDEPENDENT_CONTROLS}
        for c in controls:controls[c]['M_discrete']=m
        fields={f:'0' for f in v.INDEPENDENT_FIELDS+(v.FINE_SCALARS if setting=='fine' else ())}
        fields['S_ab']=fields['S_direct_reset']=sab
        if setting=='fine':fields['S_ab_double_global']=sab
        r['precisions'][dps]={'K':k,'fields':fields,'controls':controls}
    repair_independent(r)
    return r


def repair_independent(r):
    n=len(r['times']);mid=n//2;N=lambda x:v.number(x)
    sr=r['stored_profiles']['R'];pr=r['profiles']['R']
    for dps in v.LEVELS:
        level=r['precisions'][dps];f=level['fields'];c=level['controls']
        for control in c:
            c[control]['I_ab']=text(N(c[control]['I_mode'])+N(c[control]['I_contacts']))
            c[control]['Lg_projection']=text(N(c[control]['M_discrete'])*N(c[control]['J_Lg']))
        for field in v.INDEPENDENT_CONTROL_SCALARS:f[field]=c['24_24'][field]
        f['DeltaR']=text(N(sr[-1])-N(sr[0]));f['D_S']=text(N(f['DeltaR'])-N(f['S_ab']))
        f['D_cont']=text(N(f['DeltaR'])-N(f['I_ab']));f['E_Q']=text(N(f['I_ab'])-N(f['S_ab']))
        f['decomposition_error']='0'
        f['initial_R_reconstruction_error']=text(N(sr[0])-N(r['anchor_recomputed']['initial']['R']))
        f['endpoint_R_reconstruction_error']=text(N(sr[-1])-N(r['anchor_recomputed']['endpoint']['R']))
        f['E_operator']=text(N(f['endpoint_R_reconstruction_error'])-N(f['initial_R_reconstruction_error']))
        f['E_evolution_ledger']=text(N(f['D_cont'])-N(f['E_flow'])-N(f['E_operator']))
        f['midpoint_E_operator']=text(N(sr[mid])-N(r['anchor_recomputed']['midpoint']['R']))
        f['midpoint_density_gap']=text(N(sr[mid])-N(pr[mid]))
        f['midpoint_projection_closure']=text(N(f['midpoint_density_gap'])-N(f['midpoint_E_flow'])-N(f['midpoint_E_operator']))
        for name,a,b in [('forcing_control_gap','24_24','16_24'),('ledger_control_gap','16_24','16_16'),('joint_control_gap','24_24','16_16')]:f[name]=text(N(c[a]['I_ab'])-N(c[b]['I_ab']))
        for key in ('R','P','F'):f[key+'_profile_max_error']=text(max(abs(N(a)-N(b)) for a,b in zip(r['profiles'][key],r['stored_profiles'][key])))
        if r['setting']=='fine':f['S_native_minus_double']=text(N(f['S_ab'])-N(f['S_ab_double_global']))


def fixtures(sab='-0.00001'):
    p={k:None for k in v.PRIMARY_TOP}
    p.update(schema_version=1,route='primary_phase_green_primitive_integer_reductions',status='COMPLETED_SAVED_DATA_DIAGNOSTIC',old_metric_status='FAIL',fixed_cases=12,interval=['-4.5','-3.5'],registration_sha256='a'*64,freeze_commit='b'*40,input_manifest_sha256='c'*64,producer_sha256='d'*64,source_unchanged_during_execution=True,records=[],failures=[],fatal_failures=[],attribution_cases=[],classification='NO_GATE_SCALE_ATTRIBUTION',precision_gap_max_exact_rational='0',quadrature_control_gap_max_exact_rational='0',resources={'seconds':1,'peak_rss_kib':1000,'scope':'fabricated','limit_seconds':900,'limit_peak_rss_kib':262144},arithmetic={k:'fabricated' for k in ('phase_source','quadrature_contact_integrands','MP_scope','quad_error','moment_quantization_bound_scope','gap_universe')})
    p['arithmetic'].update(decimal_precisions=[80,100],fractional_bits=[266,333])
    for source in v.SOURCES:
        for setting in v.SETTINGS:
            case={'source':source,'setting':setting,'controls':[]}
            for order in v.PRIMARY_CONTROLS:case['controls'].append({'quadrature_order':order,'levels':{d:[primary_row(setting,k,sab) for k in v.KS] for d in v.LEVELS},'serialization_max_error_by_dps':{d:'0' for d in v.LEVELS}})
            p['records'].append(case)
    i={k:None for k in v.INDEPENDENT_TOP}
    i.update(schema_version=1,old_metric_status='FAIL',fixed_cases=12,input_manifest_sha256='c'*64,resource_limits={'wall_seconds':900,'peak_rss_kib':262144},route='independent_direct_F_nested_Duhamel',status='COMPLETED_DIAGNOSTIC_WITHOUT_RECLASSIFYING_OLD_FAIL',provenance={'freeze_commit':'b'*40,'registration_sha256':'a'*64,'registered_files_checked':2,'input_manifest_sha256':'c'*64},interval=['-4.5','-3.5'],midpoint='-4',controls=[[16,16],[16,24],[24,24]],canonical_control=[24,24],reduction_decimal_digits=[80,100],source_phase_arithmetic='fabricated',quadrature_constants='fabricated',high_precision_scope='fabricated',units={k:'fabricated' for k in ('profiles_R_P','F','baseline_r0_p0','contact_F','baseline_work_integral')},momentum_target='fabricated',rows=[independent_row(t,s,k,sab) for t in v.SETTINGS for s in v.SOURCES for k in v.KS],serialization_max_error_by_dps={d:'0' for d in v.LEVELS},native_profile_serialization_max_error='0',native_profile_serialization_scope='fabricated',limits=['fabricated'],resources={'elapsed_seconds':1,'peak_rss_kib':1000,'scope':'fabricated wrapper','limit_seconds':900,'limit_peak_rss_kib':262144},producer_sha256='d'*64)
    i['provenance'].update(full_frozen_verification={'public_freeze_commit':'b'*40,'registration_sha256':'a'*64,'freeze_receipt_sha256':'e'*64,'frozen_files':{'inputs/INPUT_MANIFEST.json':'c'*64,'independent/diagnostic_independent.py':'d'*64},'package_manifest_sha256':None,'input_capsule_verification':{'capsules_verified':4,'members_verified':76,'array_payload_values_interpreted':False,'input_manifest_sha256':'c'*64}},post_run_frozen_verification='PASS')
    i['units']['baseline_r0_p0']='a0^4 times physical renormalized baseline rho0,p0'
    receipts=[{'name':name,'physical_route':True,'status':'PASS_EXECUTION','exit_code':0,'timed_out':False,'memory_exceeded':False,'timeout_seconds':900,'memory_limit_kib':262144,'process_creation_limit':0,'executor_uid':1000,'elapsed_seconds':2,'peak_rss_kib':1100,'wait4_peak_rss_kib':1100,'observed_peak_group_rss_kib':1000} for name in ('primary','independent')]
    seal(p)
    return p,i,receipts


def seal(p):
    # Recompute producer-declaration metadata without accepting a classification
    # from it: cross_validate checks this declaration again independently.
    p['precision_gap_max_exact_rational']=str(max(v.maximum_gap(c['levels']['80'],c['levels']['100']) for case in p['records'] for c in case['controls']))
    p['quadrature_control_gap_max_exact_rational']=str(max(v.maximum_gap(case['controls'][0]['levels'][d],case['controls'][1]['levels'][d]) for case in p['records'] for d in v.LEVELS))
    _,failures,canonical,_=v.check_primary(p)
    p['classification'],p['attribution_cases']=v.classification(failures,canonical)
    p['failures']=[]
    for f in failures:
        out={k:f[k] for k in ('source','setting','dps','gate')}
        if 'K' in f:out['K']=f['K'];out['order']=int(f['control'])
        if f['gate']=='quadrature_control_gap':out['exact_gap']=f['exact_rational_value']
        elif f['gate'] in v.PRIMARY_SCIENCE_ANCHORS:
            case=next(c for c in p['records'] if (c['source'],c['setting'])==(f['source'],f['setting']))
            control=next(c for c in case['controls'] if c['quadrature_order']==int(f['control']))
            out['values']=next(r for r in control['levels'][f['dps']] if r['K']==f['K'])[f['gate']]
        else:out['value']=text(Q(f['exact_rational_value']))
        p['failures'].append(out)


def each_first_primary(p):
    return [c['levels'][d][0] for c in p['records'][0]['controls'] for d in v.LEVELS]


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():raise RuntimeError('Fresh guard receipt required')
    checks=[]
    def accepted(name,payload,classification):
        result=v.cross_validate(*payload)
        if result['classification']!=classification:raise RuntimeError(name+' wrong classification '+result['classification'])
        checks.append({'name':name,'expected':classification,'actual':result['classification']})
    def rejected(name,mutate):
        payload=fixtures();mutate(*payload)
        try:v.cross_validate(*payload)
        except (v.FatalValidation,KeyError) as e:checks.append({'name':name,'expected':'FATAL_REJECTION','actual':type(e).__name__})
        else:raise RuntimeError('Fatal mutation survived: '+name)
    accepted('strict same-case attribution witness',fixtures(),'LEDGER_ERROR_DEMONSTRATED')
    accepted('no gate-scale witness',fixtures('0'),'NO_GATE_SCALE_ATTRIBUTION')
    payload=fixtures();p,i,_=payload
    for row in each_first_primary(p):row['source_jet_profile_max_error']='0.0000003'
    seal(p);accepted('source reconstruction negative is retained execution',payload,'CONSISTENCY_FAILURE')
    payload=fixtures();p,i,_=payload
    for row in each_first_primary(p):row['forcing_jet_profile_max_error']='0.0000003'
    seal(p);accepted('forcing reconstruction negative is retained execution',payload,'CONSISTENCY_FAILURE')
    payload=fixtures();p,i,_=payload
    for row in each_first_primary(p):row['full_contact_ledger_integral']='-0.00001';repair_primary(row,'coarse')
    seal(p);accepted('contact reconstruction/continuous negative is retained',payload,'CONSISTENCY_FAILURE')
    payload=fixtures();p,i,_=payload
    for row in each_first_primary(p):row['midpoint_E_flow']='0.00001';row['midpoint_flow_triangle_bound']='0.00001';row['stored_R'][64]='0.00001';repair_primary(row,'coarse')
    i['rows'][0]['stored_profiles']['R'][64]='0.00001';i['rows'][0]['anchor_recomputed']['midpoint']['R']='0.00001'
    for d in v.LEVELS:
        for c in i['rows'][0]['precisions'][d]['controls'].values():c['midpoint_E_flow']='0.00001';c['midpoint_triangle_bound']='0.00001'
    repair_independent(i['rows'][0]);seal(p);accepted('midpoint flow negative is retained',payload,'CONSISTENCY_FAILURE')
    payload=fixtures();p,i,_=payload
    for d in v.LEVELS:i['rows'][0]['precisions'][d]['controls']['16_16']['I_mode']='0.00001'
    repair_independent(i['rows'][0]);accepted('independent ledger-order control negative is retained',payload,'CONSISTENCY_FAILURE')
    rejected('missing primary control',lambda p,i,r:p['records'][0]['controls'].pop())
    rejected('duplicate primary control',lambda p,i,r:p['records'][0]['controls'].append(copy.deepcopy(p['records'][0]['controls'][0])))
    rejected('missing independent control',lambda p,i,r:i['rows'][0]['precisions']['80']['controls'].pop('16_24'))
    rejected('missing precision context',lambda p,i,r:p['records'][0]['controls'][0]['levels'].pop('80'))
    rejected('extra precision context',lambda p,i,r:p['records'][0]['controls'][0]['levels'].update({'120':[]}))
    rejected('unknown primary scientific scalar',lambda p,i,r:p['records'][0]['controls'][0]['levels']['80'][0].update(unregistered='0'))
    rejected('unknown independent scientific scalar',lambda p,i,r:i['rows'][0]['precisions']['80']['fields'].update(unregistered='0'))
    rejected('hidden extra top-level scalar',lambda p,i,r:p.update(unregistered_science='0'))
    rejected('decimal value passed through binary64',lambda p,i,r:p['records'][0]['controls'][0]['levels']['80'][0].update(E_flow=0.0))
    rejected('nonfinite decimal',lambda p,i,r:p['records'][0]['controls'][0]['levels']['80'][0].update(E_flow='NaN'))
    rejected('infinite decimal',lambda p,i,r:i['rows'][0]['precisions']['80']['fields'].update(E_flow='Infinity'))
    rejected('definitional signed decomposition arithmetic mutation',lambda p,i,r:p['records'][0]['controls'][0]['levels']['80'][0].update(D_S='0.00002'))
    rejected('matched integral correction sign mutation',lambda p,i,r:p['records'][0]['controls'][0]['levels']['80'][0].update(I_ab='0.00001'))
    rejected('omitted fourth reconstruction term',lambda p,i,r:p['records'][0]['controls'][0]['levels']['80'][0].pop('E_reconstruction'))
    rejected('primitive arithmetic definition mutation',lambda p,i,r:p['records'][0]['controls'][0]['levels']['80'][0].update(I_A='0.00001'))
    rejected('analytic source normalization mutation',lambda p,i,r:p['records'][0]['controls'][0]['levels']['80'][0].update(analytic_source_moment='1'))
    rejected('flow triangle mutation',lambda p,i,r:i['rows'][0]['precisions']['80']['controls']['16_16'].update(triangle_bound='-1'))
    rejected('understated arithmetic precision universe',lambda p,i,r:p['records'][0]['controls'][0]['levels']['80'][0].update(moment_quantization_only_bound='0.0000000000005'))
    rejected('actual arithmetic precision gap failure',lambda p,i,r:p['records'][0]['controls'][0]['levels']['80'][0].update(moment_quantization_only_bound='0.000000000002'))
    rejected('serialization receipt mutation',lambda p,i,r:p['records'][0]['controls'][0]['serialization_max_error_by_dps'].update({'80':'1/100000000000'}))
    rejected('native serialization receipt mutation',lambda p,i,r:i.update(native_profile_serialization_max_error='1/100000000000'))
    rejected('forged author classification',lambda p,i,r:p.update(classification='NO_GATE_SCALE_ATTRIBUTION'))
    rejected('forged attribution membership',lambda p,i,r:p.update(attribution_cases=[]))
    rejected('source changed during execution',lambda p,i,r:p.update(source_unchanged_during_execution=False))
    rejected('outer timeout failure',lambda p,i,r:r[0].update(timed_out=True))
    rejected('outer memory failure',lambda p,i,r:r[0].update(peak_rss_kib=262145))
    rejected('inner forged resource report',lambda p,i,r:p['resources'].update(seconds=3))
    rejected('root executor defeats NPROC containment',lambda p,i,r:r[0].update(executor_uid=0))
    rejected('missing UID containment evidence',lambda p,i,r:r[0].pop('executor_uid'))
    rejected('NPROC containment absent',lambda p,i,r:r[0].update(process_creation_limit=1))
    rejected('surviving route descendant',lambda p,i,r:r[0].update(surviving_descendant_pids=[123]))
    rejected('RSS witness forged',lambda p,i,r:r[0].update(wait4_peak_rss_kib=1200))
    rejected('fixed time grid mutation',lambda p,i,r:p['records'][0]['controls'][0]['levels']['80'][0]['times'].__setitem__(1,'-4.49'))
    rejected('baseline units mutation',lambda p,i,r:i['units'].update(baseline_r0_p0='physical rho0'))
    rejected('boolean schema version',lambda p,i,r:p.update(schema_version=True))
    rejected('floating independent control order',lambda p,i,r:i['controls'][0].__setitem__(0,16.0))
    rejected('floating canonical control order',lambda p,i,r:i.update(canonical_control=[24.0,24.0]))
    rejected('floating declared precision',lambda p,i,r:p['arithmetic']['decimal_precisions'].__setitem__(0,80.0))
    rejected('boolean NPROC declaration',lambda p,i,r:r[0].update(process_creation_limit=False))
    rejected('independent producer provenance mutation',lambda p,i,r:i.update(producer_sha256='0'*64))
    rejected('independent capsule payload interpretation mutation',lambda p,i,r:i['provenance']['full_frozen_verification']['input_capsule_verification'].update(array_payload_values_interpreted=True))
    rejected('independent post-run verification missing',lambda p,i,r:i['provenance'].update(post_run_frozen_verification='MISSING'))
    rejected('independent registered inventory count forged',lambda p,i,r:i['provenance'].update(registered_files_checked=3))
    rejected('floating independent resource limits',lambda p,i,r:i['resource_limits'].update(wall_seconds=900.0,peak_rss_kib=262144.0))
    with tempfile.TemporaryDirectory() as td:
        path=Path(td)/'duplicate.json';path.write_text('{"a":1,"a":2}')
        try:v.read_json(path)
        except v.FatalValidation:checks.append({'name':'duplicate JSON key','expected':'FATAL_REJECTION','actual':'FatalValidation'})
        else:raise RuntimeError('Duplicate JSON key survived')
    # A pinned malicious helper must never execute when an earlier bootstrap
    # property fails. Every checkpoint below is fabricated and has no inputs.
    for mutation in ('root_symlink','code_symlink','self_byte_mismatch','helper_byte_mismatch'):
        with tempfile.TemporaryDirectory() as td:
            td=Path(td);root=td/'checkpoint';root.mkdir();code=root/'code';code.mkdir()
            entry=code/'validate_active.py';entry.write_bytes((HERE/'validate_active.py').read_bytes())
            marker=td/'IMPORTED_HELPER';helper=code/'active_integrity.py'
            helper.write_text('from pathlib import Path\nPath('+repr(str(marker))+').write_text("bad import")\n')
            registration={'files':{'code/validate_active.py':hashlib.sha256(entry.read_bytes()).hexdigest(),
                                   'code/active_integrity.py':hashlib.sha256(helper.read_bytes()).hexdigest()}}
            if mutation=='self_byte_mismatch':registration['files']['code/validate_active.py']='0'*64
            if mutation=='helper_byte_mismatch':registration['files']['code/active_integrity.py']='0'*64
            rp=root/'FULL_REGISTRATION.json';rp.write_text(json.dumps(registration,sort_keys=True))
            if mutation=='root_symlink':
                alias=td/'aliased';alias.symlink_to(root,target_is_directory=True);root=alias;entry=root/'code/validate_active.py'
            if mutation=='code_symlink':
                external=td/'external_code';code.rename(external);code.symlink_to(external,target_is_directory=True)
            command=[sys.executable]+(['-O'] if sys.flags.optimize else [])+[str(entry),
                 '--checkpoint-root',str(root),'--registration-sha256',hashlib.sha256(rp.read_bytes()).hexdigest(),
                 '--freeze-commit','b'*40,'--primary',str(td/'absent-primary'),
                 '--independent',str(td/'absent-independent'),'--execution-receipts',str(td/'absent-receipts'),
                 '--output-dir',str(td/'absent-output')]
            run=subprocess.run(command,capture_output=True,text=True,timeout=20)
            if run.returncode==0 or marker.exists():raise RuntimeError('Pre-import bootstrap mutation survived: '+mutation)
            if 'FatalValidation' not in run.stderr:raise RuntimeError('Bootstrap failed for unrelated reason: '+run.stderr)
            checks.append({'name':'pre-import '+mutation,'expected':'FATAL_BEFORE_HELPER_IMPORT','actual':'FATAL_BEFORE_HELPER_IMPORT'})
    receipt={'status':'PASS_FABRICATED_VALIDATOR_GUARDS','check_count':len(checks),'checks':checks,
             'physical_evaluations':0,'saved_scientific_values_loaded':0,'python_optimization':sys.flags.optimize,
             'validator_sha256':hashlib.sha256((HERE/'validate_active.py').read_bytes()).hexdigest(),
             'guard_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'scope':'Manufactured exact-decimal record fixtures only. No source model, mode propagator, numerical array reader or actual route data are executed.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:receipt[k] for k in ('status','check_count','physical_evaluations')}))
if __name__=='__main__':main()
