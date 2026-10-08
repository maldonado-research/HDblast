#!/usr/bin/env python3
"""Pure JSON analysis of the frozen independent real polynomial family.

No physical callback, source coefficient, phase, trajectory or saved array is
evaluated here. An externally pinned successful entry receipt binds the JSON
budget produced by the reviewed registered algorithm. Hashes identify its
code-defined polynomials; they do not export standalone coefficient vectors.
"""
from fractions import Fraction as Q
from decimal import Decimal, localcontext, ROUND_CEILING
from pathlib import Path
import argparse
import copy
import hashlib
import json
import re

SOURCES=('positive_B','signed_uB')
SCOPE='BD_PREHISTORY_TARGET_AT_FIXED_RATIONAL_PROBES'
DEGREE=112
BITS=512
DELTA=Q(1,128)
UNIT=Q(1,2**BITS)
CONFIGURATION='INDEPENDENT_DYADIC512_SOURCE112_ODE160'

def require(ok,message):
    if not ok: raise ValueError(message)

def canonical(value):
    return str(value.numerator)+'/'+str(value.denominator)

def rational(value):
    require(type(value) is str and re.fullmatch(r'-?(?:0|[1-9][0-9]*)/[1-9][0-9]*',value),
            'Canonical exact rational string required')
    parsed=Q(value)
    require(value==canonical(parsed),'Noncanonical rational')
    return parsed

def sha(raw): return hashlib.sha256(raw).hexdigest()

def load(raw):
    def unique(items):
        out={}
        for key,value in items:
            require(key not in out,'Duplicate JSON key')
            out[key]=value
        return out
    def invalid(value): raise ValueError('Nonfinite JSON value: '+value)
    return json.loads(raw,object_pairs_hook=unique,parse_constant=invalid)

def outward(value):
    return Q((value.numerator*2**BITS+value.denominator-1)//value.denominator,2**BITS)

def decimal_upper(value):
    with localcontext() as ctx:
        ctx.prec=12
        ctx.rounding=ROUND_CEILING
        return str(Decimal(value.numerator)/Decimal(value.denominator))

def geometry():
    rows=[];left=DELTA
    while left<Q(1,2):
        right=min(3*left/2,Q(1,2));cx=(left+right)/2
        half=(right-left)/2;radius=left/2
        q=1-cx+radius;t=5-cx-radius;d=1-q*q
        require(0<half<=radius/2 and q<1 and t>0 and d>0,'Invalid analytic disk')
        mb=2*(4/t**2+4*q/(t*d**2)+2/d**2+8*q*q/d**3+4*q*q/d**4)
        mz=q*mb+2*(2/t+4*q/d**2)
        parts=(right-left)*64
        children=(parts.numerator+parts.denominator-1)//parts.denominator
        rows.append({'left':left-5,'right':right-5,'center':cx-5,
                     'half':half,'radius':radius,'M':(mb,mz),'children':children})
        left=right
    require(len(rows)==11 and sum((r['right']-r['left'] for r in rows),Q(0))==Q(1,2)-DELTA,
            'Incomplete source geometry')
    return rows

def validate_models(budget,fabricated=False):
    require(type(budget) is dict,'Budget object required')
    fixed={'schema_version':1,'scope':'FABRICATED_ONLY' if fabricated else SCOPE,
           'configuration':CONFIGURATION,'fabricated_source_provider':fabricated,
           'physical_source_evaluations':0 if fabricated else 22,
           'source_callbacks':0 if fabricated else 22,
           'source_degree':112,'coefficient_bits':512,'exp_degree':200,
           'geometric_panels':11,'ode_children_per_source':37,
           'archive_arrays_decoded':0,'retained_arrays_decoded':0,
           'original_binary80_state_error':'NOT_ENCLOSED',
           'full_twelve_case_certificate':'UNRESOLVED'}
    for key,value in fixed.items():
        require(type(budget.get(key)) is type(value) and budget[key]==value,
                'Model contract differs: '+key)
    given=budget.get('source_model_rows')
    require(type(given) is list and len(given)==22,'Exactly22 source model rows required')
    geom=geometry();seen=set();checked=[]
    for row in given:
        require(type(row) is dict and row.get('source') in SOURCES,'Unknown source model')
        parent=row.get('parent')
        require(type(parent) is int and 0<=parent<11,'Invalid parent index')
        key=(row['source'],parent)
        require(key not in seen,'Duplicate source parent')
        seen.add(key)
        expected=geom[parent];sid=SOURCES.index(row['source'])
        for field,name in [('center','center'),('half_width','half'),('analytic_radius','radius')]:
            require(rational(row.get(field))==expected[name],'Changed source geometry: '+field)
        require(row.get('interval')==[canonical(expected['left']),canonical(expected['right'])],
                'Changed parent interval or missing coverage')
        require(type(row.get('children')) is int and row['children']==expected['children'],
                'Changed child partition')
        M=expected['M'][sid]
        require(rational(row.get('disk_bound'))==M,'Changed holomorphic source majorant')
        tail=rational(row.get('analytic_tail'))
        require(tail==M/2**DEGREE,'Missing or altered source Taylor tail')
        coefficient=rational(row.get('coefficient_error'))
        uniform=rational(row.get('uniform_source_error'))
        require(coefficient>=0 and uniform>=tail,'Negative or omitted model error')
        require((coefficient/UNIT).denominator==1 and (uniform/UNIT).denominator==1,
                'Model error must retain outward dyadic512 export')
        # Both fields are separately rounded outward by less than one unit.
        # The underlying exact coefficient norm is not exported. These
        # inequalities detect incompatible rounding metadata; they do not
        # reconstruct that norm from hashes.
        require(tail+coefficient-UNIT<uniform<tail+coefficient+UNIT,
                'Uniform/coefficient error exports are inconsistent')
        digest=row.get('chosen_coefficients_sha256')
        require(type(digest) is str and re.fullmatch(r'[0-9a-f]{64}',digest),
                'Canonical complete-polynomial SHA256 required')
        safe=max(uniform,tail+coefficient)
        checked.append({'source':row['source'],'parent':parent,
            'coefficient_count_from_reviewed_degree':113,
            'chosen_coefficients_sha256':digest,
            'interval':row['interval'],'width':canonical(expected['right']-expected['left']),
            'analytic_tail':canonical(tail),'outward_coefficient_error':canonical(coefficient),
            'outward_uniform_error':canonical(uniform),'safe_source_error':canonical(safe),
            'safe_parent_L1_error':canonical((expected['right']-expected['left'])*safe)})
    require(seen=={(source,index) for source in SOURCES for index in range(11)},
            'Incomplete source-parent universe')
    return sorted(checked,key=lambda row:(row['source'],row['parent']))

def build_certificate(budget,fabricated=False):
    models=validate_models(budget,fabricated)
    interior={source:sum((rational(row['safe_parent_L1_error']) for row in models
                         if row['source']==source),Q(0)) for source in SOURCES}
    d=DELTA*(2-DELTA);ell=1/(5-DELTA);H=Q(3,8)**63
    pd=2*(1-DELTA)/d**2
    rows=[]
    for source,rho in zip(SOURCES,(0,1)):
        for K in (64,128,256):
            cap=H*(pd+rho+2*ell+2*K+DELTA*(4*K*K+4*K*ell+6*ell*ell))
            total=cap+interior[source]
            density=total*(Q(K*K,252)+Q(K,294))
            pressure=total*(Q(K**3,162)+Q(5*K,1764))
            rows.append({'source':source,'K':K,'cap_W_bound':canonical(cap),
                'interior_source_and_coefficient_L1_bound':canonical(interior[source]),
                'total_W_representation_bound':canonical(total),
                'density_representation_bound':canonical(density),
                'pressure_representation_bound':canonical(pressure),
                'density_representation_upper_decimal':decimal_upper(density),
                'pressure_representation_upper_decimal':decimal_upper(pressure)})
    return {'schema_version':1,
        'status':'FABRICATED_CONTINUOUS_MODEL_TEST_ONLY' if fabricated else 'PASS_CONTINUOUS_SOURCE_REPRESENTATION_CERTIFICATE',
        'scope':'EXACT_DUHAMEL_REAL_POLYNOMIAL_FAMILY_SOURCE_AND_COEFFICIENT_CHOICE_ONLY',
        'domain':['-5/1','-9/2'],'source_degree':112,'source_parents_each':11,
        'sources':list(SOURCES),'spectral_domain':'0<k<=K, K in {64,128,256}',
        'time_domain':['-9/2','-7/2'],
        'polynomial_family':'The reviewed registered builder defines one exact real dyadic polynomial per source and parent before any momentum loop; full ordered113-coefficient lists are SHA256-identified. Coefficient vectors are not exported here.',
        'canonical_invariant':'Exact real-source/polynomial difference has c=0; no numerical U/W box or retained mode is projected.',
        'amplitude_envelope':'|delta A(k)| <= total_W_representation_bound/(2k); no uniform infrared A bound is asserted.',
        'measure_premise':'Same inherited represented Pi; only Pi>=3 used in rational norm bounds.',
        'models':models,'bounds':rows,
        'point_phase_and_ODE_arithmetic':'EXCLUDED_FROM_THIS_CERTIFICATE; separate registered nine-probe enclosures retain them',
        'all_k_numerical_mode_error':'NOT_ENCLOSED',
        'original_binary80_state_error':'NOT_ENCLOSED',
        'full_twelve_case_certificate':'UNRESOLVED','metric_calibration':'FAIL',
        'higher_dimensional_Big_Bang_origin':'NOT_ESTABLISHED','external_novelty':'NOT_ASSESSED',
        'uv_tail':'NOT_SUPPLIED','new_physical_source_callbacks':0,'saved_array_decodes':0,
        'acceptance_gate_policy':'No new fitted gate and no alteration of the inherited2e-8 full-integral gate.'}

def synthetic_budget(extra=False):
    budget={'schema_version':1,'scope':'FABRICATED_ONLY','configuration':CONFIGURATION,
        'fabricated_source_provider':True,'physical_source_evaluations':0,'source_callbacks':0,
        'source_degree':112,'coefficient_bits':512,'exp_degree':200,'geometric_panels':11,
        'ode_children_per_source':37,'archive_arrays_decoded':0,'retained_arrays_decoded':0,
        'original_binary80_state_error':'NOT_ENCLOSED','full_twelve_case_certificate':'UNRESOLVED',
        'source_model_rows':[]}
    for source in SOURCES:
        for index,p in enumerate(geometry()):
            tail=p['M'][SOURCES.index(source)]/2**112
            exact_coefficient=Q(1,3*2**512) if extra else Q(0)
            coefficient=outward(exact_coefficient)
            uniform=outward(tail+exact_coefficient)
            fixture=[canonical(Q((-1)**j,j+1)) for j in range(113)]
            digest=sha(json.dumps(fixture,separators=(',',':')).encode())
            budget['source_model_rows'].append({'source':source,'parent':index,
                'center':canonical(p['center']),'half_width':canonical(p['half']),
                'analytic_radius':canonical(p['radius']),
                'interval':[canonical(p['left']),canonical(p['right'])],
                'children':p['children'],'disk_bound':canonical(p['M'][SOURCES.index(source)]),
                'analytic_tail':canonical(tail),'coefficient_error':canonical(coefficient),
                'uniform_source_error':canonical(uniform),'chosen_coefficients_sha256':digest})
    return budget

def self_test():
    checks=[];controls=[]
    def check(ok,name):
        require(ok,name);checks.append(name)
    def reject(name,change):
        value=synthetic_budget();change(value)
        try: build_certificate(value,True)
        except (ValueError,TypeError): controls.append(name)
        else: raise ValueError('Synthetic mutant accepted: '+name)
    base=build_certificate(synthetic_budget(),True)
    extra=build_certificate(synthetic_budget(True),True)
    check(len(base['models'])==22 and len(base['bounds'])==6,'complete22 models and6 finite-band rows')
    for a,b in zip(base['models'],extra['models']):
        aa=rational(a['safe_source_error']);bb=rational(b['safe_source_error'])
        check(0<=bb-aa<2*UNIT,'upward extra coefficient rounding covered '+a['source']+str(a['parent']))
        check(aa>=rational(a['analytic_tail']),'zero coefficient error retains analytic tail '+a['source']+str(a['parent']))
    for a,b in zip(base['bounds'],extra['bounds']):
        check(rational(a['pressure_representation_bound'])<=rational(b['pressure_representation_bound']),
              'nonnegative coefficient error cannot improve pressure '+a['source']+str(a['K']))
    # Independent monomial identity proves that omitting coefficient112 changes
    # a real source moment. The protocol degree check rejects degree111; the
    # prior numerical review tests both engines retain this coefficient.
    h=Q(1,128);moment=2*h/113
    check(moment>0 and moment!=0,'source coefficient112 has nonzero exact zero-phase moment')
    full=[canonical(Q(j==112)) for j in range(113)]
    check(sha(json.dumps(full,separators=(',',':')).encode()) !=
          sha(json.dumps(full[:-1],separators=(',',':')).encode()),
          'complete polynomial hash changes under highest coefficient omission')
    reject('missing source parent',lambda v:v['source_model_rows'].pop())
    reject('duplicate source parent',lambda v:v['source_model_rows'].__setitem__(1,copy.deepcopy(v['source_model_rows'][0])))
    reject('omitted highest source degree',lambda v:v.__setitem__('source_degree',111))
    reject('coefficient precision reduced',lambda v:v.__setitem__('coefficient_bits',511))
    reject('source tail halved',lambda v:v['source_model_rows'][0].__setitem__('analytic_tail',canonical(rational(v['source_model_rows'][0]['analytic_tail'])/2)))
    reject('holomorphic majorant changed',lambda v:v['source_model_rows'][0].__setitem__('disk_bound','1/1'))
    reject('interval coverage changed',lambda v:v['source_model_rows'][0].__setitem__('interval',['-5/1','-9/2']))
    reject('negative coefficient error',lambda v:v['source_model_rows'][0].__setitem__('coefficient_error','-1/1'))
    reject('omitted uniform error',lambda v:v['source_model_rows'][0].__setitem__('uniform_source_error','0/1'))
    reject('inconsistent upward rounding',lambda v:v['source_model_rows'][0].__setitem__('uniform_source_error',canonical(outward(rational(v['source_model_rows'][0]['analytic_tail']))+3*UNIT)))
    reject('malformed polynomial hash',lambda v:v['source_model_rows'][0].__setitem__('chosen_coefficients_sha256','0'*63))
    reject('noncanonical coefficient radius',lambda v:v['source_model_rows'][0].__setitem__('coefficient_error','0/2'))
    reject('promoted original-state status',lambda v:v.__setitem__('original_binary80_state_error','ENCLOSED'))
    reject('fabricated provider mislabelled physical',lambda v:v.__setitem__('fabricated_source_provider',False))
    return {'status':'PASS_PURE_JSON_SYNTHETIC_CONTROLS','checks':checks,'check_count':len(checks),
            'rejected_controls':controls,'rejected_control_count':len(controls),
            'synthetic_scope':'Metadata/arithmetic fixtures only; no fixture polynomial is claimed to enclose the physical bump.',
            'new_physical_source_callbacks':0,'saved_array_decodes':0,
            'postprocessor_sha256':sha(Path(__file__).read_bytes())}

def pinned_json(path,digest):
    require(type(digest) is str and re.fullmatch(r'[0-9a-f]{64}',digest),'External expected SHA256 required')
    p=Path(path)
    require(p.is_file() and not p.is_symlink() and p.stat().st_size<=20*1024*1024,'Bounded regular JSON required')
    raw=p.read_bytes();require(sha(raw)==digest,'External JSON SHA256 differs')
    return load(raw),raw

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--budget')
    parser.add_argument('--budget-sha256')
    parser.add_argument('--entry-receipt')
    parser.add_argument('--entry-receipt-sha256')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    require(not args.output.exists(),'Fresh postprocessing output required')
    if args.self_test:
        require(all(v is None for v in (args.budget,args.budget_sha256,args.entry_receipt,args.entry_receipt_sha256)),
                'Self-test accepts no physical output JSON')
        result=self_test()
    else:
        require(all((args.budget,args.budget_sha256,args.entry_receipt,args.entry_receipt_sha256)),
                'Externally pinned independent budget and entry receipt required')
        entry,entry_raw=pinned_json(args.entry_receipt,args.entry_receipt_sha256)
        budget,budget_raw=pinned_json(args.budget,args.budget_sha256)
        require(type(entry) is dict and entry.get('status')=='PASS_AUTHENTICATED_OFFLINE_REGISTERED_ENTRY'
                and entry.get('route')=='independent' and entry.get('fabricated_only') is False,
                'Successful registered independent entry required')
        require(entry.get('original_binary80_state_error')=='NOT_ENCLOSED'
                and entry.get('full_twelve_case_certificate')=='UNRESOLVED','Entry scope promotion forbidden')
        require(entry.get('files',{}).get('BUDGET.json')=={'bytes':len(budget_raw),'sha256':sha(budget_raw)},
                'Entry does not bind supplied budget bytes')
        callbacks=entry.get('source_callbacks')
        require(type(callbacks) is dict and type(callbacks.get('physical_source_callbacks')) is int
                and callbacks['physical_source_callbacks']==22 and type(callbacks.get('archive_arrays_decoded')) is int
                and callbacks['archive_arrays_decoded']==0,'Entry source inventory differs')
        result=build_certificate(budget)
        events=callbacks.get('events');require(type(events) is list and len(events)==22,'Complete source journal required')
        expected={(source,canonical(p['center'])) for source in SOURCES for p in geometry()}
        observed=[]
        for index,event in enumerate(events,1):
            require(type(event) is dict and type(event.get('index')) is int and event['index']==index
                    and event.get('route')=='independent' and event.get('scope')=='UNIQUE_ANALYTIC_BD_PREHISTORY_ONLY',
                    'Registered source journal event differs')
            rational(event.get('center'));observed.append((event.get('source'),event['center']))
        require(len(set(observed))==22 and set(observed)==expected,'Incomplete registered source journal')
        result['provenance']={'independent_budget_sha256':sha(budget_raw),
            'independent_entry_receipt_sha256':sha(entry_raw),'freeze_commit':entry.get('freeze_commit'),
            'registration_sha256':entry.get('registration_sha256'),
            'postprocessor_sha256':sha(Path(__file__).read_bytes()),
            'trust_boundary':'Externally supplied JSON pins and authenticated registered-entry receipt; this postprocessor performs no network authentication.'}
    with args.output.open('x') as handle:
        handle.write(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'new_physical_source_callbacks':0,'saved_array_decodes':0}))

if __name__=='__main__': main()
