"""Standalone exact readback of every exported node and source vector.

No retained inputs are opened and no physical source/target is constructed.
The exported gzip is sufficient to recompute every finite transport summary.
"""
from fractions import Fraction as Q
from pathlib import Path
from math import factorial
import argparse
import gzip
import importlib.util
import json
import re
import sys
from registration_guard import (CAPSULES,SOURCES,authenticate,authenticate_local,
                               require,load,rational,safe_file,real_directory,sha,reauthenticate)

def _transport(root):
    path=safe_file(root,'transport/discrete_transport.py')
    spec=importlib.util.spec_from_file_location('readback_discrete_transport',path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    return module

def _read_json(directory,name):
    path=safe_file(directory,name)
    require(path.stat().st_size<=20*1024*1024,'Oversized output JSON')
    return load(path.read_bytes())

def _source_budgets(certificate):
    out=[]
    for source in SOURCES:
        rows=[r for r in certificate['rows'] if r['source']==source]
        source_l1=Q(0);source_error=Q(0)
        for row in rows:
            half=rational(row['half_width']);co=tuple(rational(q) for q in row['coefficients'])
            source_l1+=2*half*sum((abs(c)*half**j for j,c in enumerate(co)),Q(0))
            source_error+=2*half*rational(row['uniform_error'])
        n=832;x=Q(256)
        tails={'W':source_l1*x**(n+1)/factorial(n+1)/(1-x/Q(n+2)),
               'U':source_l1*x**(n+1)/(2*factorial(n+2))/(1-x/Q(n+3))}
        canon=lambda q:f'{q.numerator}/{q.denominator}'
        out.append({'source':source,'source_L1_majorant':canon(source_l1),'source_error_L1':canon(source_error),
                    'kernel_tails':{k:canon(q) for k,q in tails.items()}})
    return out

def validate_target(target,k):
    require(type(target) is dict and set(target)=={'U','W'},'Complete target rectangle required')
    radii={}
    for part in ('U','W'):
        require(type(target[part]) is dict and set(target[part])=={'real','imag'},'Target axes differ')
        for axis in ('real','imag'):
            ends=target[part][axis]
            require(type(ends) is list and len(ends)==2 and all(type(x) is int for x in ends) and ends[0]<=ends[1],
                    'Exact ordered dyadic96 target endpoints required')
        if k==0:require(target[part]['imag']==[0,0],'k0 target imaginary coordinates must be exact zero')
        radii[part]=Q(sum(target[part][axis][1]-target[part][axis][0] for axis in ('real','imag')),2*(1<<96))
        require(radii[part]<=Q(1,10**18),'Target export gate exceeded on readback')
    return radii

def _spec_hash(spec):
    return sha(json.dumps(spec,sort_keys=True,separators=(',',':'),ensure_ascii=True).encode('ascii'))

def validate_input_authentication(records,root,fabricated):
    require(type(records) is list and len(records)==4,'Four ordered input authentication records required')
    document=load(safe_file(root,'provenance/INPUT_SPEC.json').read_bytes())
    selection=['k.npy','momentum_weights.npy','observation_eta.npy','u_1.npy','w_1.npy']
    for ident,row,capsule in zip(CAPSULES,records,document['capsules']):
        require(type(row) is dict and row['capsule']==ident and ident==capsule['source']+'/'+capsule['grid'],
                'Input authentication capsule order differs')
        if fabricated:
            require(set(row)=={'capsule','fabricated','file_sha256','input_spec_sha256','authenticated_members',
                               'selected_member_names','reader_inputspec'} and row['fabricated'] is True and
                    row['authenticated_members']==19 and row['selected_member_names']==selection,
                    'Manufactured input authentication scope differs')
            spec=row['reader_inputspec'];template=capsule['reader_inputspec']
            require(type(spec) is dict and set(spec)==set(template) and spec['kind']=='npz' and
                    spec['storage_layout']==template['storage_layout'] and type(spec['file_bytes']) is int and
                    0<spec['file_bytes']<8*1024*1024 and len(spec['members'])==19,
                    'Manufactured reader spec shape/schema differs')
            for m,pinned in zip(spec['members'],template['members']):
                require(type(m) is dict and set(m)==set(pinned)-{'payload_bytes'},'Manufactured member schema differs')
                for key in ('name','descr','shape','npy_version','fortran_order','decode'):
                    require(m[key]==pinned[key], 'Manufactured member layout/selection differs')
            require(spec['file_sha256']==row['file_sha256'] and _spec_hash(spec)==row['input_spec_sha256'],
                    'Manufactured reader spec/hash linkage differs')
            for field in ('file_sha256','input_spec_sha256'):
                require(type(row[field]) is str and re.fullmatch('[0-9a-f]{64}',row[field]),
                        'Manufactured authentication digest required')
        else:
            spec={k:v for k,v in capsule['reader_inputspec'].items() if k!='members'}
            spec['members']=[{k:v for k,v in m.items() if k!='payload_bytes'} for m in capsule['reader_inputspec']['members']]
            expected={'capsule':ident,'repository_path':capsule['repository_path'],
                      'file_sha256':spec['file_sha256'],'input_spec_sha256':_spec_hash(spec),
                      'full_member_count':19,'authenticated_members':19,'selected_member_names':selection,
                      'original_archive_stream_sha256':capsule['original']['sha256']}
            require(row==expected,'Actual input authentication differs from frozen InputSpec')
    return True

def validate_output(output,candidate,fabricated=False,check_receipt=True,verified=None):
    output=real_directory(output);root=real_directory(candidate);transport=_transport(root)
    from source_certificate import validate_certificate
    data=_read_json(output,'DATA.json');summary=_read_json(output,'SCIENCE_SUMMARY.json')
    require(data==summary and safe_file(output,'DATA.json').read_bytes()==safe_file(output,'SCIENCE_SUMMARY.json').read_bytes(),
            'DATA and scientific summary differ')
    require(data['scope']==('FABRICATED_ONLY' if fabricated else 'ALL_NODE_EXACT_STORED_MINUS_CERTIFIED_BD_INCOMING_TARGET') and
            data['row_count']==12 and data['distinct_capsule_nodes']==49152 and data['prefix_node_appearances']==86016,
            'Scientific output scope/universe differs')
    require(data['represented_epsilon']==f'{transport.EPSILON.numerator}/{transport.EPSILON.denominator}' and
            data['represented_Pi']==f'{transport.PI.numerator}/{transport.PI.denominator}' and
            data['export_endpoint_exponent']==-96 and data['state_projection'] is False and
            data['diagnostic_interpolation_certifies_all_nodes'] is False,'Normalization/export/projection differs')
    for key,val in {'full_twelve_case_continuous_pressure_contacts':'UNRESOLVED','metric_calibration':'FAIL',
                    'higher_dimensional_Big_Bang_origin':'NOT_ESTABLISHED','external_novelty':'NOT_ASSESSED'}.items():
        require(data[key]==val,'Scientific limitation changed: '+key)
    expected={'physical_source_callbacks':0 if fabricated else 22,'retained_arrays_decoded':0 if fabricated else 20,
              'fabricated_arrays_decoded':20 if fabricated else 0,'selected_scalar_binary80_slots':294936,
              'actual_nodes':0 if fabricated else 49152,'fabricated_nodes':49152 if fabricated else 0,
              'target_node_evaluations':49152,'diagnostic_target_evaluations':18,'total_target_evaluations':49170,
              'physical_trajectories':0,'likelihood_evaluations':0,'unselected_member_decodes':0,
              'authenticated_capsules':4,'authenticated_members':76}
    require(data['counter']==expected and all(type(x) is int for x in data['counter'].values()),'Payload counters differ')
    validate_input_authentication(data['input_authentication'],root,fabricated)
    certificate=_read_json(output,'SOURCE_CERTIFICATE.json')
    require(all(r['fabricated'] is fabricated for r in certificate['rows']),'Source certificate mode differs')
    proof=validate_certificate(certificate,root)
    require(_source_budgets(certificate)==data['source_budgets'],'Exported source model budgets differ')
    diagnostic=_read_json(output,'DIAGNOSTIC_CROSSCHECK.json')
    require(diagnostic['target_evaluations']==18 and diagnostic['methods']==2 and
            diagnostic['rectangle_pairs_per_method']==36 and diagnostic['real_component_intersections_per_method']==72 and
            diagnostic['fabricated'] is fabricated and diagnostic['certifies_all_nodes_by_interpolation'] is False and
            len(diagnostic['rows'])==18,'Fixed diagnostic crosscheck differs')
    momenta=('0/1','1/1099511627776','1/4096','1/4','1/1','16/1','64/1','128/1','256/1')
    for index,row in enumerate(diagnostic['rows']):
        require(type(row) is dict and set(row)=={'source','momentum','target','method_intersections'} and
                row['source']==SOURCES[index//9] and row['momentum']==momenta[index%9],
                'Diagnostic row schema/identity changed')
        validate_target(row['target'],rational(row['momentum']))
        require(type(row['method_intersections']) is dict,'Diagnostic method result mapping required')
        if fabricated:require(row['method_intersections']=={},'Fabricated diagnostics cannot claim prior physical intersections')
    if not fabricated:
        require(diagnostic['total_rectangle_intersections']==72 and diagnostic['total_real_component_intersections']==144,
                'Both prior method diagnostic counts differ')
        prior={m:_read_json(root/'prior',m+'_DATA.json') for m in ('primary','independent')}
        for index,row in enumerate(diagnostic['rows']):
            require(row['source']==SOURCES[index//9] and row['momentum']==momenta[index%9], 'Diagnostic row identity changed')
            require(set(row['method_intersections'])==set(prior),'Both diagnostic methods required')
            for method,old in prior.items():
                require(type(row['method_intersections'][method]) is dict and
                        set(row['method_intersections'][method])=={'U/real','U/imag','W/real','W/imag'},
                        'Complete diagnostic component intersection schema required')
                matches=[r for r in old['whole_rows'] if r['source']==row['source'] and r['momentum']==row['momentum']]
                require(len(matches)==1,'Prior diagnostic identity differs')
                for part in ('U','W'):
                    for axis in ('real','imag'):
                        lo,hi=(Q(i,1<<96) for i in row['target'][part][axis])
                        oldlo,oldhi=(rational(matches[0][part][axis][bound]) for bound in ('lo','hi'))
                        require(max(lo,oldlo)<=min(hi,oldhi) and row['method_intersections'][method][part+'/'+axis] is True,
                                'Diagnostic intersection changed')
    else:
        require(diagnostic['total_rectangle_intersections']==0 and diagnostic['total_real_component_intersections']==0,
                'Fabricated diagnostics cannot claim physical crosschecks')
    nodefile=safe_file(output,'NODE_TARGETS.jsonl.gz')
    with nodefile.open('rb') as raw:
        header=raw.read(10)
    require(len(header)==10 and header[:4]==b'\x1f\x8b\x08\x00' and header[4:8]==b'\x00'*4,
            'Deterministic gzip header differs')
    nodes=0;expanded=0;flags=0;maxrad={'U':Q(0),'W':Q(0)};rows=[];coordinate={}
    with gzip.open(nodefile,'rb') as gz:
        for capsule in CAPSULES:
            grid=capsule.split('/')[1];plan=transport.PREFIXES[grid];n=plan[256]
            acc=transport.PrefixAccumulator();coordinate[capsule]=[];bycount={count:cutoff for cutoff,count in plan.items()}
            for index in range(n):
                raw=gz.readline(65537);expanded+=len(raw)
                require(raw.endswith(b'\n') and len(raw)<=65536 and expanded<=128*1024*1024,
                        'Missing/oversized expanded node row')
                row=load(raw)
                require(type(row) is dict and set(row)=={'capsule','index','k','weight','stored_raw_u','stored_raw_w',
                                                        'negative_zero_flags','target_dyadic96'} and
                        row['capsule']==capsule and type(row['index']) is int and row['index']==index,
                        'Exported node identity/schema differs')
                k,w=rational(row['k']),rational(row['weight'])
                require(0<k<256 and w>0,'Exported exact k/weight invalid')
                for cutoff,count in plan.items():
                    require((index<count and k<cutoff) or (index>=count and k>cutoff),'Exported prefix mask differs')
                u,ww=tuple(rational(q) for q in row['stored_raw_u']),tuple(rational(q) for q in row['stored_raw_w'])
                require(len(u)==len(ww)==2,'Exported raw complex state shape differs')
                negative=row['negative_zero_flags']
                require(type(negative) is list and len(negative)==6 and all(type(x) is bool for x in negative),
                        'Signed-zero metadata differs')
                require(all(not flag or q==0 for flag,q in zip(negative,(k,w,*u,*ww))),
                        'Negative-zero flag attached to nonzero represented component')
                flags+=sum(negative)
                target=row['target_dyadic96']
                radii=validate_target(target,k)
                for part,radius in radii.items():maxrad[part]=max(maxrad[part],radius)
                error=transport.from_saved_target(k,w,tuple(q/transport.EPSILON for q in u),
                                                  tuple(q/transport.EPSILON for q in ww),target)
                acc.add(error);coordinate[capsule].append((k,w));nodes+=1
                if index+1 in bycount:rows.append(acc.snapshot(capsule,Q(bycount[index+1])))
            require(acc.count==n,'Exported full capsule incomplete')
        require(gz.read(1)==b'','Unexpected extra exported nodes')
    require(nodes==49152 and transport.encode_exact(rows)==data['rows'],'All-node exact transport readback differs')
    canon=lambda q:f'{q.numerator}/{q.denominator}'
    require({k:canon(q) for k,q in maxrad.items()}==data['max_complete_target_export_L1_radius'] and
            flags==data['negative_zero_component_flags'],'Observed export radius/zero flags differ')
    equality={}
    for grid in ('coarse','fine'):
        a,b=coordinate['positive_B/'+grid],coordinate['signed_uB/'+grid]
        equality[grid]={'k_values_equal':all(x[0]==y[0] for x,y in zip(a,b)),
                        'weight_values_equal':all(x[1]==y[1] for x,y in zip(a,b)),
                        'compared_node_pairs':len(a),'basis':'EXACT_DECODED_VALUE_PAIRS_AFTER_SELECTED_DECODE'}
    require(equality==data['source_coordinate_value_equality'],'Value-level cross-source coordinate comparison differs')
    source_attempts=safe_file(output,'SOURCE_ATTEMPTS.jsonl').read_bytes().splitlines()
    decode_attempts=safe_file(output,'DECODE_ATTEMPTS.jsonl').read_bytes().splitlines()
    require(len(source_attempts)==(0 if fabricated else 22) and len(decode_attempts)==(0 if fabricated else 20),
            'Actual attempt journals differ')
    if not fabricated:
        for i,raw in enumerate(source_attempts):
            event=load(raw);row=certificate['rows'][i]
            require(event=={'index':i+1,'event':'real_source_taylor','source':row['source'],
                           'center':row['center'],'scope':'UNIQUE_ANALYTIC_BD_PREHISTORY_ONLY'},
                    'Source attempt universe/order differs')
        sequence=[(c,m) for c in CAPSULES for m in ('observation_eta.npy','k.npy','momentum_weights.npy','u_1.npy','w_1.npy')]
        for i,(raw,(c,m)) in enumerate(zip(decode_attempts,sequence)):
            require(load(raw)=={'index':i+1,'capsule':c,'member':m,'event':'selected_retained_array_decode_attempt'},
                    'Decode attempt universe/order differs')
    if check_receipt:
        receipt=_read_json(output,'ENTRY_RECEIPT.json')
        require(type(receipt) is dict and set(receipt)=={'status','fabricated_only','registration_sha256',
                'public_go_sha256','freeze_commit','uid','runtime','denied_syscalls','counter','files'},
                'Exact entry receipt schema required')
        status='PASS_FABRICATED_STORED_COMPARISON' if fabricated else 'PASS_AUTHENTICATED_STORED_COMPARISON'
        require(receipt['status']==status and receipt['fabricated_only'] is fabricated and receipt['counter']==expected,
                'Entry receipt status/scope/counters differ')
        require(receipt['registration_sha256']==sha(safe_file(root,'FULL_REGISTRATION.json').read_bytes()),
                'Entry receipt registration linkage differs')
        require(type(receipt['uid']) is int and receipt['uid']>0,'Entry receipt nonroot identity differs')
        require(type(receipt['denied_syscalls']) is list and
                {'socket','socketpair','connect','clone','execve','io_uring_setup','io_uring_enter','io_uring_register'}<=
                set(receipt['denied_syscalls']),'Entry receipt kernel denial policy incomplete')
        require(type(receipt['runtime']) is dict and receipt['runtime'].get('packages')==
                {'python-flint':'0.9.0','sympy':'1.14.0','mpmath':'1.3.0'} and
                receipt['runtime'].get('platform')=='linux','Entry receipt runtime requirements differ')
        if fabricated:
            require(receipt['public_go_sha256'] is None and receipt['freeze_commit'] is None,
                    'Fabricated entry receipt cannot claim public authorization')
        else:
            require(type(verified) is dict and verified.get('fabricated_only') is False and
                    verified['root']==root and verified['registration_sha256']==receipt['registration_sha256'] and
                    verified['public_go_sha256']==receipt['public_go_sha256'] and
                    verified['receipt']['freeze_commit']==receipt['freeze_commit'],
                    'Actual entry receipt requires current externally authenticated PUBLIC_GO linkage')
        names={'SOURCE_CERTIFICATE.json','DIAGNOSTIC_CROSSCHECK.json','NODE_TARGETS.jsonl.gz',
               'DATA.json','SCIENCE_SUMMARY.json','RESOURCE.json','SOURCE_ATTEMPTS.jsonl',
               'DECODE_ATTEMPTS.jsonl','OUTPUT_READBACK.json'}
        require(type(receipt['files']) is dict and set(receipt['files'])==names,
                'Complete exact entry payload inventory required')
        for name,item in receipt['files'].items():
            require(type(item) is dict and set(item)=={'bytes','sha256'} and type(item['bytes']) is int and
                    0<=item['bytes']<=128*1024*1024 and type(item['sha256']) is str and
                    re.fullmatch('[0-9a-f]{64}',item['sha256']),'Invalid entry payload pin')
            path=safe_file(output,name)
            require(path.stat().st_size==item['bytes'] and sha(path.read_bytes())==item['sha256'],
                    'Entry output hash differs: '+name)
    return {'status':'PASS_STANDALONE_EXACT_ALL_NODE_READBACK','fabricated_only':fabricated,
            'nodes_verified':nodes,'prefix_cases_verified':12,'prefix_node_appearances':86016,
            'source_certificate':proof,'retained_arrays_opened':0,'physical_source_constructed':False,
            'independent_readback_scope':'SERIALIZED_TARGET_WIDTHS_SOURCE_BUDGETS_AND_FINITE_TRANSPORT_AGGREGATES',
            'original_array_to_export_and_target_truth':'RELY_ON_AUTHENTICATED_PINNED_WORKER_AND_REVIEWED_ANALYTIC_PROOF',
            'max_complete_target_export_L1_radius':{k:canon(q) for k,q in maxrad.items()}}

def compare_replays(first,second):
    first,second=real_directory(first),real_directory(second)
    stable=('DATA.json','SCIENCE_SUMMARY.json','SOURCE_CERTIFICATE.json','DIAGNOSTIC_CROSSCHECK.json',
            'NODE_TARGETS.jsonl.gz','SOURCE_ATTEMPTS.jsonl','DECODE_ATTEMPTS.jsonl','OUTPUT_READBACK.json')
    for name in stable:
        require(safe_file(first,name).read_bytes()==safe_file(second,name).read_bytes(),'Replay exact bytes differ: '+name)
    return {'status':'PASS_EXACT_NORMAL_OPTIMIZED_REPLAY','scientific_files_exact':list(stable),
            'excluded_resource_fields_only':['RESOURCE.json wall_seconds','RESOURCE.json peak_rss_kib'],
            'DATA_exact':True,'SCIENCE_SUMMARY_exact':True}

def main():
    p=argparse.ArgumentParser();p.add_argument('--candidate',required=True);p.add_argument('--output',required=True)
    p.add_argument('--registration-sha256',required=True);p.add_argument('--public-go');p.add_argument('--public-go-sha256')
    p.add_argument('--fabricated-only',action='store_true');p.add_argument('--compare-output');a=p.parse_args()
    if a.fabricated_only:
        require(a.public_go is None and a.public_go_sha256 is None,'Fabricated readback cannot claim PUBLIC_GO')
        v=authenticate_local(a.candidate,a.registration_sha256)
    else:
        require(a.public_go is not None and a.public_go_sha256 is not None,'Actual readback requires external PUBLIC_GO')
        v=authenticate(a.candidate,a.registration_sha256,a.public_go,a.public_go_sha256)
    result=validate_output(a.output,v['root'],a.fabricated_only,verified=v)
    if a.compare_output:result['replay']=compare_replays(a.output,a.compare_output)
    reauthenticate(v)
    print(json.dumps(result,sort_keys=True,separators=(',',':')))

if __name__=='__main__':
    main()
