"""Streaming exact serialized endpoint readback; no source evaluation or Arb.

This verifies exported widths, canonical correlation and all24 finite aggregates.
It relies on the authenticated worker/reviewed theorem for original-array truth
and validated flow truth, neither inferred from replay or interval self-agreement.
"""
import sys
if __name__=='__main__':
    if not sys.flags.isolated:raise ValueError('isolated Python -I required for standalone readback')
    import argparse
    import hashlib
    import json
    import os
    from pathlib import Path
    import re
    import stat
    import types
    def _standalone_capture(path,pin=None,limit=1<<20):
        with os.fdopen(os.open(path,os.O_RDONLY|os.O_NOFOLLOW|os.O_NONBLOCK),'rb') as stream:
            before=os.fstat(stream.fileno())
            if not stat.S_ISREG(before.st_mode) or before.st_nlink!=1 or not 0<before.st_size<=limit:
                raise ValueError('bounded unaliased readback bootstrap file required')
            raw=stream.read(limit+1);after=os.fstat(stream.fileno())
        if len(raw)!=before.st_size or (before.st_size,before.st_mtime_ns,before.st_ctime_ns)!=(
                after.st_size,after.st_mtime_ns,after.st_ctime_ns):
            raise ValueError('readback bootstrap file changed')
        if pin is not None and (len(raw)!=pin['bytes'] or hashlib.sha256(raw).hexdigest()!=pin['sha256']):
            raise ValueError('readback captured source byte pin mismatch')
        return raw
    def _standalone_module(name,path,raw):
        module=types.ModuleType(name);module.__file__=str(path);module.__package__=''
        sys.modules[name]=module;exec(compile(raw,str(path),'exec'),module.__dict__)
        return module
    def _standalone_pairs(items):
        out={}
        for key,value in items:
            if key in out:raise ValueError('duplicate readback registration key')
            out[key]=value
        return out
    _parser=argparse.ArgumentParser()
    _parser.add_argument('output',type=Path)
    _parser.add_argument('--registration',type=Path,required=True)
    _parser.add_argument('--registration-sha256',required=True)
    _args=_parser.parse_args();_root=Path(__file__).resolve().parents[1]
    _raw=_standalone_capture(_args.registration,limit=8<<20)
    if not re.fullmatch('[0-9a-f]{64}',_args.registration_sha256) or hashlib.sha256(_raw).hexdigest()!=_args.registration_sha256:
        raise ValueError('readback registration SHA mismatch')
    _registration=json.loads(_raw,object_pairs_hook=_standalone_pairs,
        parse_constant=lambda _:(_ for _ in ()).throw(ValueError('nonfinite readback registration')))
    _guard_path=_root/'execution/registration_guard.py';_aggregate_path=_root/'engine/endpoint_aggregate.py'
    _guard_raw=_standalone_capture(_guard_path,_registration['files']['execution/registration_guard.py'])
    _aggregate_raw=_standalone_capture(_aggregate_path,_registration['files']['engine/endpoint_aggregate.py'])
    _guard=_standalone_module('registration_guard',_guard_path,_guard_raw)
    # Read-only output validation authenticates source; it neither authorizes
    # retained inputs/source callbacks nor executes a numerical trajectory.
    _guard.authenticate(_root,_args.registration,_args.registration_sha256,True)
    _standalone_module('endpoint_aggregate',_aggregate_path,_aggregate_raw)
    _standalone_output=_args.output
from fractions import Fraction as Q
from pathlib import Path
import gzip
import hashlib
import json
import sys
sys.dont_write_bytecode=True
from registration_guard import require,load
from endpoint_aggregate import node_terms,Accumulator,SCALE

CAPSULES=tuple(s+'/'+g for s in ('positive_B','signed_uB') for g in ('coarse','fine'))
ENDPOINTS=(Q(-4),Q(-7,2))
PREFIXES={'coarse':{64:2048,128:4096,256:8192},'fine':{64:4096,128:8192,256:16384}}

def canonical(value):
    if type(value) is Q:return str(value.numerator)+'/'+str(value.denominator)
    if type(value) is dict:return {k:canonical(v) for k,v in value.items()}
    if type(value) in (list,tuple):return [canonical(v) for v in value]
    return value

def exact_pair(pair):
    require(type(pair) is list and len(pair)==2 and all(type(x) is str for x in pair),'two exact ratio strings required')
    return tuple(Q(x) for x in pair)

def validate(directory):
    directory=Path(directory);data=load((directory/'DATA.json').read_bytes())
    require(data['scope']=='SAVED_ENDPOINT_MINUS_EXACT_FLOW_FROM_EXACT_REPRESENTED_INCOMING_STATE' and
            type(data['manufactured']) is bool,'exact endpoint scope/mode required')
    manufactured=data['manufactured'];certificate=load((directory/'SOURCE_CERTIFICATE.json').read_bytes())
    require(certificate['manufactured'] is manufactured and len(certificate['rows'])==128 and
            certificate['source_count']==2 and certificate['cell_count_per_source']==64 and certificate['degree']==24,
            'complete independent source certificate schema required')
    for index,row in enumerate(certificate['rows']):
        require(row['source']==('positive_B','signed_uB')[index//64] and row['cell']==index%64 and
                Q(row['center'])==Q(-9,2)+Q(2*(index%64)+1,128),'source order/geometry changed')
        require(Q(row['half'])==Q(1,128) and Q(row['radius'])==Q(1,8) and Q(row['cauchy_majorant'])==64,
                'source geometry/majorant changed')
        co=row['coefficients'];require(type(co) is list and len(co)==25 and all(type(c) is str for c in co),'25 exact source ratios required')
        require(hashlib.sha256(json.dumps(co,separators=(',',':')).encode()).hexdigest()==row['coefficient_sha256'],
                'source coefficient hash differs')
        tail,error,ce=Q(row['analytic_tail']),Q(row['source_error']),Q(row['coefficient_error'])
        rounding=Q(row['source_error_rounding']);raw_error=tail+ce;unit=1<<512
        rounded=Q((raw_error.numerator*unit+raw_error.denominator-1)//raw_error.denominator,unit)
        require(ce>=0 and error==rounded and 0<=rounding<Q(1,unit) and error==tail+ce+rounding and
                tail==(Q(0) if manufactured else Q(1,15*2**90)),
                'source uncertainty missing/weakened/double counted')
    expected=[];maxrad={'U':Q(0),'W':Q(0)};nodes=comparisons=negative_zero_count=0
    capsule_index=0;node_index=0;acc={t:Accumulator() for t in ENDPOINTS}
    decoded_bytes=0
    with gzip.open(directory/'NODE_CERTIFICATES.jsonl.gz','rb') as stream:
        while raw:=stream.readline(65537):
            decoded_bytes+=len(raw)
            require(len(raw)<=65536 and raw.endswith(b'\n') and decoded_bytes<=512*(1<<20),
                    'bounded complete node JSON line required')
            row=load(raw);require(capsule_index<4,'too many capsule nodes')
            capsule=CAPSULES[capsule_index];grid=capsule.split('/')[1];plan=PREFIXES[grid]
            require(row['capsule']==capsule and row['index']==node_index,'capsule/node order changed')
            k,weight=Q(row['k']),Q(row['weight'])
            for K,count in plan.items():require((node_index<count and k<K) or (node_index>=count and k>K),'exact prefix mask changed')
            iu,iw=exact_pair(row['incoming_U']),exact_pair(row['incoming_W'])
            flags=row['negative_zero_flags']
            require(type(flags) is list and len(flags)==14 and all(type(x) is bool for x in flags) and
                    flags[0] is False and flags[1] is False,'14 exact decoder sign flags required')
            normalized=[*iu,*iw]
            for t in ENDPOINTS:
                item=row['endpoints'][str(t)]
                normalized.extend((*exact_pair(item['saved_U']),*exact_pair(item['saved_W'])))
            require(all(not flag or value==0 for flag,value in zip(flags[2:],normalized)),
                    'negative-zero flag attached to nonzero component')
            negative_zero_count+=sum(flags)
            require(set(row['endpoints'])=={'-4','-7/2'},'two fixed endpoint rows required')
            for t in ENDPOINTS:
                item=row['endpoints'][str(t)];target=item['target_dyadic96']
                require(set(target)=={'U','W'},'two full target components required')
                for part in ('U','W'):
                    require(set(target[part])=={'real','imag'},'two full target axes required')
                    for axis in ('real','imag'):
                        pair=target[part][axis]
                        require(type(pair) is list and len(pair)==2 and all(type(x) is int for x in pair) and pair[0]<=pair[1],
                                'ordered outward dyadic96 target endpoints required')
                    radius=Q(sum(target[part][a][1]-target[part][a][0] for a in ('real','imag')),2*SCALE)
                    require(radius<=Q(1,10**18),'complete target width gate failed')
                    maxrad[part]=max(maxrad[part],radius)
                terms=node_terms(k,weight,t,exact_pair(item['saved_U']),exact_pair(item['saved_W']),target,iu,iw)
                require(item['delta_U']==canonical(terms['delta_U']) and item['delta_W']==canonical(terms['delta_W']) and
                        item['canonical_drift']==canonical(terms['canonical_drift']),'saved minus target/canonical bridge differs')
                acc[t].add(k,weight,terms);comparisons+=1
            nodes+=1;node_index+=1
            bycount={c:K for K,c in plan.items()}
            if node_index in bycount:expected.extend(acc[t].snapshot(capsule,bycount[node_index],t) for t in ENDPOINTS)
            if node_index==plan[256]:
                capsule_index+=1;node_index=0;acc={t:Accumulator() for t in ENDPOINTS}
    require(capsule_index==4 and node_index==0 and nodes==49152 and comparisons==98304,'complete all-node endpoint universe required')
    require(data['nodes']==nodes and data['endpoint_comparisons']==comparisons and data['rows']==canonical(expected),
            'serialized finite aggregates differ from independent exact readback')
    summary=load((directory/'SCIENCE_SUMMARY.json').read_bytes())
    require(summary['maximum_complete_endpoint_export_L1_radius']==canonical(maxrad) and summary['nodes']==nodes and
            summary['endpoint_comparisons']==comparisons and summary['decoded_real_slots']==688152 and
            summary['negative_zero_flags']==negative_zero_count,'scientific summary differs')
    source_attempts=(directory/'SOURCE_ATTEMPTS.jsonl').read_bytes().splitlines()
    decode_attempts=(directory/'DECODE_ATTEMPTS.jsonl').read_bytes().splitlines()
    require(len(source_attempts)==(0 if manufactured else 128) and len(decode_attempts)==(0 if manufactured else 36),
            'source/decode actual attempt counts differ')
    return {'status':'PASS_STANDALONE_ENDPOINT_READBACK','manufactured':manufactured,'nodes':nodes,
            'endpoint_comparisons':comparisons,'prefix_cases':24,'source_coefficient_vectors':128,
            'maximum_complete_endpoint_export_L1_radius':canonical(maxrad),
            'source_evaluations':0,'retained_arrays_opened':0,
            'truth_scope':'AUTHENTICATED_WORKER_AND_REVIEWED_ANALYTIC_PROOF; READBACK_CHECKS_SERIALIZATION_AND_FINITE_AGGREGATES'}

if __name__=='__main__':print(json.dumps(validate(_standalone_output),sort_keys=True))
