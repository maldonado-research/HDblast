"""Independent exact arithmetic review of actual exported endpoint certificates.

No production module, Arb, source builder, retained archive or decoder is used.
Inverse epsilon normalization and binary80 representability are checked from
exports only. Original stored value identity and validated target-flow truth
remain grounded in the authenticated producer, decoder and analytic proofs;
this review performs no second retained decode or physical reconstruction.
"""
from fractions import Fraction as Q
from pathlib import Path
import gzip
import hashlib
import json
import re
import sys
import time

BASE=Path(__file__).resolve().parent
PI=Q(14488038916154245685,4611686018427387904)
EPS=Q(3777893186295716171,37778931862957161709568)
SCALE=1<<96
def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def load(raw):
    def pairs(items):
        d={}
        for k,v in items:
            require(k not in d,'duplicate key');d[k]=v
        return d
    return json.loads(raw,object_pairs_hook=pairs,parse_constant=lambda value:(_ for _ in ()).throw(ValueError('nonfinite JSON')))
def canonical(value):
    if type(value) is Q:return f'{value.numerator}/{value.denominator}'
    if type(value) in (list,tuple):return [canonical(v) for v in value]
    if type(value) is dict:return {k:canonical(v) for k,v in value.items()}
    return value
def extent(interval):return max(abs(interval[0]),abs(interval[1]))
def support(coefficients,box,constant=Q(0)):
    lows=[];highs=[]
    for coefficient,interval in zip(coefficients,box):
        endpoints=[coefficient*x for x in interval]
        lows.append(min(endpoints));highs.append(max(endpoints))
    return (constant+sum(lows,Q(0)),constant+sum(highs,Q(0)))
def overlap(a,b):
    result=(max(a[0],b[0]),min(a[1],b[1]))
    require(result[0]<=result[1],'empty proved interval intersection')
    return result
def exact(value):
    require(type(value) is str and re.fullmatch(r'-?(?:0|[1-9][0-9]*)/[1-9][0-9]*',value) is not None,
            'canonical exact numerator/denominator string required')
    q=Q(value);require(value==f'{q.numerator}/{q.denominator}','noncanonical exact ratio')
    return q

def represented_binary80(value,negative_zero=False):
    require(type(value) is Q and type(negative_zero) is bool,'exact exported scalar/sign required')
    require(not negative_zero or value==0,'negative zero flag on nonzero value')
    if value==0:return
    numerator=abs(value.numerator);denominator=value.denominator
    require(denominator&(denominator-1)==0,'unnormalized represented scalar is not dyadic')
    # Separate the integer's trailing powers of two, then bound the odd
    # significand and exponent directly. No binary80 bytes are read or decoded.
    trailing=(numerator & -numerator).bit_length()-1
    odd=numerator>>trailing;exponent=trailing-(denominator.bit_length()-1)
    require(odd.bit_length()<=64 and exponent>=-16445 and odd.bit_length()-1+exponent<=16383,
            'export is outside finite canonical binary80 representability')

def exact_pair(value):
    require(type(value) is list and len(value)==2,'complex pair required')
    return tuple(exact(x) for x in value)

def empty():
    return {'count':0,'last':Q(0),'sums':{k:Q(0) for k in ('weight_sum','weighted_U_L1_upper',
      'weighted_W_L1_upper','R_triangle_upper','P_triangle_upper','canonical_drift_weighted_signed',
      'canonical_drift_weighted_abs')},'signed_intervals':{'R_error':(Q(0),Q(0)),'P_error':(Q(0),Q(0))},
      'maxima':{k:Q(0) for k in ('U_L1_upper','W_L1_upper','first_order_wronskian_drift_abs')}}

def run(directory,outname):
    started=time.monotonic();directory=Path(directory)
    data=load((directory/'DATA.json').read_bytes())
    require(data['manufactured'] is False and data['scope']=='SAVED_ENDPOINT_MINUS_EXACT_FLOW_FROM_EXACT_REPRESENTED_INCOMING_STATE',
            'actual exact represented-anchor endpoint scope required')
    summaries={x['case']:x for x in data['rows']}
    require(type(data['rows']) is list and len(data['rows'])==len(summaries)==24,'all24 unique aggregates required')
    capsules=[s+'/'+g for s in ('positive_B','signed_uB') for g in ('coarse','fine')]
    capsule_index=index=nodes=comparisons=checked=0;sign_flags=negative_zeros=stored_components=represented_scalars=0
    decoded_bytes=0;computed_rows=[];last_k=Q(0)
    totals={t:empty() for t in (Q(-4),Q(-7,2))};radius={'U':Q(0),'W':Q(0)}
    with gzip.open(directory/'NODE_CERTIFICATES.jsonl.gz','rb') as stream:
        while raw:=stream.readline(65537):
            decoded_bytes+=len(raw)
            require(len(raw)<=65536 and raw.endswith(b'\n') and decoded_bytes<=512*(1<<20),'bounded complete node line required')
            row=load(raw);require(capsule_index<4,'extra capsule')
            capsule=capsules[capsule_index];source,grid=capsule.split('/')
            n=8192 if grid=='coarse' else 16384;plan={n//4:64,n//2:128,n:256}
            k,weight=exact(row['k']),exact(row['weight'])
            require(row['capsule']==capsule and type(row['index']) is int and row['index']==index and
                    last_k<k<256 and weight>0,'capsule/node ordering or momentum/weight invalid')
            last_k=k
            for count,cutoff in plan.items():
                require((index<count and k<cutoff) or (index>=count and k>cutoff),'exact cutoff/prefix mask differs')
            flags=row['negative_zero_flags']
            require(type(flags) is list and len(flags)==14 and all(type(flag) is bool for flag in flags) and
                    flags[0] is False and flags[1] is False,'14 binary80 sign flags required')
            represented_binary80(k,flags[0]);represented_binary80(weight,flags[1]);represented_scalars+=2
            sign_flags+=14;negative_zeros+=sum(flags)
            incoming_u,incoming_w=exact_pair(row['incoming_U']),exact_pair(row['incoming_W'])
            normalized=[*incoming_u,*incoming_w]
            require(set(row['endpoints'])=={'-4','-7/2'},'two fixed endpoints required')
            for t in (Q(-4),Q(-7,2)):
                item=row['endpoints'][str(t)]
                normalized.extend((*exact_pair(item['saved_U']),*exact_pair(item['saved_W'])))
            for component,flag in zip(normalized,flags[2:]):
                represented_binary80(component*EPS,flag)
                stored_components+=1;represented_scalars+=1
            for t in (Q(-4),Q(-7,2)):
                item=row['endpoints'][str(t)];saved_u=exact_pair(item['saved_U']);saved_w=exact_pair(item['saved_W'])
                require(set(item)=={'saved_U','saved_W','target_dyadic96','delta_U','delta_W','canonical_drift'},'endpoint fields differ')
                require(set(item['target_dyadic96'])=={'U','W'},'target fields differ')
                boxes={}
                for part,saved in (('U',saved_u),('W',saved_w)):
                    require(set(item['target_dyadic96'][part])=={'real','imag'},'target axes differ')
                    boxes[part]=[];width=0
                    for j,axis in enumerate(('real','imag')):
                        lo,hi=item['target_dyadic96'][part][axis]
                        require(type(lo) is int and type(hi) is int and lo<=hi,'bad interval')
                        boxes[part].append((saved[j]-Q(hi,SCALE),saved[j]-Q(lo,SCALE)))
                        width+=hi-lo
                    require(Q(width,2*SCALE)<=Q(1,10**18),'complete per-node target width gate fails')
                    radius[part]=max(radius[part],Q(width,2*SCALE))
                    require(item['delta_'+part]==canonical(boxes[part]),'saved-minus-target serialization differs')
                du,dw=boxes['U'],boxes['W'];L=-1/t;measure=weight*k*k/(2*PI*PI)
                c=(saved_u[0]-incoming_u[0])-(saved_w[1]-incoming_w[1])/(2*k)
                require(exact(item['canonical_drift'])==c,'canonical drift differs')
                overlap(du[0],support((1/(2*k),),(dw[1],),c))
                A=(k+3*L*L/(2*k),k/3-L*L/(2*k));B=-L/(2*k)
                direct=[support((measure*a,measure*B,-measure/2),(du[0],dw[0],dw[1])) for a in A]
                correlated=[support((measure*B,measure*(a/(2*k)-Q(1,2))),(dw[0],dw[1]),measure*a*c) for a in A]
                r,p=(overlap(a,b) for a,b in zip(direct,correlated))
                upperu=sum(map(extent,du),Q(0));upperw=sum(map(extent,dw),Q(0))
                values={'weight_sum':weight,'weighted_U_L1_upper':measure*upperu,'weighted_W_L1_upper':measure*upperw,
                        'R_triangle_upper':extent(r),'P_triangle_upper':extent(p),
                        'canonical_drift_weighted_signed':measure*c,'canonical_drift_weighted_abs':measure*abs(c)}
                acc=totals[t];acc['count']+=1;acc['last']=k
                for key,value in values.items():acc['sums'][key]+=value
                for key,value in (('R_error',r),('P_error',p)):
                    acc['signed_intervals'][key]=tuple(a+b for a,b in zip(acc['signed_intervals'][key],value))
                for key,value in (('U_L1_upper',upperu),('W_L1_upper',upperw),('first_order_wronskian_drift_abs',2*EPS*abs(c))):
                    acc['maxima'][key]=max(acc['maxima'][key],value)
                comparisons+=1
            nodes+=1;index+=1
            if index in plan:
                for t,acc in totals.items():
                    expected=summaries[f'{capsule}/{plan[index]}/{t}']
                    require(expected['node_count']==acc['count'] and Q(expected['last_k'])==acc['last'],'prefix counter differs')
                    for key in ('sums','signed_intervals','maxima'):
                        require(expected[key]==canonical(acc[key]),f'aggregate differs: {expected["case"]}/{key}')
                    require(expected['capsule']==capsule and expected['K']==str(plan[index]) and expected['endpoint']==str(t) and
                            expected['incoming_state_error']=='ZERO_BY_EXACT_REPRESENTED_ANCHOR_DEFINITION' and
                            expected['contact_difference']=='EXACT_ZERO_UNDER_IDENTICAL_COMPLETE_CONTACTS' and
                            expected['full_momentum_and_continuous_time_integrals']=='UNRESOLVED','aggregate scope differs')
                    computed_rows.append({'case':expected['case'],'nodes':acc['count'],'last_k':canonical(acc['last']),
                                          **{key:canonical(acc[key]) for key in ('sums','signed_intervals','maxima')}})
                    checked+=1
            if index==n:
                capsule_index+=1;index=0;last_k=Q(0);totals={t:empty() for t in (Q(-4),Q(-7,2))}
    require(capsule_index==4 and index==0 and nodes==49152 and comparisons==98304 and checked==24,'incomplete fixed universe')
    require(max(radius.values())<=Q(1,10**18),'export gate fails')
    summary=load((directory/'SCIENCE_SUMMARY.json').read_bytes())
    require(data['nodes']==nodes and data['endpoint_comparisons']==comparisons and summary['nodes']==nodes and
            summary['endpoint_comparisons']==comparisons and summary['decoded_real_slots']==sign_flags+24==688152 and
            summary['negative_zero_flags']==negative_zeros and summary['maximum_complete_endpoint_export_L1_radius']==canonical(radius),
            'actual scalar/sign/width coverage differs')
    require(summary['status']=='FINITE_RETAINED_ENDPOINT_ERRORS_ENCLOSED' and summary['selected_arrays']==36 and
            summary['source_callbacks']==128 and summary['incoming_error_counted'] is False and
            summary['prior_metric_calibration']=='FAIL_UNCHANGED' and
            summary['original_unsaved_continuous_solver_trajectory']=='UNRECOVERABLE_FROM_RETAINED_MODE_SNAPSHOTS' and
            summary['continuous_momentum_pressure_contact_time_UV_certificate']=='UNRESOLVED','scientific scope changed')
    result={'status':'PASS_INDEPENDENT_ACTUAL_EXPORTED_ALLNODE_ARITHMETIC','nodes':nodes,'comparisons':comparisons,
      'prefixes':checked,'normalized_stored_scalar_components_checked':stored_components,
      'exported_scalar_binary80_representability_checks':represented_scalars,'negative_zero_flags_checked':sign_flags,
      'negative_zero_flags_true':negative_zeros,'producer_reported_total_decoded_real_slots':summary['decoded_real_slots'],
      'observation_time_values_exported_for_independent_check':False,'maximum_radius':canonical(radius),
      'input_output_files':{name:hashlib.sha256((directory/name).read_bytes()).hexdigest()
                           for name in ('DATA.json','NODE_CERTIFICATES.jsonl.gz','SCIENCE_SUMMARY.json')},
      'retained_arrays_opened':0,'physical_source_calls':0,'production_modules_imported':0,
      'original_decoded_value_identity_independently_verified':False,
      'target_flow_enclosures_independently_reconstructed':False,
      'truth_basis':'Exact serialized arithmetic checked independently; original values and target containment rely on the public byte-pinned decoder/worker and reviewed source/flow proofs.',
      'inverse_normalization_check':'Each exported normalized component times the exact represented epsilon is finite canonical binary80 representable; this does not independently authenticate its archived value.',
      'wall_seconds':time.monotonic()-started,'computed_finite_aggregates':computed_rows}
    with (BASE/outname).open('x') as out:json.dump(result,out,sort_keys=True,indent=2);out.write('\n')
    print(json.dumps({key:value for key,value in result.items() if key!='computed_finite_aggregates'},sort_keys=True))

if __name__=='__main__':run(*sys.argv[1:])
