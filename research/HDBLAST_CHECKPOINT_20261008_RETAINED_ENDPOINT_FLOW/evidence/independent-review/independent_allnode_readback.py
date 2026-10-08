"""Read manufactured node certificates with separately written exact operators.

Imports no production module, Arb, source, or retained decoder. The fixture's
closed-form stored values validate member selection and epsilon normalization.
This checks serialized arithmetic, not the truth of a target-flow enclosure.
"""
from fractions import Fraction as Q
from pathlib import Path
import gzip
import hashlib
import json
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
    return json.loads(raw,object_pairs_hook=pairs)
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
def state(sid,step,index,divisor):
    return (Q((sid+1)*step*(index%17+1),divisor)/EPS,
            Q(-(sid+1)*step*(index%13+1),divisor)/EPS)
def empty():
    return {'count':0,'last':Q(0),'sums':{k:Q(0) for k in ('weight_sum','weighted_U_L1_upper',
      'weighted_W_L1_upper','R_triangle_upper','P_triangle_upper','canonical_drift_weighted_signed',
      'canonical_drift_weighted_abs')},'signed_intervals':{'R_error':(Q(0),Q(0)),'P_error':(Q(0),Q(0))},
      'maxima':{k:Q(0) for k in ('U_L1_upper','W_L1_upper','first_order_wronskian_drift_abs')}}

def run(directory,outname):
    started=time.monotonic();directory=Path(directory)
    data=load((directory/'DATA.json').read_bytes())
    require(data['manufactured'] is True,'this control reads fabricated output only')
    summaries={x['case']:x for x in data['rows']}
    require(len(summaries)==24,'all24 aggregates required')
    capsules=[s+'/'+g for s in ('positive_B','signed_uB') for g in ('coarse','fine')]
    capsule_index=index=nodes=comparisons=checked=0;sign_flags=0
    totals={t:empty() for t in (Q(-4),Q(-7,2))};radius={'U':Q(0),'W':Q(0)}
    with gzip.open(directory/'NODE_CERTIFICATES.jsonl.gz','rb') as stream:
        for raw in stream:
            row=load(raw);require(capsule_index<4,'extra capsule')
            capsule=capsules[capsule_index];source,grid=capsule.split('/')
            sid=0 if source=='positive_B' else 1;den=64 if grid=='coarse' else 128
            n=8192 if grid=='coarse' else 16384;plan={n//4:64,n//2:128,n:256}
            k=Q(2*index+1,den);weight=Q(2,den)
            require(row['capsule']==capsule and row['index']==index and Q(row['k'])==k and Q(row['weight'])==weight,
                    'fixture momentum/weight/node ordering differs')
            require(type(row['negative_zero_flags']) is list and len(row['negative_zero_flags'])==14 and
                    all(flag is False for flag in row['negative_zero_flags']),
                    'manufactured nonzero scalar carries an incorrect sign flag')
            sign_flags+=14
            incoming_u=state(sid,1,index,1<<30);incoming_w=state(sid,1,index,1<<31)
            require(tuple(map(Q,row['incoming_U']))==incoming_u and tuple(map(Q,row['incoming_W']))==incoming_w,
                    'incoming normalization/selection differs')
            for step,t in ((2,Q(-4)),(3,Q(-7,2))):
                item=row['endpoints'][str(t)];saved_u=state(sid,step,index,1<<30);saved_w=state(sid,step,index,1<<31)
                require(tuple(map(Q,item['saved_U']))==saved_u and tuple(map(Q,item['saved_W']))==saved_w,
                        'later selection or normalization differs')
                boxes={}
                for part,saved in (('U',saved_u),('W',saved_w)):
                    boxes[part]=[];width=0
                    for j,axis in enumerate(('real','imag')):
                        lo,hi=item['target_dyadic96'][part][axis]
                        require(type(lo) is int and type(hi) is int and lo<=hi,'bad interval')
                        boxes[part].append((saved[j]-Q(hi,SCALE),saved[j]-Q(lo,SCALE)))
                        width+=hi-lo
                    radius[part]=max(radius[part],Q(width,2*SCALE))
                    require(item['delta_'+part]==canonical(boxes[part]),'saved-minus-target serialization differs')
                du,dw=boxes['U'],boxes['W'];L=-1/t;measure=weight*k*k/(2*PI*PI)
                c=(saved_u[0]-incoming_u[0])-(saved_w[1]-incoming_w[1])/(2*k)
                require(Q(item['canonical_drift'])==c,'canonical drift differs')
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
                    checked+=1
            if index==n:
                capsule_index+=1;index=0;totals={t:empty() for t in (Q(-4),Q(-7,2))}
    require(capsule_index==4 and index==0 and nodes==49152 and comparisons==98304 and checked==24,'incomplete fixed universe')
    require(max(radius.values())<=Q(1,10**18),'export gate fails')
    summary=load((directory/'SCIENCE_SUMMARY.json').read_bytes())
    require(summary['decoded_real_slots']==sign_flags+24==688152 and summary['negative_zero_flags']==0,
            'complete manufactured scalar/sign coverage differs')
    result={'status':'PASS_INDEPENDENT_MANUFACTURED_ALLNODE_READBACK','nodes':nodes,'comparisons':comparisons,
      'prefixes':checked,'normalized_stored_scalar_components_checked':nodes*12,
      'negative_zero_flags_checked':sign_flags,'decoded_real_slots':sign_flags+24,'maximum_radius':canonical(radius),
      'input_output_files':{name:hashlib.sha256((directory/name).read_bytes()).hexdigest()
                           for name in ('DATA.json','NODE_CERTIFICATES.jsonl.gz')},
      'retained_arrays_opened':0,'physical_source_calls':0,'production_modules_imported':0,
      'wall_seconds':time.monotonic()-started}
    with (BASE/outname).open('x') as out:json.dump(result,out,sort_keys=True,indent=2);out.write('\n')
    print(json.dumps(result,sort_keys=True))

if __name__=='__main__':run(*sys.argv[1:])
