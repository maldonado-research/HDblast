#!/usr/bin/env python3
"""Independent saved-result audit; never import a source or producer.

Only standard-library exact rational arithmetic and saved execution artifacts
are used. This augments, rather than replaces, the frozen strict reviewer.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import isfinite
from pathlib import Path
import re

SOURCES=('positive_B','signed_uB')
KS=(Q(0),Q(1,2**40),Q(1,2**12),Q(1,4),Q(1),Q(16),Q(64),Q(128),Q(256))
CHANNELS=('M0','Mexp','Mu')
GATE=Q(1,10**26)
UNIT=Q(1,2**512)
SOURCE_TAIL=Q(1,15*2**90)


def require(value,message):
    if not value:raise ValueError(message)


def load(raw):
    def unique(items):
        d={}
        for k,v in items:
            require(k not in d,'Duplicate JSON key')
            d[k]=v
        return d
    def invalid(value):raise ValueError('Invalid JSON constant '+value)
    return json.loads(raw,object_pairs_hook=unique,parse_constant=invalid)


def rational(value):
    require(type(value) is str and len(value)<=4096 and
            re.fullmatch(r'-?(?:0|[1-9][0-9]*)/[1-9][0-9]*',value),
            'Bounded canonical rational required')
    q=Q(value)
    require(value==canonical(q),'Noncanonical or nonreduced rational')
    return q


def canonical(q):return str(q.numerator)+'/'+str(q.denominator)


def rect(value):
    require(type(value) is dict and set(value)=={'real','imag'},'Exact rectangle fields')
    out={}
    for axis in ('real','imag'):
        r=value[axis]
        require(type(r) is dict and set(r)=={'lo','hi'},'Exact endpoint fields')
        lo,hi=rational(r['lo']),rational(r['hi'])
        require(lo<=hi,'Reversed interval')
        out[axis]=(lo,hi)
    return out


def radius(value):return sum(((hi-lo)/2 for lo,hi in value.values()),Q(0))


def expected_rows(field):
    if field=='panel_rows':
        return [(s,j,k,(Q(-9,2)+Q(j,64),Q(-9,2)+Q(j+1,64)))
                for s in SOURCES for j in range(64) for k in KS]
    if field=='whole_rows':return [(s,None,k,(Q(-9,2),Q(-7,2))) for s in SOURCES for k in KS]
    if field=='source_work_panel_rows':
        return [(s,j,None,(Q(-9,2)+Q(j,64),Q(-9,2)+Q(j+1,64))) for s in SOURCES for j in range(64)]
    return [(s,None,None,(Q(-9,2),Q(-7,2))) for s in SOURCES]


def inspect_payload(p,route):
    require(set(p)=={'schema_version','scope','configuration','archive_arrays_decoded',
                    'panel_rows','whole_rows','source_work_panel_rows','source_work_whole_rows'},
            'Payload field inventory differs')
    require(type(p['schema_version']) is int and p['schema_version']==1 and
            type(p['archive_arrays_decoded']) is int and p['archive_arrays_decoded']==0,
            'Payload schema/archive scope differs')
    require(p['scope']=='UNIFORM_ANALYTIC_MODEL_PLUS_FIXED_RATIONAL_MOMENT_PROBES',
            'Payload scope differs')
    config='PRIMARY_ARB256_SOURCE24_PHASE96' if route=='primary' else 'INDEPENDENT_DYADIC512_SOURCE24_ODE96'
    require(p['configuration']==config,'Method configuration differs')
    result={};maximum=Q(0)
    for field in ('panel_rows','whole_rows','source_work_panel_rows','source_work_whole_rows'):
        expected=expected_rows(field);rows=p[field]
        require(type(rows) is list and len(rows)==len(expected),'Row count differs '+field)
        channels={}
        for index,(row,(source,panel,k,bounds)) in enumerate(zip(rows,expected)):
            metadata={'source','interval'}|({'panel'} if panel is not None else set())|({'momentum'} if k is not None else set())
            science={'moments','total_absolute_radii'} if k is not None else {'moment','total_absolute_radius'}
            require(type(row) is dict and set(row)==metadata|science,'Row fields differ')
            require(row['source']==source and tuple(rational(v) for v in row['interval'])==bounds,
                    'Source/interval membership differs')
            if panel is not None:require(type(row['panel']) is int and row['panel']==panel,'Panel membership/order differs')
            if k is not None:
                require(rational(row['momentum'])==k,'Probe membership/order differs')
                require(set(row['moments'])==set(CHANNELS) and set(row['total_absolute_radii'])==set(CHANNELS),
                        'Moment channel inventory differs')
                items=[(name,row['moments'][name],row['total_absolute_radii'][name]) for name in CHANNELS]
            else:items=[('Lg',row['moment'],row['total_absolute_radius'])]
            for name,value,declared in items:
                r=rect(value);d=rational(declared)
                require(d==radius(r),'Declared radius differs from exported L1 half-width')
                require(d>=0,'Negative declared radius')
                if name in ('M0','Lg') or k==0:
                    require(r['imag']==(Q(0),Q(0)),'Real-target imaginary projection differs')
                channels[(index,name)]=r
                if 'whole' in field:maximum=max(maximum,d)
        result[field]=channels
    for sid,source in enumerate(SOURCES):
        whole_m0=[]
        for kid,k in enumerate(KS):
            w=result['whole_rows'][(sid*9+kid,'M0')]['real']
            whole_m0.append(w)
            local=[result['panel_rows'][(sid*64*9+j*9+kid,'M0')]['real'] for j in range(64)]
            summed=(sum(lo for lo,_ in local),sum(hi for _,hi in local))
            require(max(w[0],summed[0])<=min(w[1],summed[1]),'Whole M0/additive source intervals disjoint')
            if k==0:
                me=result['whole_rows'][(sid*9+kid,'Mexp')]['real']
                require(max(w[0],me[0])<=min(w[1],me[1]),'Zero-phase M0/Mexp identity inconsistent')
        require(max(lo for lo,_ in whole_m0)<=min(hi for _,hi in whole_m0),
                'Whole M0 momentum-independent target inconsistent')
        for j in range(64):
            boxes=[result['panel_rows'][(sid*64*9+j*9+kid,'M0')]['real'] for kid in range(9)]
            require(max(lo for lo,_ in boxes)<=min(hi for _,hi in boxes),'Panel M0 momentum-independent target inconsistent')
            me=result['panel_rows'][(sid*64*9+j*9,'Mexp')]['real'];m0=boxes[0]
            require(max(me[0],m0[0])<=min(me[1],m0[1]),'Panel zero-phase M0/Mexp identity inconsistent')
        work=result['source_work_whole_rows'][(sid,'Lg')]['real']
        local=[result['source_work_panel_rows'][(sid*64+j,'Lg')]['real'] for j in range(64)]
        summed=(sum(lo for lo,_ in local),sum(hi for _,hi in local))
        require(max(work[0],summed[0])<=min(work[1],summed[1]),'Whole/additive Lg intervals inconsistent')
    return result,maximum


def sha(raw):return hashlib.sha256(raw).hexdigest()


def inspect_saved(directory,route):
    raws={}
    for name in ('OUTPUT.json','EXECUTION.json','SOURCE_ATTEMPTS.jsonl','child.log'):
        p=directory/name
        require(p.is_file() and not p.is_symlink(),'Missing/symlink saved artifact '+name)
        raws[name]=p.read_bytes()
    env=load(raws['OUTPUT.json']);outer=load(raws['EXECUTION.json'])
    require(env['route']==outer['route']==route,'Saved route differs')
    require(outer['status']=='PASS_BOUNDED_ROUTE_EXECUTION' and outer['exit_code']==0 and
            outer['timed_out'] is False and outer['fabricated_only'] is False,
            'Completed actual guarded execution required')
    require(type(outer['uid']) is int and outer['uid']>0,'Nonroot execution required')
    require(outer['outer_limits']=={'wall_seconds':120,'peak_rss_kib':131072},'Fixed resource budget differs')
    require(type(outer['wall_seconds']) in (float,int) and isfinite(outer['wall_seconds']) and
            0<outer['wall_seconds']<=120 and type(outer['peak_rss_kib_wait4']) is int and
            0<outer['peak_rss_kib_wait4']<=131072,'Authoritative resource gate failed')
    require(outer['output_sha256']==sha(raws['OUTPUT.json']) and outer['output_bytes']==len(raws['OUTPUT.json']),
            'Saved output pin/size differs')
    require(outer['child_log_sha256']==sha(raws['child.log']) and
            outer['source_attempt_journal_sha256']==sha(raws['SOURCE_ATTEMPTS.jsonl']) and
            outer['source_attempt_journal_bytes']==len(raws['SOURCE_ATTEMPTS.jsonl']),
            'Saved log/journal pin differs')
    chronology=env['execution']['source_chronology']
    require(chronology['route']==route and chronology['source_bundle_constructions']==128 and
            type(chronology['source_bundle_constructions']) is int and
            type(chronology['archive_arrays_decoded']) is int and
            chronology['archive_arrays_decoded']==0,'Actual source/archive counters differ')
    expected=[{'source':s,'center':str(Q(-9,2)+Q(2*j+1,128))} for s in SOURCES for j in range(64)]
    require(chronology['events']==expected,'Exact source event universe/order differs')
    lines=raws['SOURCE_ATTEMPTS.jsonl'].splitlines()
    require(raws['SOURCE_ATTEMPTS.jsonl'].endswith(b'\n') and len(lines)==128,'Persisted journal completeness differs')
    for index,(line,event) in enumerate(zip(lines,expected),1):
        require(load(line)=={'index':index,'route':route,**event},'Persisted permitted attempt differs')
    evidence=env['method_evidence']
    if route=='independent':
        require(type(evidence['physical_source_evaluations']) is int and evidence['physical_source_evaluations']==128 and
                evidence['fabricated_source_provider'] is False and
                type(evidence['retained_arrays_decoded']) is int and type(evidence['primary_helper_imports']) is int and
                evidence['retained_arrays_decoded']==evidence['primary_helper_imports']==0,
                'Independent actual-source flags differ')
        for field,value in [('source_degree',24),('mode_degree',96),('scalar_exp_degree',200),('dyadic_bits',512)]:
            require(type(evidence[field]) is int and evidence[field]==value,'Independent degree/precision differs')
        require(evidence['source_forcing_disk_bound']=='64/1' and evidence['Lg_disk_bound']=='32/1',
                'Independent analytic source majorants differ')
        require(len(evidence['source_model_rows'])==128 and len(evidence['mode_error_rows'])==1152,
                'Independent complete error/source universe differs')
        source_errors={};work_errors={}
        for index,row in enumerate(evidence['source_model_rows']):
            source=SOURCES[index//64];panel=index%64
            require(row['source']==source and type(row['panel']) is int and row['panel']==panel and
                    rational(row['center'])==Q(-9,2)+Q(2*panel+1,128),'Independent source row membership differs')
            for name,tail,dest in [('forcing',SOURCE_TAIL,source_errors),('Lg',SOURCE_TAIL/2,work_errors)]:
                part=row[name];extra=rational(part['coefficient_error_uniform']);full=rational(part['uniform_error'])
                require(rational(part['analytic_tail'])==tail and extra>=0 and full==tail+extra,
                        'Independent source error decomposition differs')
                dest[(source,panel)]=full
        mode_keys={'source_ME_error','source_Mu_error','local_ME_defect_error','local_Mu_defect_error',
                   'prefix_ME_defect_increment','prefix_Mu_defect_increment','inherited_w_radius',
                   'inherited_u_radius','chosen_w_point_shift','chosen_u_point_shift',
                   'w_radius_upward_rounding','u_radius_upward_rounding'}
        mode_rows={}
        expected=((s,k,j) for s in SOURCES for k in KS for j in range(64))
        for row,(source,k,panel) in zip(evidence['mode_error_rows'],expected):
            require(row['source']==source and rational(row['momentum'])==k and type(row['panel']) is int and
                    row['panel']==panel,'Independent defect row order/membership differs')
            for key in mode_keys:require(rational(row[key])>=0,'Negative error evidence')
            for key,target in [('source_ME_error',source_errors[(source,panel)]/64),
                               ('source_Mu_error',source_errors[(source,panel)]/(2*64**2))]:
                reported=rational(row[key]);require(target<=reported<target+UNIT,'Source propagation/outward error differs')
            if panel==0:require(row['inherited_w_radius']==row['inherited_u_radius']=='0/1','Nonzero incoming state uncertainty')
            mode_rows[(source,k,panel)]={key:rational(row[key]) for key in mode_keys}
        for source in SOURCES:
            for k in KS:
                for j in range(64):
                    e=mode_rows[(source,k,j)]
                    upper_w=sum((e[key] for key in ('inherited_w_radius','source_ME_error',
                                    'prefix_ME_defect_increment','chosen_w_point_shift','w_radius_upward_rounding')),Q(0))
                    upper_u=e['inherited_u_radius']+e['inherited_w_radius']/64+sum((e[key] for key in
                        ('source_Mu_error','prefix_Mu_defect_increment','chosen_u_point_shift','u_radius_upward_rounding')),Q(0))
                    if j<63:
                        nxt=mode_rows[(source,k,j+1)]
                        require(nxt['inherited_w_radius']<=upper_w and nxt['inherited_u_radius']<=upper_u,
                                'Absolute state-error carry omits a contribution')
                    else:
                        out=env['payload']['whole_rows'][SOURCES.index(source)*9+KS.index(k)]
                        factor=1 if k==0 else 2
                        require(rational(out['total_absolute_radii']['Mexp'])<=factor*(upper_w+UNIT) and
                                rational(out['total_absolute_radii']['Mu'])<=factor*(upper_u+UNIT),
                                'Whole rectangle error exceeds independently summed propagation budget')
            m0row=env['payload']['whole_rows'][SOURCES.index(source)*9]
            workrow=env['payload']['source_work_whole_rows'][SOURCES.index(source)]
            m0budget=sum((source_errors[(source,j)]/64 for j in range(64)),Q(0))+65*UNIT
            workbudget=sum((work_errors[(source,j)]/64 for j in range(64)),Q(0))+65*UNIT
            require(rational(m0row['total_absolute_radii']['M0'])<=m0budget and
                    rational(workrow['total_absolute_radius'])<=workbudget,'Whole M0/Lg source/error budget differs')
    else:
        model=evidence['analytic_model']
        require(type(model['source_degree']) is int and model['source_degree']==24 and
                type(model['phase_degree']) is int and model['phase_degree']==96 and
                model['source_majorant']=='64/1' and model['Lg_majorant']=='32/1' and
                model['source_complex_disk_radius']=='1/8' and model['panel_half_width']=='1/128' and
                rational(model['source_uniform_tail'])==SOURCE_TAIL and rational(model['Lg_uniform_tail'])==SOURCE_TAIL/2,
                'Primary analytic model parameters differ')
    parsed,maximum=inspect_payload(env['payload'],route)
    return env,outer,parsed,maximum,{name:sha(raw) for name,raw in raws.items()}


def outward_decimal(q,digits=36,lower=True):
    scale=10**digits
    n=q.numerator*scale//q.denominator if lower else -((-q.numerator*scale)//q.denominator)
    sign='-' if n<0 else '';n=abs(n)
    whole,part=divmod(n,scale)
    return sign+str(whole)+'.'+str(part).rjust(digits,'0')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--primary',type=Path,required=True)
    parser.add_argument('--independent',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    require(not args.output.exists(),'Fresh review receipt required')
    p,po,pr,pmax,ppins=inspect_saved(args.primary,'primary')
    i,io,ir,imax,ipins=inspect_saved(args.independent,'independent')
    for key in ('registration_sha256','remote_receipt_sha256'):
        require(po[key]==io[key],'Routes used different registration/readback')
    require(p['execution']['freeze_commit']==i['execution']['freeze_commit'],'Routes used different freeze')
    intersections=[];maximum_hull=Q(0);disjoint=[]
    for field in pr:
        require(set(pr[field])==set(ir[field]),'Cross-route universe differs')
        for key in pr[field]:
            a,b=pr[field][key],ir[field][key];intersection={};hull={}
            for axis in ('real','imag'):
                lo=max(a[axis][0],b[axis][0]);hi=min(a[axis][1],b[axis][1])
                if lo>hi:disjoint.append({'field':field,'row_channel':list(key),'axis':axis})
                intersection[axis]=(lo,hi)
                hull[axis]=(min(a[axis][0],b[axis][0]),max(a[axis][1],b[axis][1]))
            if 'whole' in field:maximum_hull=max(maximum_hull,radius(hull))
            if 'whole' in field:
                intersections.append({'field':field,'row_index':key[0],'channel':key[1],
                    'intersection':{axis:{'lo':canonical(lo),'hi':canonical(hi),
                        'decimal_lo_outward':outward_decimal(lo,lower=True),
                        'decimal_hi_outward':outward_decimal(hi,lower=False)}
                        for axis,(lo,hi) in intersection.items()}})
    passed=not disjoint and max(pmax,imax,maximum_hull)<=GATE
    receipt={'status':'PASS_INDEPENDENT_SAVED_OPERATOR_RESULT_REVIEW' if passed else 'SCIENTIFIC_NEGATIVE_OR_UNRESOLVED',
        'all_component_intersections_nonempty':not disjoint,'disjoint_components':disjoint,
        'maximum_primary_whole_L1_radius':canonical(pmax),
        'maximum_independent_whole_L1_radius':canonical(imax),
        'maximum_pair_whole_hull_L1_radius':canonical(maximum_hull),'gate':canonical(GATE),
        'registration_sha256':po['registration_sha256'],'remote_receipt_sha256':po['remote_receipt_sha256'],
        'freeze_commit':p['execution']['freeze_commit'],
        'artifact_sha256':{'primary':ppins,'independent':ipins},
        'resources':{'primary':{'wall_seconds':po['wall_seconds'],'peak_rss_kib_wait4':po['peak_rss_kib_wait4']},
                     'independent':{'wall_seconds':io['wall_seconds'],'peak_rss_kib_wait4':io['peak_rss_kib_wait4']}},
        'whole_intersections':intersections,'new_source_callbacks_in_review':0,'retained_arrays_decoded_in_review':0,
        'scope':'Uniform central-interval analytic source/phase model bounds plus computed enclosures at nine exact rational moment probes; full continuous-frequency numerical width, pressure/contacts/ledger and physical cause remain unresolved',
        'frozen_proof_and_library_premises':'CONDITIONAL_NOT_REPROVED_BY_SAVED_RESULT_AUDIT',
        'reviewer_source_sha256':sha(Path(__file__).read_bytes())}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:receipt[k] for k in ('status','all_component_intersections_nonempty','maximum_primary_whole_L1_radius','maximum_independent_whole_L1_radius','maximum_pair_whole_hull_L1_radius')}))


if __name__=='__main__':main()
