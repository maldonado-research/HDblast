#!/usr/bin/env python3
"""Fabricated mathematical controls and full-universe resource rehearsal."""
from fractions import Fraction as Q
from pathlib import Path
import argparse,json,time,resource,copy
from flint import arb,acb,ctx
import entire_target_engine as e

def qr(q):return str(q.numerator)+'/'+str(q.denominator)

def contains(rect,z):
    for key,part in (('real',z.real),('imag',z.imag)):
        lo,hi=(Q(i,e.EXPORT_SCALE) for i in rect[key])
        if not (lo<=e.exact_q(part.lower()) and e.exact_q(part.upper())<=hi):return False
    return True

def fixture(sign=1,scale=Q(1),constant=False,highest=False):
    rows=[]
    for g in e.geometry():
        co=[scale*Q(sign**j,(j+1)*2**j)/g['half']**j for j in range(113)]
        if constant:co=[scale]+[Q(0)]*112
        if highest:co=[Q(0)]*112+[scale/g['half']**112]
        rows.append({'center':g['center'],'half':g['half'],
                     'coefficients':co,'uniform_error':Q(0)})
    return rows

def controls():
    checks=[];mutants=[]
    def check(ok,name):
        e.require(ok,name);checks.append(name)
    def reject(name,fn):
        try:fn()
        except ValueError:mutants.append(name)
        else:raise ValueError('Accepted mutant: '+name)
    ctx.prec=e.BITS
    constant=e.EntireTarget(fixture(constant=True));d=Q(63,128)
    for k in (Q(0),Q(1,2**40),Q(1,4),Q(1),Q(16),Q(128),Q(256)):
        result=constant.evaluate(k)
        if not k:
            u=acb(e.av(-d*d/2));w=acb(e.av(-d))
        else:
            z=acb(0,e.av(2*k));phase=(z*e.av(d)).exp()
            w=-(phase-1)/z;u=-(phase-1-z*e.av(d))/(z*z)
        check(contains(result['U'],u),'constant finite-source exact U at '+qr(k))
        check(contains(result['W'],w),'constant finite-source exact W at '+qr(k))
    highest=e.EntireTarget(fixture(highest=True));r=highest.evaluate(Q(0))
    wm=-d/113
    um=-sum((g['half']*2*Q(1,113)*(-Q(9,2)-g['center']) for g in e.geometry()),Q(0))
    check(contains(r['W'],acb(e.av(wm))),'highest coefficient112 exact W moment retained')
    check(contains(r['U'],acb(e.av(um))),'highest coefficient112 exact U moment retained')
    coords=e.stored_state_coordinates(Q(1,8),(Q(5),Q(-7)),(Q(11),Q(13)),Q(1,16))
    check(coords['c']==Q(-752),'nonzero exact stored canonical defect retained')
    check(coords['first_order_wronskian_defect']==Q(-94),'epsilon scaling retained in Wronskian defect')
    check(coords['kA_saved']==(Q(104),Q(-88)),'stable kA sign and half factor')
    rows=fixture();bad=copy.deepcopy(rows);bad[0]['coefficients'].pop()
    reject('degree111 source model',lambda:e.EntireTarget(bad))
    bad2=copy.deepcopy(rows);bad2[0]['uniform_error']=Q(-1)
    reject('negative source certificate radius',lambda:e.EntireTarget(bad2))
    bad3=copy.deepcopy(rows);bad3[0]['center']+=Q(1,100)
    reject('changed exact panel geometry',lambda:e.EntireTarget(bad3))
    reject('missing source parent',lambda:e.EntireTarget(rows[:-1]))
    reject('node exceeds fixed momentum band',lambda:constant.evaluate(Q(257)))
    reject('nonrational node',lambda:constant.evaluate(1.0))
    reject('gate cannot be relaxed by export',lambda:e.export_rectangle(acb(arb(0,1)),Q(0)))
    reject('negative analytic tail',lambda:e.export_rectangle(acb(0),Q(-1)))
    ctx.prec=511
    reject('arithmetic precision change',lambda:constant.evaluate(Q(1)))
    ctx.prec=e.BITS
    return {'passed_controls':checks,'passed_control_count':len(checks),
            'rejected_mutations':mutants,'rejected_mutation_count':len(mutants)}

def benchmark():
    started=time.monotonic();summaries=[];count=0;maxr=Q(0)
    for source,sign in [('fabricated_positive',1),('fabricated_alternating',-1)]:
        t0=time.monotonic();model=e.EntireTarget(fixture(sign,Q(2**24)))
        t1=time.monotonic()
        for n in (8192,16384):
            for j in range(n):
                k=Q(256*(2*j+1),2*n)
                out=model.evaluate(k)
                for name in ('U','W'):
                    rect=out[name]
                    rad=Q(sum(rect[a][1]-rect[a][0] for a in ('real','imag')),2*e.EXPORT_SCALE)
                    maxr=max(maxr,rad)
                count+=1
        t2=time.monotonic()
        summaries.append({'source':source,'source_L1_majorant':qr(model.source_l1),
                          'build_seconds':t1-t0,'evaluate_and_export_seconds':t2-t1})
    e.require(count==49152,'All fabricated nodes must execute')
    return {'node_count':count,'complex_rectangles':2*count,
            'largest_exact_exported_L1_radius':qr(maxr),
            'source_summaries':summaries,'wall_seconds':time.monotonic()-started,
            'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();e.require(not args.output.exists(),'Fresh output required')
    result={'scope':'FABRICATED_ONLY','physical_source_callbacks':0,'saved_array_decodes':0,
            'bits':e.BITS,'momentum_degree':e.MOMENTUM_DEGREE,'source_degree':e.SOURCE_DEGREE,
            'export_bits':e.EXPORT_BITS,'fixed_exported_L1_gate':qr(e.EXPORTED_L1_GATE),
            'controls':controls(),'benchmark':benchmark()}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS_FABRICATED_CONTROLS_AND_49152_NODES',
                      'output':str(args.output),'wall_seconds':result['benchmark']['wall_seconds']}))
