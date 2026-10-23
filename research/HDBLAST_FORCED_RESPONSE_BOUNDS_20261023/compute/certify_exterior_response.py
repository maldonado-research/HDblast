#!/usr/bin/env python3
"""Whole-cell interval enclosure of an actual forced-response mean on one cone.

This does not compute the global waveform or replace the selected forcing by
its comparison source. All23 bulk sectors enter the norm bound for the full
coupling correction; the shell contribution vanishes by domain of dependence.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import sys
import time

MANIFEST_SHA = '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2'
X0 = Q(1,4)
R0 = Q(5,4)
SIGMA = Q(203,1000000)

def need(ok, message):
    if not ok:
        raise ValueError(message)

def pin_source(root):
    raw = (root/'MANIFEST.json').read_bytes()
    need(hashlib.sha256(raw).hexdigest() == MANIFEST_SHA, 'source manifest pin')
    manifest = json.loads(raw)
    for row in manifest['files']:
        name = row['path']
        need(not Path(name).is_absolute() and '..' not in Path(name).parts, 'source member path')
        data = (root/name).read_bytes()
        need(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'], 'source member '+name)
    return {'manifest_sha256':MANIFEST_SHA,'files_authenticated':len(manifest['files'])+1}

def sqrt_upper(x):
    z = math.isqrt(x.numerator//x.denominator)+1
    need(Q(z*z) >= x, 'square root bound')
    return Q(z)

def forcing_bound(caps, angular):
    total = Q(0)
    blocks = {}
    for field,basis,multiplicity in [('yy',angular.SCALARS,1),('y0',angular.SCALARS,2),('00',angular.SCALARS,1),('ya',angular.GRADS,2),('0a',angular.GRADS,2),('ab',angular.TENSORS,1),('scalar',angular.SCALARS,1)]:
        prefix = 'scalar_forcing' if field == 'scalar' else 'metric_forcing/'+field
        c = [caps[prefix+'/'+str(i)] for i in range(len(basis))]
        gram = [[angular.inner(a,b) for b in basis] for a in basis]
        square = multiplicity*sum((abs(gram[i][k])*c[i]*c[k] for i in range(len(c)) for k in range(len(c))),Q(0))
        total += square
        blocks[field] = {'multiplicity':multiplicity,'gram':[[str(v) for v in r] for r in gram],'squared_upper':str(square)}
    return {'F_H_upper':str(sqrt_upper(total)),'squared_upper':str(total),'blocks':blocks}

def run(source, n):
    need(n>=1, 'positive partition size')
    own_pin=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    authentication = pin_source(source)
    sys.path.insert(0,str(source/'numeric'))
    sys.path.insert(0,str(source/'projection'))
    import actual_first_jets as jets
    import fast_dyadic
    import source_coefficients as coefficients
    import angular_generic as angular
    need(tuple(angular.inner(s,angular.SCALARS[0]) for s in angular.SCALARS)==(Q(1),Q(0),Q(0)), 'constant scalar slot equals probability sphere mean')
    accepted=json.loads((source/'numeric/ACTUAL_SOURCES_NORMAL.json').read_text())
    global_F=Q(accepted['global_cover']['summary']['positive_bundle_bulk_forcing_H_bound']['F_H_upper'])
    need(global_F==17387513156573,'accepted global full23 forcing norm')
    maps=source/'projection/GEOMETRIC_COEFFICIENT_MAPS.json'
    jets.inputs()
    h=X0/n
    rows=[]
    total_area=Q(0)
    integral_lo=integral_hi=Q(0)
    caps={}
    started=time.monotonic()
    for i in range(n):
        for k in range(i,2*n-i):
            xlo,xhi=i*h,(i+1)*h
            rlo,rhi=R0-X0+k*h,R0-X0+(k+1)*h
            area=h*h/Q(2) if k in (i,2*n-i-1) else h*h
            need(xlo>=0 and xhi<=X0 and rlo-xhi>=Q(3,4), 'controlled exterior bounding box')
            m,s,b=jets.actual_jet(jets.I(xlo,xhi),jets.I(rlo,rhi),mode='actual')
            convert=lambda d:{key:fast_dyadic.D(value) if hasattr(value,'lo') and hasattr(value,'hi') else value for key,value in d.items()}
            values=coefficients.evaluate_sources(convert(m),convert(s),convert(b),names=('metric_forcing','scalar_forcing'),path=maps)
            normalized={}
            for name,fields in values.items():
                fields=fields if isinstance(fields,dict) else {'':fields}
                for field,seq in fields.items():
                    prefix=name+('/'+field if field else '')
                    normalized[prefix]=[v.row() for v in seq]
                    for index,value in enumerate(seq):
                        key=prefix+'/'+str(index)
                        caps[key]=max(caps.get(key,Q(0)),value.mag())
            lo,hi=values['metric_forcing']['yy'][0].lo,values['metric_forcing']['yy'][0].hi
            integral_lo+=area*lo/Q(2)
            integral_hi+=area*hi/Q(2)
            total_area+=area
            rows.append({'i':i,'k':k,'x':[str(xlo),str(xhi)],'r':[str(rlo),str(rhi)],'intersection':'abs(r-5/4)<=1/4-x','exact_area':str(area),'sources':normalized})
        print(json.dumps({'row':i+1,'rows':n,'cells':len(rows),'wall_seconds':round(time.monotonic()-started,3)}),flush=True)
    need(len(rows)==n*(n+1) and total_area==X0*X0,'exact gap-free cone partition and area')
    hb=forcing_bound(caps,angular)
    Fcap=min(Q(hb['F_H_upper']),global_F)
    qC=X0*SIGMA
    BC=X0*Fcap/(1-qC)
    # W-S = K(A W_x+B W_r+D0 W); value gain is x0^2/2.
    value_error=X0*X0*SIGMA*BC/Q(2)
    response=[integral_lo-value_error,integral_hi+value_error]
    epsilon=Q(1,1000000)
    tplus=X0+2*epsilon
    need(tplus<3 and R0-tplus>0 and R0+tplus<8,'enlarged exterior cone in controlled domain')
    qplus=tplus*SIGMA
    Bplus=tplus*global_F/(1-qplus)
    variation=2*epsilon*Bplus
    rectangle_response=[response[0]-variation,response[1]+variation]
    source_post=pin_source(source)
    need(authentication==source_post,'source bytes unchanged')
    need(own_pin==hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'evaluator bytes unchanged')
    return {'status':'CERTIFIED_CAUSAL_RESTRICTED_ACTUAL_FORCED_MEAN_INTERVAL','source_authentication':authentication,'evaluator_sha256':own_pin,
        'domain':{'apex_x':str(X0),'apex_r':str(R0),'cone':'0<=x<=1/4,abs(r-5/4)<=1/4-x','initial_interval':['1','3/2'],'min_r_minus_x_on_cone':'1','backward_closed':True,'shell_and_reflection_reached':False},
        'partition':{'n':n,'h':str(h),'cells':len(rows),'area':str(total_area),'proof':'For x-row i, k=i..2n-i-1; two outer cells each have clipped triangular area h^2/2 and all inner cells h^2. Rectangles bound the true cone intersections; interiors disjoint and union equals cone.'},
        'actual_source_cell_enclosures':rows,'forcing_H_bound':hb,'accepted_global_F_H_upper':str(global_F),'F_H_upper_used':str(Fcap),
        'forcing_coefficient_caps':{k:str(v) for k,v in caps.items()},
        'full_coupling':{'A_upper':'1/10000','B_upper':'1/10000','D0_upper':'1/1000000','sigma_upper':str(SIGMA),'K_X_operator_upper':str(X0),'q_C_upper':str(qC),'W_X_C_upper':str(BC),'all23_bulk_sectors_retained':True,'Robin_term_zero_only_by_backward_domain_of_dependence':True},
        'evaluated_Duhamel_mean_interval':[str(integral_lo),str(integral_hi)],
        'coupling_point_value_error_upper':str(value_error),
        'actual_W_yy_probability_mean_interval':[str(v) for v in response],
        'strict_positive':response[0]>0,'strict_negative':response[1]<0,
        'continuous_rectangle_corollary':{'x':[str(X0-epsilon),str(X0+epsilon)],'r':[str(R0-epsilon),str(R0+epsilon)],
            'epsilon':str(epsilon),'enlarged_cone_apex_x':str(tplus),'enlarged_cone_apex_r':str(R0),
            'enlarged_cone_min_r_minus_x':str(R0-tplus),'global_F_H_upper_used':str(global_F),
            'q_enlarged_cone_upper':str(qplus),'W_gradient_H_upper':str(Bplus),
            'point_to_rectangle_l1_distance_upper':str(2*epsilon),'mean_variation_upper':str(variation),
            'actual_W_yy_probability_mean_interval':[str(v) for v in rectangle_response],
            'strict_positive':rectangle_response[0]>0,
            'proof':'The enlarged backward cone contains the complete rectangle. The actual W has both mean derivatives bounded by Bplus there. Integrate the gradient along the straight segment inside that convex cone; L1 displacement is at most2epsilon.'},
        'normalization':'Actual metric forcing=-2ell^2rho^2QEred. yy scalar basis constant coefficient equals probability S3 mean. K has factor1/2. Sources and area use scaled x,r.',
        'limits':['A pointwise enclosure of the full auxiliary forced increment, not a global evaluated waveform or exact physical quadratic solution by itself.','The full23-state coupling is enclosed analytically; its particular mean is not assumed to decouple.','No physical force, invariant observable, nonlinear family, cosmological-origin or novelty conclusion.','All original actual C3 uncertainty and source interval dependencies are retained; no comparator substitution or source derivative bounds.']}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--source-root',type=Path,required=True)
    p.add_argument('--n',type=int,default=4)
    p.add_argument('--require-positive',action='store_true')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    need(not args.output.exists(),'refuse overwrite output')
    result=run(args.source_root.resolve(),args.n)
    need(not args.require_positive or result['strict_positive'],'mean sign not separated from zero')
    need(not args.require_positive or result['continuous_rectangle_corollary']['strict_positive'],'rectangle mean sign not separated from zero')
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['status','evaluated_Duhamel_mean_interval','coupling_point_value_error_upper','actual_W_yy_probability_mean_interval','strict_positive','strict_negative']},indent=2))

if __name__=='__main__':
    main()
