#!/usr/bin/env python3
"""Own model construction + local/whole operator route, fabricated input only."""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import resource
import sys
import time
import generic_operator as op
import source_models as sm


def require(ok,message):
    if not ok: raise RuntimeError(message)


def formal_checks():
    checks=[]
    for d in ((Q(3,2),Q(-2,7),Q(-1,5)),
              (Q(10001,8192),Q(-29,256),Q(-1)),
              (Q(10003,8192),Q(31,256),Q(-1))):
        inv=sm.reciprocal_coefficients(d,26)
        for n in range(27):
            value=sum((d[j]*inv[n-j] for j in range(min(n,len(d)-1)+1)),Q(0))
            require(value==(1 if n==0 else 0),'formal reciprocal coefficient')
        checks.append('independent exact reciprocal identity: '+str(d))
        f=(Q(-1,5),)+tuple(-a for a in inv[1:])
        beta=sm.normalized_exp_coefficients(f,26)
        for n in range(1,27):
            require(n*beta[n]==sum((j*f[j]*beta[n-j] for j in range(1,n+1)),Q(0)),
                    'formal exponential differential coefficient')
        checks.append('independent exact normalized exponential differential identity: '+str(d))
    for q in (Q(0),Q(-1,5),Q(-3,7),Q(-1)):
        lo,hi=sm.exp_rational_interval(q)
        lo2,hi2=sm.exp_rational_interval(q,220)
        require(lo<=lo2<=hi2<=hi,'independent exact alternating enclosure nesting')
        checks.append('fake scalar exponential enclosure nested: '+str(q))
    try:
        sm.registered_source_model(Q(0),'not_a_registered_source',None)
    except ValueError as e:
        require('freeze' in str(e),'registered source blocked before source evaluation')
    else: raise RuntimeError('unauthenticated registered builder executed')
    checks.append('registered source callback rejected before any source construction')
    return checks


def fabricated_models():
    """Explicit finite polynomials times exp(q); no B or zB callback.

    d/f are intentionally unrelated to the registered source formula. The
    declared fabricated target is this finite polynomial, with a uniform
    extra source uncertainty. There is no unproved analytic truncation.
    """
    start=time.monotonic();models=[[],[]];work_models=[[],[]];half=Q(1,128);stable=hashlib.sha256()
    for source in range(2):
        for panel in range(64):
            d=(Q(10001+panel,8192),Q(((-1)**source)*(29+panel),256),Q(-1,5))
            inv=sm.reciprocal_coefficients(d,26)
            q=-Q(panel+1,300+panel)
            exponent=(q,)+tuple(-a for a in inv[1:])
            beta=sm.normalized_exp_coefficients(exponent,26)
            lo0,hi0=sm.exp_rational_interval(q)
            if source==0:
                gamma=beta[:25]
            else:
                gamma=tuple(Q(panel+1,101)*beta[n]+(beta[n-1] if n else 0)
                            for n in range(25))
            fake_geometry_center=Q(-13,7)
            ls=tuple(-(-Q(1))**n/fake_geometry_center**(n+1) for n in range(25))
            work_gamma=tuple(sum((ls[j]*gamma[n-j] for j in range(n+1)),Q(0))
                             for n in range(25))
            for multiplier,tail,collection,label in (
                (gamma,Q(1,15*2**90),models,'forcing'),
                (work_gamma,Q(1,15*2**91),work_models,'fake_geometry_work')):
                coeff=[];radii=[]
                for a in multiplier:
                    lo,hi=(a*lo0,a*hi0) if a>=0 else (a*hi0,a*lo0)
                    point=sm.dyadic_point((lo+hi)/2)
                    coeff.append(point);radii.append(max(point-lo,hi-point))
                source_error=tail+sum((r*half**n for n,r in enumerate(radii)),Q(0))
                left=op.centered_to_left(coeff,half)
                collection[source].append((left,source_error))
                stable.update(json.dumps({'label':label,'coefficients':[op.rational(a) for a in left],
                                          'source_error':op.rational(source_error)},
                                         sort_keys=True,separators=(',',':')).encode())
    return models,work_models,{'forcing_models':128,'work_models':128,'source_degree':24,'source_coefficient_bits':512,
                   'wall_seconds':time.monotonic()-start,
                   'chosen_model_sha256':stable.hexdigest(),
                   'scope':'fabricated finite polynomial times own enclosed exp(q); no registered source function'}


def fabricated_work_benchmark(models):
    start=time.monotonic();stable=hashlib.sha256();local=0;whole=0;max_radius=Q(0)
    for source in range(2):
        total=Q(0);radius=Q(0)
        for coeff,error in models[source]:
            center=op.exact_integral(coeff,Q(1,64));r=error/Q(64)
            encoded={'lower':op.rational(center-r),'upper':op.rational(center+r)}
            stable.update(json.dumps(encoded,sort_keys=True,separators=(',',':')).encode())
            total+=center;radius=op.outward_radius(radius+r)
            local+=1
        stable.update(json.dumps({'lower':op.rational(total-radius),'upper':op.rational(total+radius)},
                                 sort_keys=True,separators=(',',':')).encode())
        max_radius=max(max_radius,radius);whole+=1
    return {'local_rows':local,'whole_rows':whole,'wall_seconds':time.monotonic()-start,
            'stable_serialized_work_sha256':stable.hexdigest(),
            'maximum_whole_work_radius':op.rational(max_radius),
            'scope':'fabricated own geometry-work polynomial; not registered Lg'}


def main():
    if len(sys.argv)!=2: raise SystemExit('Usage: fabricated_complete_route.py FRESH_RECEIPT.json')
    target=Path(sys.argv[1])
    if target.exists(): raise SystemExit('Fresh receipt required')
    started=time.monotonic()
    formal=formal_checks();checks,mutations=op.fabricated_tests()
    models,work_models,model_receipt=fabricated_models()
    operator_receipt=op.fabricated_benchmark(models)
    work_receipt=fabricated_work_benchmark(work_models)
    total=time.monotonic()-started
    receipt={'status':'PASS_FABRICATED_COMPLETE_INDEPENDENT_OPERATOR_ONLY',
             'python_optimization':sys.flags.optimize,'formal_model_checks':formal,
             'operator_checks':checks,'mutation_controls':mutations,
             'model_benchmark':model_receipt,'operator_benchmark':operator_receipt,
             'work_benchmark':work_receipt,
             'total_wall_seconds':total,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
             'physical_source_evaluations':0,'retained_arrays_decoded':0,'primary_helper_imports':0,
             'registered_source_builder_success_calls':0,
             'direct_pressure_full_twelve_case_certificate':'UNRESOLVED',
             'source_pins':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in (Path(__file__),Path(op.__file__),Path(sm.__file__))}}
    target.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':receipt['status'],'total_wall_seconds':total,
                     'model_wall_seconds':model_receipt['wall_seconds'],
                     'operator_wall_seconds':operator_receipt['wall_seconds'],
                     'peak_rss_kib':receipt['peak_rss_kib'],'local_rows':operator_receipt['cases'],
                     'whole_rows':operator_receipt['whole_interval_cases'],
                     'work_panel_rows':work_receipt['local_rows'],'work_whole_rows':work_receipt['whole_rows']}))


if __name__=='__main__':main()
