#!/usr/bin/env python3
"""Fabricated arrays only; has no paths to saved physical arrays or outcomes."""
from decimal import Decimal
from fractions import Fraction
import argparse
import hashlib
import json
from pathlib import Path
import resource
import time
import mpmath
import numpy as np
import ledger_primary as lp


def fabricate(source, setting):
    LD=np.longdouble; CD=np.clongdouble
    n=8192 if setting=='coarse' else 16384
    count=577 if setting=='coarse' else 1153
    dt=LD(1)/(128 if setting=='coarse' else 256)
    k=(np.arange(n,dtype=LD)+LD('.5'))*LD(256)/n
    weight=np.full(n,LD(256)/n,dtype=LD)
    pi=np.arccos(LD(-1)); eps=LD('.0001')
    amplitude=LD('0.00000001')*(1 if source=='positive_B' else -1)
    d=(amplitude*(np.cos(k)+CD(1j)*np.sin(3*k)))/(1+k*k)
    # Unrestricted real c deliberately survives; this is not a normalized projection.
    c=(LD('0.000000000000001')*np.cos(2*k)+CD(1j)*LD('0.0000000000001'))/(1+k)
    u=c+d; w=CD(2j)*k*d
    E=np.exp(CD(2j)*k)
    ub=u+np.expm1(CD(2j)*k)*d; wb=E*w
    times=LD(-6)+np.arange(count,dtype=LD)*dt
    values=np.zeros((count,3,9),dtype=LD)
    ledger=np.zeros((count,3),dtype=LD)
    mu=weight*k*k/(2*pi*pi)/(2*k*eps)
    prefixes=[int(np.count_nonzero(k<cutoff)) for cutoff in lp.CUTOFFS]
    for j,t in enumerate(times):
        L=-LD(1)/t
        phase=np.exp(CD(2j)*k*(t+LD('2.5')))
        R=mu*((2*k*k+3*L*L)*c.real+((3*L*L*d-L*w)*phase).real)
        P=mu*((2*k*k/3-L*L)*c.real+((-(4*k*k/3+L*L)*d-L*w)*phase).real)
        F=2*mu*(3*L**3*c.real+((2*k*k*L+3*L**3)*d*phase+L*L*w*phase).real)
        for b,end in enumerate(prefixes):
            values[j,b,3]=np.sum(R[:end],dtype=LD)
            values[j,b,4]=np.sum(P[:end],dtype=LD)
            ledger[j,b]=np.sum(F[:end],dtype=LD)
    geometry=np.zeros((count,10),dtype=LD)
    geometry[:,0]=geometry[:,1]=-LD(1)/times
    return {'k':k,'momentum_weights':weight,'u_4':u.astype(CD),'w_4':w.astype(CD),
            'u_5':ub.astype(CD),'w_5':wb.astype(CD),'history_eta':times,
            'observation_eta':np.array((-5.5,-4.5,-4,-3.5,-2.5,-1.5),dtype=LD),
            'history_values':values,'history_ledger_integrand':ledger,
            'history_baseline_contact':np.zeros((count,3),dtype=LD),
            'history_source_jet':np.zeros((count,6),dtype=LD),
            'history_forcing_jet':np.zeros((count,4),dtype=LD),'history_geometry':geometry,
            'history_cutoffs':np.array(lp.CUTOFFS,dtype=np.int64),
            'history_quantity_names':np.array(lp.QUANTITIES),
            'history_geometry_names':np.array(lp.GEOMETRY)}


def pure_checks():
    checks=[]
    def check(name,ok):
        lp.require(ok,name)
        checks.append(name)
    for p in (1,3,266,333):
        for n in (-13,-12,-11,-5,-3,-1,0,1,3,5,11,12,13):
            check(f'ties_even_p{p}_n{n}',lp.round_shift(n,p)==round(Fraction(n,1<<p)))
    for dps,p in lp.LEVELS:
        ctx=mpmath.mp.clone();ctx.dps=dps
        for raw in ('0','1','-1','0.0001','1.000000000000000000867361737988403547205962240695953369140625'):
            value=np.longdouble(raw)
            x=lp.represented(ctx,value)
            check(f'LD_ratio_{dps}_{raw}',lp.fraction_of_mpf(x)==Fraction(*value.as_integer_ratio()))
        for name,value in (('negative_zero',np.longdouble('-0')),('tiny',np.finfo(np.longdouble).smallest_subnormal),
                           ('largest',np.finfo(np.longdouble).max)):
            x=lp.represented(ctx,value)
            check(f'LD_extreme_{dps}_{name}',lp.fraction_of_mpf(x)==Fraction(*value.as_integer_ratio()))
        # Exact phase recurrence vs separately evaluated MP exponential, with unprojected complex seed.
        for kval in ('0.00390625','64','200','255.99609375'):
            k=ctx.mpf(kval);omega=2*k;step=ctx.mpf(1)/256
            alpha=ctx.mpc('0.001234567890123456789','-0.000987654321987654321')
            rot=ctx.exp(ctx.j*omega*step)
            er,ei=lp.quantize(rot.real,p),lp.quantize(rot.imag,p)
            zr,zi=lp.quantize(alpha.real,p),lp.quantize(alpha.imag,p)
            maximum=ctx.zero
            for j in range(257):
                direct=alpha*ctx.exp(ctx.j*omega*j*step)
                maximum=max(maximum,abs(ctx.mpc(lp.dyadic(ctx,zr,p),lp.dyadic(ctx,zi,p))-direct))
                zr,zi=lp.round_shift(zr*er-zi*ei,p),lp.round_shift(zr*ei+zi*er,p)
            check(f'direct_MP_phase_{dps}_{kval}',maximum<ctx.mpf('1e-70'))
    x=np.longdouble(1)+np.longdouble(2)**-60
    check('binary64_downcast_detection',float(x)==1 and x!=1)
    try: lp.represented(mpmath.mp.clone(),float(x))
    except ValueError: check('reject_prior_binary64_downcast',True)
    else: raise ValueError('Prior binary64 loss accepted')
    # Native prefix/reset orders differ under cancellation: do not replace one by the other.
    h=np.array([np.longdouble('1e30'),0,0,0,0,4,5,6,7],dtype=np.longdouble)
    global_reset=lp.simpson_to(h,8,np.longdouble('.125'))-lp.simpson_to(h,4,np.longdouble('.125'))
    direct_reset=lp.simpson_to(h[4:],4,np.longdouble('.125'))
    check('global_prefix_contract_is_distinct',global_reset!=direct_reset)
    mutations=0
    for bad in ('missing_mode','wrong_dtype','nonzero_source','wrong_profile_time','wrong_quantity_name','wrong_cutoff'):
        a=fabricate('positive_B','coarse')
        if bad=='missing_mode': del a['u_5']
        elif bad=='wrong_dtype': a['k']=a['k'].astype(np.float64)
        elif bad=='nonzero_source': a['history_source_jet'][-1,0]=1
        elif bad=='wrong_profile_time': a['history_eta'][-1]=np.longdouble('-1.49')
        elif bad=='wrong_quantity_name': a['history_quantity_names'][3]='q'
        elif bad=='wrong_cutoff': a['history_cutoffs'][0]=65
        try: lp.validate_arrays(a,'coarse')
        except ValueError: mutations+=1
        else: raise ValueError('Mutation accepted '+bad)
    def minimal_rows():
        return [{'K':k,'R_profile_max_error':'0','P_profile_max_error':'0','D_cont':'0',
                 'E_flow':'0','decomposition_error':'0','recurrence_primitive_gap':'0'} for k in lp.CUTOFFS]
    r80=minimal_rows();r100=minimal_rows()
    r80[0]['R_profile_max_error']=r100[0]['R_profile_max_error']='0.0000002000000000000000000001'
    _,negative,fatal=lp.gap_and_gates(r80,r100)
    check('physical_consistency_failure_is_negative_not_fatal',bool(negative) and not fatal)
    r80=minimal_rows();r100=minimal_rows()
    r80[0]['R_profile_max_error']=r100[0]['R_profile_max_error']=lp.GATES['profile']
    _,negative,fatal=lp.gap_and_gates(r80,r100)
    check('exact_profile_gate_boundary_is_accepted',not negative and not fatal)
    r80=minimal_rows();r100=minimal_rows()
    r80[0]['recurrence_primitive_gap']=r100[0]['recurrence_primitive_gap']='0.0000000000011'
    _,negative,fatal=lp.gap_and_gates(r80,r100)
    check('shared_primitive_bias_fails_each_precision_despite_zero_gap',not negative and len(fatal)==2)
    r80=minimal_rows();r100=minimal_rows()
    r80[0]['decomposition_error']=r100[0]['decomposition_error']='0.000000000002'
    _,negative,fatal=lp.gap_and_gates(r80,r100)
    check('shared_arithmetic_closure_failure_is_fatal',not negative and len(fatal)==2)
    return {'checks':checks,'checks_count':len(checks),'rejected_mutations':mutations,
            'physical_input_files_opened':0,'physical_source_evaluations':0}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--full-shape',action='store_true')
    args=parser.parse_args()
    lp.require(not args.output.exists(),'Fresh receipt required')
    started=time.monotonic()
    source_before=lp.sha(Path(lp.__file__))
    test_before=lp.sha(Path(__file__))
    checks=pure_checks()
    result={'status':'PASS_SYNTHETIC_PURE_PREFLIGHT','pure_checks':checks}
    if args.full_shape:
        main_started=time.monotonic()
        run=lp.run_arrays(fabricate,main_started)
        lp.require(not run['failures'] and not run['fatal_failures'],'Full fabricated diagnostic failed')
        result['full_shape_run']=run
        result['status']='PASS_FULL_SHAPE_SYNTHETIC_PREFLIGHT'
    result['total_preflight_seconds']=time.monotonic()-started
    result['peak_rss_kib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    lp.require(lp.sha(Path(lp.__file__))==source_before and lp.sha(Path(__file__))==test_before,
               'Synthetic source changed during test; receipt cannot bind a completed source')
    result['source_sha256']=source_before
    result['test_sha256']=test_before
    result['physical_saved_data_loaded']=False
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'seconds':result['total_preflight_seconds'],
                      'peak_rss_kib':result['peak_rss_kib'],'output':str(args.output)}))


if __name__=='__main__':main()
