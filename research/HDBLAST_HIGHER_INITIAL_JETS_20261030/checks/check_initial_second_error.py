#!/usr/bin/env python3
"""Continuous higher-initial-error caps and an affine-extension endpoint defect.

The fixed smoothed nominal data are unchanged.  Every error comparison is an
exact rational interval or explicitly bounded convolution remainder.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode=True

CORE_PIN='f2d1bfef82220fc1e385a3bdb6de35487b1ef28c4caee37c8dddb6a2b68cb782'


def need(ok,message):
    if not ok:raise RuntimeError(message)


def pin(path,digest):
    data=path.read_bytes()
    need(digest is None or hashlib.sha256(data).hexdigest()==digest,'INPUT_HASH:'+str(path))
    return data


def encode(z):
    if isinstance(z,F):return str(z)
    if isinstance(z,dict):return {k:encode(v) for k,v in z.items()}
    if isinstance(z,(list,tuple)):return [encode(v) for v in z]
    return z


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',required=True,type=Path)
    ap.add_argument('--cells',default=64,type=int)
    ap.add_argument('--endpoint-only',action='store_true')
    args=ap.parse_args()
    need(64<=args.cells<=1024 and args.cells%64==0,'mesh divisible by64')
    core=Path(__file__).with_name('initial_second_core.py')
    raw=pin(core,CORE_PIN)
    spec=importlib.util.spec_from_file_location('higher_initial_local',core)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    h=mod.HigherInitial(args.repo_root);c=h.ctx;I=h.I;ell=h.ell
    base=h.root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial'
    report=json.loads(pin(base/'results/INITIAL_REPORT.json',mod.INITIAL_PIN))
    profiles=json.loads(pin(base/'results/INITIAL_PROFILES.json',mod.PROFILES_PIN))
    labels=list(c.op.LABELS);weights=c.weights
    need(report['labels']==profiles['labels']==labels and
         list(map(F,profiles['gram']))==weights,'same full23 initial data')
    du,dv,R=F(1,2**30),F(1,2048),F(7,4)
    C1,C2=F(65,64),F(66,64)
    need(2*du<min(C2-C1,R-C2) and du<1,'separated curvature jumps')
    need(report['status']=='PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER', 'nominal derivative certificate')
    M2u,M3u,M2v=[list(map(F,row)) for row in report['nominal_second_derivative_component_caps']]
    mu1u,mu1v=35*du/128,35*dv/128

    # One-sided nominal curvature values. These are piecewise/a.e. second
    # derivatives, since the nominal scalar displacement is only C1 at joins.
    jumps=[]
    for x,left,right in ((C1,'zero','cutoff'),(C2,'cutoff','plain'),(R,'plain','affine')):
        lu,_=h.nominal(I(x),left)
        L=[2*ell**2*z.c[2] for z in lu]
        if right=='affine':U=[I(0)]*23
        else:
            ru,_=h.nominal(I(x),right);U=[2*ell**2*z.c[2] for z in ru]
        D=[b-a for a,b in zip(L,U)]
        jumps.append({'r':x,'left_branch':left,'right_branch':right,
                      'left_second_derivative':[[z.lo,z.hi] for z in L],
                      'right_second_derivative':[[z.lo,z.hi] for z in U],
                      'jump_abs_component_caps':[z.mag() for z in D]})
    maxjump=[max(row['jump_abs_component_caps'][i] for row in jumps) for i in range(23)]
    smoothing_u2=[M3u[i]*mu1u+maxjump[i]/2 for i in range(23)]

    # Exact same P1 velocity nodes; their endpoint radii are bounded by each
    # cell's authenticated interpolation cap. They are not assumed exact
    # samples of the analytic nominal velocity.
    nodes=list(map(F,profiles['nodes']));vv=[list(map(F,row)) for row in profiles['profiles']['v']]
    derivative_P1=[F(0)]*23
    need(len(report['interpolation_cells'])==len(nodes)-1,'P1 partition')
    for k,cell in enumerate(report['interpolation_cells']):
        left,right=nodes[k:k+2];step=right-left
        need(list(map(F,cell['r']))==[left,right] and step>0,'P1 exact cell cover')
        M2=list(map(F,cell['second_derivative_component_caps'][2]))
        radius=list(map(F,cell['component_error_caps'][2]))
        for i in range(23):
            # The cell cap is >= both endpoint radii. The derivative deviation
            # from its exact secant is <= step*M2/2 by integral averaging.
            derivative_P1[i]=max(derivative_P1[i],step*M2[i]/2+2*radius[i]/step)
    # For this proof only, extend the analytic nominal velocity by its own
    # tangent line beyond R (and zero through its flat inner collar). Its
    # derivative is globally Lipschitz with piecewise cap M2v. This auxiliary
    # function is not substituted for the candidate's P1 velocity. The latter
    # extends with the final P1 slope; the derivative discrepancy from the
    # auxiliary tangent beyond R is bounded by the same final-cell secant
    # and nodal-radius charge above. Positivity of the kernel therefore gives
    # |v_nom'(r)-(k_dv*P1v)'(r)| <= Dv + M2v E|s|, including near R.
    smoothing_v1=[M2v[i]*mu1v for i in range(23)]
    need(nodes[-1]==R and nodes[1]==C1 and all(z==0 for row in vv[:2] for z in row),
         'registered zero first cell and affine final-slope extension')

    # Endpoint proof exploits the even kernel and affine tangent extension.
    # The outer half has zero curvature; the inner half contributes N2(R)/2
    # with remainder bounded by half of its first absolute moment.
    pu,pv=h.paired(I(R));nu,nv=h.nominal(I(R),'plain')
    endpoint=[]
    lower_square=F(0);nonzero=[]
    for i,label in enumerate(labels):
        N2=2*ell**2*nu[i].c[2]
        D2=sum((2*ell**2*pu[i].d[k].c[2] for k in range(3)),I(0))
        actual=N2+D2
        limit=N2/2+D2
        remainder=M3u[i]*35*du/256
        fixed=limit+I(-remainder,remainder)
        smooth=N2/2+I(-remainder,remainder)
        minimum=F(0) if fixed.lo<=0<=fixed.hi else min(abs(fixed.lo),abs(fixed.hi))
        lower_square+=weights[i]*minimum**2
        if minimum:nonzero.append(label)
        endpoint.append({'label':label,'nominal_left_u_rr':[N2.lo,N2.hi],
            'paired_actual_minus_nominal_u_rr':[D2.lo,D2.hi],
            'actual_u_rr':[actual.lo,actual.hi],
            'smoothed_nominal_u_rr_at_fixed_width':[smooth.lo,smooth.hi],
            'actual_minus_smoothed_u_rr_at_fixed_width':[fixed.lo,fixed.hi],
            'actual_minus_half_nominal_curvature_limit':[limit.lo,limit.hi],
            'fixed_width_half_kernel_remainder_abs_upper':remainder})
    need(endpoint[12]['label']=='w:0','registered scalar mean')
    mean=endpoint[12]['actual_minus_smoothed_u_rr_at_fixed_width']
    need(mean[1]<-F('1.6245e33') and lower_square>F('1.6245e33')**2,
         'strict actual endpoint curvature-error lower bound')
    need(-F('1.624543e33')<mean[0]<=mean[1]<-F('1.624538e33'),
         'readable signed endpoint curvature-error bracket')
    need(endpoint[12]['actual_minus_half_nominal_curvature_limit'][1]<-F('1.6245e33'),
         'nonzero kernel-width limit')

    continuous={}
    if not args.endpoint_only:
        caps=[[[F(0)]*23 for _ in range(2)] for _ in range(3)]
        cells=[]
        for k in range(args.cells,7*args.cells//4):
            left,right=F(k,args.cells),F(k+1,args.cells)
            u,v=h.paired(I(left,right))
            errors=[[[ (2*ell**2*z.d[ch].c[2]).mag() for z in u],
                     [ (ell*z.d[ch].c[1]).mag() for z in v]] for ch in range(3)]
            for ch in range(3):
                for jet in range(2):
                    caps[ch][jet]=[max(a,b) for a,b in zip(caps[ch][jet],errors[ch][jet])]
            cells.append({'r':[left,right],'channel_H_caps':[[c.norm(row) for row in channel] for channel in errors]})
        need(cells[0]['r'][0]==1 and cells[-1]['r'][1]==R and
             all(a['r'][1]==b['r'][0] for a,b in zip(cells,cells[1:])), 'whole active continuous cover')
        Hchannels=[[c.norm(row) for row in channel] for channel in caps]
        total_u2=sum((row[0] for row in Hchannels),F(0))+c.norm(smoothing_u2)
        total_v1=sum((row[1] for row in Hchannels),F(0))+c.norm(derivative_P1)+c.norm(smoothing_v1)
        need(total_u2**2>lower_square,'continuous upper and endpoint lower are consistent')
        need(total_u2<F('8.064e33') and total_v1<F('9.672e24'),
             'strict readable continuous higher-initial-error caps')
        continuous={'partition':cells,'continuous_cells':len(cells),
            'paired_channel_component_caps':caps,'paired_channel_H_caps':Hchannels,
            'actual_minus_smoothed_u_rr_H_upper':total_u2,
            'actual_minus_smoothed_P1_v_r_H_upper':total_v1,
            'readable_strict_u_rr_H_upper':'8.064e33',
            'readable_strict_P1_v_r_H_upper':'9.672e24'}
    h.authenticate();pin(core,hashlib.sha256(raw).hexdigest())
    result={'status':'PASS','method':'order5 paired initial jets plus jump-aware displacement convolution and rounded-node P1 velocity-gradient error',
        'labels':labels,'gram':weights,'ell':ell,'footprint':['0',R],
        'fixed_displacement_kernel_width':du,'fixed_velocity_kernel_width':dv,
        'initial_data_unchanged':True,'candidate':'same fixed W1 as PR97-99',
        'piecewise_nominal_M2u':M2u,'piecewise_nominal_M3u':M3u,'piecewise_nominal_M2v':M2v,
        'curvature_jumps':jumps,'displacement_curvature_smoothing_component_caps':smoothing_u2,
        'P1_velocity_derivative_component_error_caps':derivative_P1,
        'velocity_gradient_convolution_component_caps':smoothing_v1,
        'endpoint':{'r':R,'rows':endpoint,'fixed_width_nonzero_component_labels':nonzero,
            'fixed_width_error_H_norm_squared_lower':lower_square,
            'strict_fixed_width_scalar_mean_error_abs_lower':'1.6245e33',
            'strict_fixed_width_scalar_mean_signed_bracket':['-1.624543e33','-1.624538e33'],
            'strict_width_limit_scalar_mean_error_abs_lower':'1.6245e33',
            'width_limit':'actual_u_rr(R)-nominal_u_rr(R-)/2 for even normalized kernels and fixed affine tangent extension'},
        'continuous':continuous,'complete_closures_authenticated_before_and_after':True,
        'input_sha256':{'initial_producer':mod.PRODUCER_PIN,'initial_report':mod.INITIAL_PIN,'initial_profiles':mod.PROFILES_PIN,'local_order5_core':hashlib.sha256(raw).hexdigest()},
        'scope':'higher derivatives of initial DATA ERROR for fixed auxiliary numerical candidate; affine extension gives genuine nonzero endpoint curvature mismatch, not physical-model obstruction or nonlinear-spacetime result; no complete C2 error gate or full Einstein residual certified'}
    print(json.dumps(encode(result),sort_keys=True))


if __name__=='__main__':main()
