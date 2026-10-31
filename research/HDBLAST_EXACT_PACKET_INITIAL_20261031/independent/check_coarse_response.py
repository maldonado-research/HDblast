#!/usr/bin/env python3
"""Fresh Fraction response ledger using only independent displacement caps.

No new candidate producer is imported or read. Old velocity is recomputed
absolutely from the authenticated P1 table; its complete error is retained.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True

FILES = {
 'c2':('HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference/C2.json','85bf4133d931d93401b339dc508553ae488bf91f43f75bd7e736432cc856a923'),
 'combined':('HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference/COMBINED.json','3e4fb0174df7e88c566346a1b68aedbe4e731fbd062de11eafa28a2cae0a4989'),
 'absolute':('HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/CANDIDATE.json','44d05eadf32770e152ba15e8482f2d290a8baafdcf59f20ac281841960b76c1f'),
 'response':('HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030/reference/RESPONSE.json','ef097736c337003f0e2d578a830d14824772065cc36f950932819f524c6fbe84'),
 'initial':('HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/INITIAL.json','bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207'),
 'transfer':('HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/TRANSFER.json','4d0b83994270f6e2d3aa5433887a0d4b978acbac372bf28197050345de2881c3'),
 'time':('HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/TIME.json','6bf1b2fbb2b6f6a9df2f2a24bcac142bfe881437fa5638bcff20b20c7b3b509d'),
 'fx':('HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json','63ed8e3a7c0544db93de0d3bcbdba4652eb4d9a6bcd8f434327e785cd9986e60'),
 'profiles':('HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/results/INITIAL_PROFILES.json','555f911f6c336862dfacb6aa56ef5179ac06701812d593b5bfbbf6ea9bd4d5d1'),
}
INDEPENDENT_PIN='fc9ef547bd4e52f846b59a153b8552b483c7f53250d4d818e3d1ad5baffcda22'


def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(x):
    if isinstance(x,Q):
        return str(x)
    if isinstance(x,dict):
        return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [encode(v) for v in x]
    return x


def sqrtup(q):
    need(q>=0,'NEGATIVE_SQUARE')
    if q==0:
        return Q(0)
    grid=1<<128
    return Q(math.isqrt((q*grid*grid).__ceil__())+1,grid)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source-root','--root',dest='root',type=Path,required=True)
    ap.add_argument('--initial-report','--independent',dest='independent',type=Path,
                    default=Path(__file__).resolve().parent.parent/'reference/INDEPENDENT_INITIAL.json')
    ap.add_argument('--independent-sha256',default=INDEPENDENT_PIN)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    need(not args.output.exists(),'OUTPUT_MUST_BE_ABSENT')
    need(args.independent_sha256==INDEPENDENT_PIN,'FROZEN_PORTABLE_INDEPENDENT_REPORT_PIN')
    self_sha=sha(Path(__file__))
    docs={};pins={}
    for key,(rel,digest) in FILES.items():
        p=args.root/'research'/rel
        need(not p.is_symlink() and p.is_file(),'INPUT_TYPE:'+key)
        need(sha(p)==digest,'INPUT_HASH:'+key)
        docs[key]=json.loads(p.read_text());pins[key]=digest
    need(sha(args.independent)==args.independent_sha256,'INDEPENDENT_INPUT_HASH')
    own=json.loads(args.independent.read_text())
    need(own['full_coverage_certificate_in_this_report'] is True and
         own['shared_packet_and_E_approximation_difference_exactly_zero'] is True,
         'INDEPENDENT_FULL_DISPLACEMENT_CERTIFICATE_ROLE')
    cover=own['independent_full_displacement_coverage']
    need(cover['initial_gate_proved_by_this_mesh']==[True,True,True],'INITIAL_GATES_ALREADY_PROVED')
    labels=own['label_order'];weights=list(map(Q,own['Gram']))
    c2,transfer,response,time,ab,initial,profile,fx=[docs[k] for k in
        ('c2','transfer','response','time','absolute','initial','profiles','fx')]
    need(labels==c2['labels']==transfer['labels']==profile['labels'] and
         own['Gram']==c2['Gram']==transfer['Gram']==profile['gram'] and
         len(labels)==23 and min(weights)>0,'COMPLETE_SAME_POSITIVE_FULL23_NORM')
    norm=lambda caps:sqrtup(sum((w*z*z for w,z in zip(weights,caps)),Q(0)))
    act={k:Q(v) for k,v in response['actual_coefficient_H_norm_caps'].items()}
    frozen={k:Q(v) for k,v in transfer['frozen_H_operator_caps'].items()}
    diff={k:Q(v) for k,v in response['actual_coefficient_error_H_norm_caps'].items()}
    der={k:Q(v) for k,v in time['global']['H_operator_caps'].items()}
    need(all(act[k]==frozen[k]+diff[k] for k in act),'ACTUAL_OPERATOR_SPLIT')
    need(der['B_x']==der['M_x']==0 and c2['boundary_error']=={'g':'0','g_x':'0','M_x':'0'},
         'FIXED_ACTUAL_BOUNDARY_AND_B')
    need(fx['scope']['actual_original_full23_F_x_bound'] is True and
         fx['source_coordinate']=='x=t/ell','ORIGINAL_FULL_UNFILTERED_FX_ROLE')
    T=Q(3,2);bd=act['M']/4;P=1/(1-bd);q=Q(response['retained_coupled_X_operator_cap'])
    bulk=act['A']+act['B']+Q(3,4)*act['D']
    bulk0=frozen['A']+frozen['B']+Q(3,4)*frozen['D']
    diffbulk=diff['A']+diff['B']+Q(3,4)*diff['D']
    comm=der['A_x']+der['B_x']+Q(3,4)*der['D_x']
    need(0<bd<1 and T*bulk+bd<=q==Q(3409,8000000)<1,'BOTH_CAUSAL_INVERSES')
    need(P==Q(transfer['P_boundary_inverse_X_cap']),'SAME_ACTUAL_ROBIN_INVERSE')
    FH=Q(ab['source_caps']['F_H']);Fx=Q(fx['genuine_scaled_F_x_H_cap']);Gx=Q(ab['source_caps']['G_x_H']);G=Q(487122)
    need(FH==17387513156573 and Fx==905405545954284 and Gx==58807538,'UNCHANGED_ORIGINAL_SOURCE_CAPS')
    KF=Q(52162547172519,2)+G
    olderr={k:Q(v) for k,v in c2['initial_error_H_caps'].items()}
    oldabs={k:Q(v) for k,v in ab['initial_absolute_H_caps'].items()}
    need(olderr['v0_r']==Q(initial['continuous']['actual_minus_smoothed_P1_v_r_H_upper']),
         'UNCHANGED_COMPLETE_VELOCITY_GRADIENT_ERROR')
    oldBs=Q(ab['free_wave_X_cap']);oldD=Q(transfer['direct_reproof']['residual_majorants']['W_minus_B0_whole_cone_X_cap'])
    oldC=Q(transfer['direct_reproof']['new_actual_bulk_correction_X_cap']);oldE=Q(c2['complete_original_C1_error_X_cap'])

    def differentiated(absolute,errors,Bs,D,C,E):
        qfree=absolute['u_rr']+absolute['v_r']
        f00=frozen['A']*absolute['v']+frozen['B']*absolute['u_r']+frozen['D']*absolute['u']
        f0x=(frozen['A']+frozen['B'])*qfree+frozen['D']*Bs
        Dx=P*(f00+T*f0x+act['M']*Bs+Gx)
        r00=diff['A']*absolute['v']+diff['B']*absolute['u_r']+diff['D']*absolute['u']+FH
        r0x=(act['A']+act['B'])*Dx+(act['D']+comm)*D+(diff['A']+diff['B'])*qfree+(diff['D']+comm)*Bs+Fx
        Cx=P*(r00+T*r0x)
        F1=(act['A']+act['B'])*Cx+(act['D']+comm)*C
        w=errors['e0_rr']+act['A']*errors['v0']+act['B']*errors['e0_r']+act['D']*errors['e0']
        J=max((errors['v0']+T*w)/3,errors['v0_r']+w)
        Z=(J+T*(F1+comm*E))/(1-q)
        return {'free_second':qfree,'f00':f00,'f0x':f0x,'D0xX':Dx,'r00':r00,'r0x':r0x,'CxX':Cx,'F1':F1,'w0':w,'J':J,'xx':Z,'xr':Z,'rr':Z+bulk*E+bulk*C}
    replay=differentiated(oldabs,olderr,oldBs,oldD,oldC,oldE)
    need(replay['F1']==Q(c2['original_reduced_defect_time_derivative_H_cap']),
         'REPRODUCE_RETAINED_COMPLETE_DIFFERENTIATED_LEDGER')
    need(all(replay[k]==Q(c2['original_target_second_coordinate_error_H_caps'][k])
         for k in ('xx','xr','rr')),'REPRODUCE_ALL_THREE_OLD_EXACT_C2_CAPS')
    # Exact complete P1 velocity, final-secant continuation, and positive
    # even convolution. Componentwise extrema enclose its full footprint.
    nodes=list(map(Q,profile['nodes']));vv=[list(map(Q,row)) for row in profile['profiles']['v']]
    delta=Q(1,2048)
    need(nodes[0]==0 and nodes[1]==Q(65,64) and nodes[-1]==Q(7,4)
         and nodes[1]-delta>1 and all(v==0 for row in vv[:2] for v in row),
         'UNCHANGED_VELOCITY_FLAT_COLLAR_AND_FOOTPRINT')
    slopes=[[(y-x)/(right-left) for x,y in zip(vv[j],vv[j+1])]
            for j,(left,right) in enumerate(zip(nodes,nodes[1:]))]
    vcaps=[max(max(abs(row[i]) for row in vv),abs(vv[-1][i]+delta*slopes[-1][i])) for i in range(23)]
    vrcaps=[max(abs(row[i]) for row in slopes) for i in range(23)]
    absolute={k:Q(v) for k,v in zip(('u','u_r','u_rr'),cover['candidate_absolute_H_caps'])}
    absolute.update(v=norm(vcaps),v_r=norm(vrcaps))
    errors={k:Q(v) for k,v in zip(('e0','e0_r','e0_rr'),cover['error_H_caps'])}
    errors.update(v0=olderr['v0'],v0_r=olderr['v0_r'])
    Bs=max((absolute['u']+T*absolute['v'])/3,absolute['u_r']+absolute['v'])
    I=max((errors['e0']+T*errors['v0'])/3,errors['e0_r']+errors['v0'])
    D=P*((T*bulk0+bd)*Bs+G)
    R=T*(diffbulk*Bs+bulk*D)+KF
    C=P*R
    E=(I+(q-bd)/(1-bd)*R)/(1-q)
    result=differentiated(absolute,errors,Bs,D,C,E)
    gate={'C1':E<=oldE}
    gate.update({k:result[k]<=Q(c2['original_target_second_coordinate_error_H_caps'][k])/2
                 for k in ('xx','xr','rr')})
    need(all(gate.values()),'FULL_RESPONSE_GATES_WITH_OWN_INDEPENDENT_DISPLACEMENT_CAPS')
    # Meaningful arithmetic controls: these have explicit nonzero exact
    # witnesses and cannot silently erase retained initial/source terms.
    controls={'omit_f00_initial_acceleration':P*result['f00'],
              'omit_r00_initial_acceleration':P*result['r00'],
              'omit_original_Fx':T*T*(act['A']+act['B'])*P*Fx/(1-q),
              'omit_velocity_gradient_error':errors['v0_r'],
              'absolute_error_u_rr_interchange':absolute['u_rr']-errors['e0_rr']}
    need(all(v>0 for v in controls.values()),'EXPLICIT_NONZERO_OMISSION_WITNESSES')
    for key,(rel,digest) in FILES.items():
        need(sha(args.root/'research'/rel)==digest,'INPUT_CHANGED:'+key)
    need(sha(args.independent)==args.independent_sha256 and sha(Path(__file__))==self_sha,
         'INDEPENDENT_INPUT_AND_AUDITOR_UNCHANGED')
    report={'status':'PASS_COARSE_INDEPENDENT_COMPLETE_RESPONSE_LEDGER',
            'candidate':'W1_exact_packet','no_new_candidate_producer_imports_or_numbers':True,
            'complete_absolute_initial_caps':absolute,'complete_initial_error_caps':errors,
            'independent_displacement_input_sha256':args.independent_sha256,
            'source_pins_before_and_after':pins,'auditor_sha256_before_and_after':self_sha,
            'velocity_component_absolute_caps':vcaps,'velocity_component_radial_caps':vrcaps,
            'fresh_global_ledger':{'Bs':Bs,'I':I,'D0X':D,'R':R,'CX':C,'E':E,'F0':bulk*C},
            'complete_differentiated_ledger':result,'full_response_gate':gate,
            'old_ledger_reproduced_exactly':True,'old_C1_threshold':oldE,
            'old_C2_half_thresholds':{k:Q(c2['original_target_second_coordinate_error_H_caps'][k])/2 for k in ('xx','xr','rr')},
            'nonzero_omission_witnesses':controls,
            'analytic_premises':'Same actual F,Fx,G,M and operators; componentwise C2 displacement, C1 velocity; zero wedge; common first Cauchy jet; physical-side regularity and compatible corner. r1(0)=0, correction initial accelerations retained.',
            'scope':'Conditional complete physical-side coordinate error upper bounds. Coarse independent displacement caps, independently recomputed unchanged velocity absolute caps, and retained complete unchanged velocity errors. No pointwise physical ratio or global true-error reduction claimed.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(encode(report),indent=2,sort_keys=True)+'\n')
    print(report['status'])


if __name__=='__main__':
    main()
