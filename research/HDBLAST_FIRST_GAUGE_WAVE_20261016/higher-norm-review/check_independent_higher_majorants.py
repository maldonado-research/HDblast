#!/usr/bin/env python3
"""Independent exact majorants via finite Neumann sums plus a scalar tail.

No producer module is imported and no matrix inverse is computed. This checks
all twelve displayed order caps from the authenticated coefficient report.
The accompanying review supplies the differential and regularity reasoning.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

COEFFICIENT_SHA='4009c46f6d63583b4d161eebbb620077467b7561a3c4b5728d5293655d340a53'
CANDIDATE_SHA='f246bc6840ea36f123c7c53da28a9c04925954187e0f06df25add779ce083c6d'
T=Q(3)
THETA=T/2
STEPS=16


def need(ok,message):
    if not ok:
        raise ValueError(message)


def dot(A,v):
    return [sum((a*b for a,b in zip(row,v)),Q(0)) for row in A]


def add(*vs):
    return [sum(col,Q(0)) for col in zip(*vs)]


def mul(c,v):
    return [c*x for x in v]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--coefficients',type=Path,required=True)
    parser.add_argument('--candidate',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--inject',choices=['initial_velocity','commutator','radial_recovery'])
    args=parser.parse_args()
    need(not args.output.exists(),'output already exists')
    raw=args.coefficients.read_bytes()
    need(hashlib.sha256(raw).hexdigest()==COEFFICIENT_SHA,'coefficient authentication failed')
    c=json.loads(raw)
    candidate_raw=args.candidate.read_bytes()
    need(hashlib.sha256(candidate_raw).hexdigest()==CANDIDATE_SHA,'candidate authentication failed')
    candidate=json.loads(candidate_raw)
    need(candidate['status']=='PASS_CONTINUOUS_RAW_GAUGE_ALL_MIXED_C3_CAPS','candidate status')
    def mat(name,i=0,j=0):
        return [[Q(entry[f'{i},{j}']['abs_upper']) for entry in row]
                for row in c['matrix_mixed_derivative_intervals'][name]]
    def combined(i=0,j=0):
        A,B,D=[mat(k,i,j) for k in ('A','B','D')]
        return [[A[k][l]+B[k][l]+THETA*D[k][l] for l in range(3)] for k in range(3)]
    A,B,D=[mat(name) for name in ('A','B','D')]
    L=combined()
    M=[[T*x for x in row] for row in L]
    q=max(sum(row,Q(0)) for row in M)
    need(0<q<1,'Neumann convergence')
    tails=[]
    def resolvent_upper(rhs):
        v=[Q(0)]*3
        for _ in range(STEPS):
            v=add(rhs,dot(M,v))
        tail=q**STEPS*max(rhs)/(1-q)
        tails.append(tail)
        return add(v,[tail]*3)
    def source(i,j):
        return [Q(c['source_strict_absolute_mixed_derivative_caps'][f'{i},{j}']),Q(0),Q(0)]
    S,Sx,Sr,Sxx=[source(i,j) for i,j in [(0,0),(1,0),(0,1),(2,0)]]
    N0=resolvent_upper(mul(T,S))
    N1=resolvent_upper(add(mul(2,S),mul(T,add(dot(combined(1,0),N0),Sx))))
    H0=S
    H0r=Sr
    H1=add(dot(A,S),Sx)
    # Validate the producer's actual free data separately from our upper solve.
    publishedH1=list(map(Q,candidate['initial_majorants']['w_initial_velocity_A_at_zero_S0_plus_Sx0']))
    if args.inject=='initial_velocity':
        H1=add(dot(mat('A',1,0),S),Sx)
    need(publishedH1==H1,'WRONG_SECOND_TIME_INITIAL_VELOCITY')
    freeW=[max(2*H0[k]/T+2*H1[k],H0r[k]+H1[k]) for k in range(3)]
    publishedFree=list(map(Q,candidate['initial_majorants']['w_free_weighted_C1']))
    need(publishedFree==freeW,'WRONG_SECOND_TIME_FREE_BOUND')
    forcingW=add(mul(2,dot(combined(1,0),N1)),dot(combined(2,0),N0),Sxx)
    N2=resolvent_upper(add(freeW,mul(T,forcingW)))
    publishedN=[list(map(Q,candidate['weighted_C1_bounds'][key])) for key in ('U','U_x','U_xx')]
    for exact,upper in zip(publishedN,[N0,N1,N2]):
        need(all(a<=b for a,b in zip(exact,upper)),'producer resolvent exceeds independent bound')
    # Producer's recurrence must use twice the first coefficient commutator.
    pN0,pN1,_=publishedN
    commfactor=1 if args.inject=='commutator' else 2
    publishedR2=list(map(Q,candidate['commutator_forcing_majorants']['for_Uxx']))
    directR2=add(mul(commfactor,dot(combined(1,0),pN1)),dot(combined(2,0),pN0),Sxx)
    need(publishedR2==directR2,'WRONG_SECOND_COMMUTATOR_MULTIPLICITY')
    J={(0,0):mul(THETA,N0),(1,0):N0,(0,1):N0,
       (2,0):[min(a,THETA*b) for a,b in zip(N1,N2)],
       (1,1):N1,(3,0):N2,(2,1):N2}
    J[0,2]=add(J[2,0],dot(A,J[1,0]),dot(B,J[0,1]),dot(D,J[0,0]),S)
    Fx=add(dot(mat('A',1,0),J[1,0]),dot(A,J[2,0]),
           dot(mat('B',1,0),J[0,1]),dot(B,J[1,1]),
           dot(mat('D',1,0),J[0,0]),dot(D,J[1,0]),Sx)
    Fr=add(dot(mat('A',0,1),J[1,0]),dot(A,J[1,1]),
           dot(mat('B',0,1),J[0,1]),dot(B,J[0,2]),
           dot(mat('D',0,1),J[0,0]),dot(D,J[0,1]),Sr)
    J[1,2]=add(J[3,0],Fx)
    J[0,3]=add(J[2,1],Fr)
    keys=set(J)
    need(keys=={(i,j) for i in range(4) for j in range(4-i)},'ten derivative classes')
    publishedJ={tuple(map(int,k.split(','))):list(map(Q,v))
                for k,v in candidate['scaled_mixed_partial_caps'].items()}
    need(set(publishedJ)==keys,'candidate ten derivative classes')
    # An independent exact reconstruction check using the candidate's raw jets.
    pFr=add(dot(mat('A',0,1),publishedJ[1,0]),dot(A,publishedJ[1,1]),
            dot(mat('B',0,1),publishedJ[0,1]),dot(B,publishedJ[0,2]),
            dot(mat('D',0,1),publishedJ[0,0]),dot(D,publishedJ[0,1]),
            [Q(0)]*3 if args.inject=='radial_recovery' else Sr)
    need(publishedJ[0,3]==add(publishedJ[2,1],pFr),'WRONG_THIRD_RADIAL_RECOVERY')
    caps=candidate['strict_readable_scaled_total_order_caps_by_channel']
    checked=[]
    for channel,name in enumerate(('f','b','d')):
        for order in range(4):
            upper=max(v[channel] for key,v in J.items() if sum(key)==order)
            cap=Q(caps[name][order])
            need(upper<cap,'independent strict cap failure:'+name+':'+str(order))
            checked.append({'channel':name,'order':order,'strict_cap':str(cap),'verified':True})
    for key,value in J.items():
        need(all(a<=b for a,b in zip(publishedJ[key],value)),
             'published mixed bound exceeds independent Neumann bound:'+str(key))
        publishedPhysical=list(map(Q,candidate['unscaled_coordinate_mixed_partial_caps'][f'{key[0]},{key[1]}']))
        need(publishedPhysical==mul(Q(10**6)**sum(key),publishedJ[key]),'physical scaling mismatch')
    result={'status':'PASS_INDEPENDENT_NEUMANN_C3_MAJORANT_AUDIT',
            'coefficient_sha256':COEFFICIENT_SHA,
            'candidate_sha256':hashlib.sha256(args.candidate.read_bytes()).hexdigest(),
            'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'method':'16 nonnegative Neumann iterates plus q^16*||rhs||infinity/(1-q) tail; no matrix inverse or producer imports',
            'q':str(q),'number_of_mixed_derivative_classes':10,
            'number_of_component_bounds':30,'strict_order_caps_verified':checked,
            'source_and_initial_checks':['A0*S0+Sx0','max(2S0/T+2H1,S0r+H1)','second commutator factor2','rr/xrr/rrr recovery','ell to negative total derivative order'],
            'scope':'Finite rational bound audit; the accompanying review checks the underlying C3 construction and derivative identities.'}
    with args.output.open('x') as handle:
        json.dump(result,handle,indent=2,sort_keys=True);handle.write('\n')
    print(result['status'])


if __name__=='__main__':
    try:
        main()
    except (ValueError,OSError) as exc:
        print(str(exc),file=sys.stderr);sys.exit(2)
