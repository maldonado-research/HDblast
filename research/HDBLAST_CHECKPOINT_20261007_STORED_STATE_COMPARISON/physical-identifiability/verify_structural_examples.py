#!/usr/bin/env python3
"""Exact manufactured identifiability examples; no physical inputs or callbacks."""
from fractions import Fraction as Q
from pathlib import Path
import argparse,hashlib,json,sys

def matrix(rows):return [[Q(x) for x in row] for row in rows]
def mv(a,x):return [sum((v*w for v,w in zip(row,x)),Q(0)) for row in a]
def mm(a,b):return [[sum((aij*b[j][k] for j,aij in enumerate(row)),Q(0)) for k in range(len(b[0]))] for row in a]
def rank(a):
    a=[row[:] for row in a];i=0
    for j in range(len(a[0])):
        candidate=next((r for r in range(i,len(a)) if a[r][j]),None)
        if candidate is None:continue
        a[i],a[candidate]=a[candidate],a[i]
        pivot=a[i][j];a[i]=[v/pivot for v in a[i]]
        for r in range(len(a)):
            if r!=i:
                factor=a[r][j];a[r]=[v-factor*w for v,w in zip(a[r],a[i])]
        i+=1
        if i==len(a):break
    return i

def run():
    checks=[]
    def check(name,predicate):
        if not predicate:raise ValueError(name)
        checks.append(name)
    # Four measured linear responses, two effective amplitudes, three proposed
    # microscopic parameters. Nothing here evaluates the HDBLAST source.
    x=matrix([[1,0],[0,1],[1,1],[2,-1]])
    a=matrix([[1,0,1],[0,1,-1]])
    theta=list(map(Q,[2,3,5]));null=list(map(Q,[-1,1,1]))
    theta2=[v+7*w for v,w in zip(theta,null)]
    response=mm(x,a)
    check('two independent effective profiles are recoverable',rank(x)==2)
    check('three microscopic parameters give only rank two',rank(response)==2)
    check('exact nonzero microscopic null direction',any(null) and mv(a,null)==[Q(0),Q(0)])
    check('different microscopic parameter vectors give identical four responses',theta!=theta2 and mv(response,theta)==mv(response,theta2))
    nuisance=mv(a,theta)
    check('unrestricted two-amplitude 4D rival exactly replicates proposed bulk response',mv(x,nuisance)==mv(response,theta))
    data=list(map(Q,[3,-1,2,4]));weights=list(map(Q,[1,2,3,4]))
    def likelihood_quadratic(mean):return sum((w*(y-z)**2 for w,y,z in zip(weights,data,mean)),Q(0))
    check('same covariance Gaussian likelihood exponent is identical',likelihood_quadratic(mv(x,nuisance))==likelihood_quadratic(mv(response,theta)))
    ward=matrix([[1,1,-1,0]])
    different=[nuisance[0]+1,nuisance[1]]
    check('Ward row annihilates entire response family',mm(ward,x)==[[Q(0),Q(0)]])
    check('incorrect effective source passes same Ward row',mv(x,different)!=mv(x,nuisance) and mv(ward,mv(x,different))==[Q(0)])
    collapsed=matrix([[1,2],[3,6]])
    check('two profiles need not be distinguishable by chosen measurement',rank(collapsed)==1 and mv(collapsed,[Q(2),Q(-1)])==[Q(0),Q(0)])
    lam=Q(2,7);joint=[lam,6*lam]
    check('restricted joint coefficient relation is falsifiable',joint[1]-6*joint[0]==0 and Q(1)-6*joint[0]!=0)
    check('4D rival with free couplings can reproduce restricted joint relation',mv(matrix([[1,0],[0,1]]),joint)==joint)
    check('restricted zero-second-channel rival cannot reproduce nonzero joint signal',joint[1]!=0)
    first,second=Q(1,5),Q(2,7)
    check('common endpoint leading phase relation holds in manufactured variables',3*first-3*first==0)
    check('different endpoint histories invalidate imposed three-to-one phase relation',3*first-3*second!=0)
    return {'status':'PASS_EXACT_MANUFACTURED_STRUCTURAL_EXAMPLES','scope':'finite rational linear algebra only; no physical parameter inference',
            'checks':checks,'check_count':len(checks),'source_callbacks':0,'saved_array_reads':0,'saved_array_decodes':0,
            'trajectory_callbacks':0,'data_likelihood_evaluations':0,'external_novelty':'NOT_ASSESSED','higher_dimensional_cause':'NOT_ESTABLISHED'}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():raise ValueError('Fresh output required')
    result=run();science=json.dumps(result,sort_keys=True,separators=(',',':')).encode()
    result={'scientific':result,'scientific_sha256':hashlib.sha256(science).hexdigest(),'python_optimization':sys.flags.optimize,
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    with args.output.open('x') as handle:json.dump(result,handle,indent=2,sort_keys=True);handle.write('\n')
    print(json.dumps({'status':result['scientific']['status'],'checks':result['scientific']['check_count'],'scientific_sha256':result['scientific_sha256']}))
