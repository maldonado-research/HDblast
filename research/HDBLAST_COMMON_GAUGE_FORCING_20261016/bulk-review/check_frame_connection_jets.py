#!/usr/bin/env python3
"""Exact witnesses for intrinsic sphere curvature and the rank-three slot."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random
import sys


def need(ok,msg):
    if not ok:raise ValueError(msg)


def run(mutation):
    # Sphere chart n=(x1,x2,x3,sqrt(1-|x|^2)) at its origin. The sphere
    # metric equals delta and Gamma=0, but d_l Gamma^k_ij=delta_kl delta_ij.
    # For Y=x1^2-x2^2, its coordinate fourth derivatives vanish whereas
    # its covariant fourth derivative does not.
    sphere=range(3)
    delta=lambda i,j:Q(i==j)
    Hess=[[Q(2 if i==j==0 else -2 if i==j==1 else 0) for j in sphere] for i in sphere]
    dGamma=lambda l,k,i,j:delta(k,l)*delta(i,j)
    sphere_checks=0
    for d,c,a,b in product(sphere,repeat=4):
        coordinate= -sum((dGamma(c,k,a,b)*Hess[d][k]+dGamma(d,k,a,b)*Hess[c][k]
                          +dGamma(d,k,c,a)*Hess[k][b]+dGamma(d,k,c,b)*Hess[a][k]
                          for k in sphere),Q(0))
        if mutation=='constant_normal_connection':coordinate=Q(0)
        intrinsic= -delta(c,a)*Hess[d][b]-delta(c,b)*Hess[d][a]-2*delta(a,b)*Hess[d][c]
        need(coordinate==intrinsic,'SPHERE_CONNECTION_DERIVATIVE_REQUIRED')
        sphere_checks+=1
    need(-delta(0,0)*Hess[0][0]-delta(0,0)*Hess[0][0]-2*delta(0,0)*Hess[0][0]==-8,
         'nonzero fourth sphere derivative witness')
    I=range(5)
    rng=random.Random(161064)
    get=lambda:Q(rng.randint(-3,3),rng.randint(1,3))
    h={}
    dh={};ddh={}
    for a in I:
        for b in range(a,5):h[a,b]=h[b,a]=get()
    for e in I:
        for a in I:
            for b in range(a,5):dh[e,a,b]=dh[e,b,a]=get()
    for e,f in product(I,repeat=2):
        for a in I:
            for b in range(a,5):ddh[e,f,a,b]=ddh[e,f,b,a]=get()
    a,beta,ay,beta0=Q(1,9),Q(2,17),Q(-1,10000),Q(3,20000)
    def makeK(aa,bb):
        k={(A,B,E):Q(0) for A,B,E in product(I,repeat=3)}
        k[1,0,1]=aa;k[1,1,0]=aa
        for i in (2,3,4):
            k[i,0,i]=aa;k[i,1,i]=bb;k[i,i,0]=-aa;k[i,i,1]=bb
        return k
    K=makeK(a,beta)
    dk=[]
    for e in I:dk.append(makeK(ay,-a*beta) if e==0 else makeK(0,beta0) if e==1 else makeK(0,0))
    nh={(B,C,D):dh[B,C,D]-sum((K[B,C,E]*h[E,D]+K[B,D,E]*h[C,E] for E in I),Q(0))
        for B,C,D in product(I,repeat=3)}
    # Direct two-stage conversion: differentiate the rank-two correction,
    # then act with K on every slot of the resulting rank-three tensor.
    nn={};count=0
    for A,B,C,D in product(I,repeat=4):
        dn=ddh[A,B,C,D]-sum((dk[A][B,C,E]*h[E,D]+K[B,C,E]*dh[A,E,D]
                    +dk[A][B,D,E]*h[C,E]+K[B,D,E]*dh[A,C,E] for E in I),Q(0))
        direct=dn-sum((K[A,B,E]*nh[E,C,D]+K[A,C,E]*nh[B,E,D]+K[A,D,E]*nh[B,C,E]
                       for E in I),Q(0))
        slot=Q(0) if mutation=='missing_derivative_slot' else Q(1)
        expanded=ddh[A,B,C,D]-sum((dk[A][B,C,E]*h[E,D]+dk[A][B,D,E]*h[C,E]
                    +K[B,C,E]*dh[A,E,D]+K[B,D,E]*dh[A,C,E]
                    +slot*K[A,B,E]*dh[E,C,D]+K[A,C,E]*dh[B,E,D]+K[A,D,E]*dh[B,C,E]
                    for E in I),Q(0))
        expanded+=sum((slot*K[A,B,E]*(K[E,C,F]*h[F,D]+K[E,D,F]*h[C,F])
                    +K[A,C,E]*(K[B,E,F]*h[F,D]+K[B,D,F]*h[E,F])
                    +K[A,D,E]*(K[B,C,F]*h[F,E]+K[B,E,F]*h[C,F])
                    for E,F in product(I,repeat=2)),Q(0))
        need(direct==expanded,'METRIC_SECOND_JET_DERIVATIVE_SLOT_REQUIRED')
        nn[A,B,C,D]=direct;count+=1
    k0=max(sum(abs(K[A,B,E]) for E in I) for A,B in product(I,repeat=2))
    k1=max(sum(abs(dk[A][B,C,E]) for E in I) for A,B,C in product(I,repeat=3))
    h0=max(map(abs,h.values()));d1=max(map(abs,dh.values()));d2=max(map(abs,ddh.values()))
    bound=d2+5*k0*d1+(2*k1+6*k0*k0)*h0
    need(all(abs(v)<=bound for v in nn.values()),'universal second-jet majorant')
    need(k0==a+beta and k1==max(abs(ay)+a*beta,beta0),'warped K row maxima')
    return {'status':'PASS_INTRINSIC_SPHERE_AND_FULL_DERIVATIVE_SLOT_WITNESSES',
            'sphere_fourth_derivative_comparisons':sphere_checks,'metric_second_jet_comparisons':count,
            'nonzero_sphere_normal_witness':'D_1^4 Y=-8 although Gamma=0 and partial_1^4Y=0 at chart origin',
            'second_jet_majorant':'D2+5K0D1+(2K1+6K0^2)H0',
            'scope':'Exact finite witnesses for the universal tensor identities proved in the companion review; not an actual-field numerical bound.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--inject',choices=['constant_normal_connection','missing_derivative_slot'])
    a=p.parse_args();need(not a.output.exists(),'output already exists')
    r=run(a.inject);r['checker_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'])


if __name__=='__main__':
    try:main()
    except (ValueError,OSError) as e:print(str(e),file=sys.stderr);sys.exit(2)
