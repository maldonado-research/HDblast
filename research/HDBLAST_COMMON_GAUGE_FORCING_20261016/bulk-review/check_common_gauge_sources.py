#!/usr/bin/env python3
"""Independent rational two-jet audit of natural and reduced affine Hessians.

Direct side expands inverse metrics, Christoffels and coordinate residuals.
Recipe side uses background-covariant connection differences. Generic curved,
nonnormal jets test connection derivatives; a nontrivial harmonic affine jet
tests the actual C1=0 simplification. Finite witnesses supplement the proof.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random
import sys

N=5
IND=range(N)
PAIRS=list(product(IND,repeat=2))
TRIPLES=list(product(IND,repeat=3))
QUADS=list(product(IND,repeat=4))
ERRORS={'connection_factor':'CONNECTION_HESSIAN',
        'background_connection_derivative':'COVARIANT_CONNECTION_DERIVATIVE',
        'ricci_product_factor':'NATURAL_EINSTEIN_HESSIAN',
        'stress_factor':'NATURAL_EINSTEIN_HESSIAN',
        'potential_cross_factor':'NATURAL_EINSTEIN_HESSIAN',
        'scalar_connection':'SCALAR_HESSIAN',
        'harmonic_quadratic_contraction':'QUADRATIC_GAUGE_VECTOR',
        'reduction_symmetrization':'REDUCED_EINSTEIN_HESSIAN',
        'reduction_lowering_cross':'REDUCED_EINSTEIN_HESSIAN'}


def need(ok,msg):
    if not ok:raise ValueError(msg)


def mat():return [[F(0) for _ in IND] for _ in IND]
def mm(a,b):return [[sum((a[i][k]*b[k][j] for k in IND),F(0)) for j in IND] for i in IND]
def scale(c,a):return [[c*v for v in row] for row in a]
def transpose(a):return [list(row) for row in zip(*a)]


def inverse(a):
    rows=[list(a[i])+[F(i==j) for j in IND] for i in IND]
    for j in IND:
        pivot=next((k for k in range(j,N) if rows[k][j]),None)
        need(pivot is not None,'singular witness metric')
        rows[j],rows[pivot]=rows[pivot],rows[j]
        c=rows[j][j];rows[j]=[x/c for x in rows[j]]
        for i in IND:
            if i!=j:
                c=rows[i][j];rows[i]=[x-c*y for x,y in zip(rows[i],rows[j])]
    return [row[N:] for row in rows]


def symmetric(get):
    a=mat()
    for i in IND:
        for j in range(i,N):a[i][j]=a[j][i]=get()
    return a


def secondjets(get):
    out=[[None for _ in IND] for _ in IND]
    for i in IND:
        for j in range(i,N):out[i][j]=out[j][i]=symmetric(get)
    return out


def witness(seed,harmonic):
    rng=random.Random(seed)
    get=lambda:F(rng.randint(-3,3),rng.randint(1,3))
    eta=mat()
    for i in IND:eta[i][i]=F(-1 if i==1 else 1)
    if harmonic:
        g=eta;dg=[mat() for _ in IND];ddg=[[mat() for _ in IND] for _ in IND]
        h=symmetric(get);dh=[mat() for _ in IND];ddh=[[mat() for _ in IND] for _ in IND]
        for target in (dh[1],ddh[1][1]):
            for i in (0,2,3,4):
                for j in (0,2,3,4):
                    if j>=i:target[i][j]=target[j][i]=get()
            target[4][4]=-sum(target[i][i] for i in (0,2,3))
    else:
        L=mat()
        for i in IND:
            L[i][i]=F(1)
            for j in range(i+1,N):L[i][j]=get()/4
        g=mm(transpose(L),mm(eta,L))
        dg=[symmetric(get) for _ in IND];ddg=secondjets(get)
        h=symmetric(get);dh=[symmetric(get) for _ in IND];ddh=secondjets(get)
    return g,dg,ddg,h,dh,ddh,[get() for _ in IND],symmetric(get),get(),[get() for _ in IND],symmetric(get),[get() for _ in range(4)]


def audit_one(seed,harmonic,mutation):
    g,dg,ddg,h,dh,ddh,phi1,phi2,chi,chi1,chi2,V=witness(seed,harmonic)
    gi=inverse(g)
    # Direct coefficient series: eps^2 coefficient is half the Hessian.
    gm=[g,h,mat()];dgm=[dg,dh,[mat() for _ in IND]]
    ddgm=[ddg,ddh,[[mat() for _ in IND] for _ in IND]]
    im=[gi]
    for k in (1,2):im.append(scale(-1,mm(gi,mm(h,im[k-1]))))
    dim=[]
    for k in range(3):
        dim.append([scale(-1,[[sum((im[i][a][c]*dgm[j][e][c][d]*im[k-i-j][d][b]
                         for i in range(k+1) for j in range(k-i+1)
                         for c in IND for d in IND),F(0)) for b in IND] for a in IND]) for e in IND])
    Gamma=[];dGamma=[]
    for k in range(3):
        Gamma.append({(a,b,c):sum((im[i][a][d]*(dgm[k-i][b][c][d]+dgm[k-i][c][b][d]-dgm[k-i][d][b][c])
                            for i in range(k+1) for d in IND),F(0))/2 for a,b,c in TRIPLES})
        dGamma.append({(e,a,b,c):sum((
            dim[i][e][a][d]*(dgm[k-i][b][c][d]+dgm[k-i][c][b][d]-dgm[k-i][d][b][c])+
            im[i][a][d]*(ddgm[k-i][e][b][c][d]+ddgm[k-i][e][c][b][d]-ddgm[k-i][e][d][b][c])
            for i in range(k+1) for d in IND),F(0))/2 for e,a,b,c in QUADS})
    ric2={(a,b):sum((dGamma[2][c,c,a,b]-dGamma[2][b,c,a,c] for c in IND),F(0))+
         sum((Gamma[i][c,c,d]*Gamma[2-i][d,a,b]-Gamma[i][c,b,d]*Gamma[2-i][d,a,c]
              for i in range(3) for c in IND for d in IND),F(0)) for a,b in PAIRS}
    vseries=[V[0],V[1]*chi,V[2]*chi*chi/2]
    Edirect={(a,b):2*(ric2[a,b]-chi1[a]*chi1[b]-F(2,3)*sum(vseries[k]*gm[2-k][a][b] for k in range(3))) for a,b in PAIRS}
    pgrad=[phi1,chi1,[F(0)]*N];phess=[phi2,chi2,mat()]
    Hcoord=[]
    for k in range(3):
        Hcoord.append({(a,b):phess[k][a][b]-sum((Gamma[i][c,a,b]*pgrad[k-i][c]
                                  for i in range(k+1) for c in IND),F(0)) for a,b in PAIRS})
    Sdirect=2*sum((im[k][a][b]*Hcoord[2-k][a,b] for k in range(3) for a,b in PAIRS),F(0))-V[3]*chi*chi
    # Direct relative vector, lowering and g-covariant derivative series.
    Cser=[[F(0)]*N];dCser=[{(e,a):F(0) for e,a in PAIRS}]
    for k in (1,2):
        Cser.append([sum((im[i][b][c]*Gamma[k-i][a,b,c] for i in range(k) for b,c in PAIRS),F(0)) for a in IND])
        dCser.append({(e,a):sum((dim[i][e][b][c]*Gamma[k-i][a,b,c]+im[i][b][c]*dGamma[k-i][e,a,b,c]
                              for i in range(k) for b,c in PAIRS),F(0)) for e,a in PAIRS})
    Cdown=[];dCdown=[]
    for k in range(3):
        Cdown.append([sum((gm[i][a][b]*Cser[k-i][b] for i in range(k+1) for b in IND),F(0)) for a in IND])
        dCdown.append({(e,a):sum((dgm[i][e][a][b]*Cser[k-i][b]+gm[i][a][b]*dCser[k-i][e,b]
                                  for i in range(k+1) for b in IND),F(0)) for e,a in PAIRS})
    nablaC2={(a,b):dCdown[2][a,b]-sum((Gamma[k][c,a,b]*Cdown[2-k][c]
                                     for k in range(3) for c in IND),F(0)) for a,b in PAIRS}
    RedDirect={(a,b):Edirect[a,b]-nablaC2[a,b]-nablaC2[b,a] for a,b in PAIRS}

    # Independent background-covariant recipe from metric two-jets.
    G=Gamma[0];dG=dGamma[0]
    H={(a,b,c):dh[a][b][c]-sum((G[d,a,b]*h[d][c]+G[d,a,c]*h[b][d] for d in IND),F(0)) for a,b,c in TRIPLES}
    dH={(e,a,b,c):ddh[e][a][b][c]-sum((
        (0 if mutation=='background_connection_derivative' else dG[e,d,a,b])*h[d][c]+G[d,a,b]*dh[e][d][c]+
        (0 if mutation=='background_connection_derivative' else dG[e,d,a,c])*h[b][d]+G[d,a,c]*dh[e][b][d]
        for d in IND),F(0)) for e,a,b,c in QUADS}
    P={(a,b,c):sum((gi[a][d]*(H[b,c,d]+H[c,b,d]-H[d,b,c]) for d in IND),F(0))/2 for a,b,c in TRIPLES}
    pP={(e,a,b,c):sum((dim[0][e][a][d]*(H[b,c,d]+H[c,b,d]-H[d,b,c])+
         gi[a][d]*(dH[e,b,c,d]+dH[e,c,b,d]-dH[e,d,b,c]) for d in IND),F(0))/2 for e,a,b,c in QUADS}
    NP={(e,a,b,c):pP[e,a,b,c]+sum((G[a,e,d]*P[d,b,c]-G[d,e,b]*P[a,d,c]-G[d,e,c]*P[a,b,d]
                                               for d in IND),F(0)) for e,a,b,c in QUADS}
    hu=mm(gi,h);hi=mm(hu,gi)
    Hup=[mm(gi,[[H[e,a,b] for b in IND] for a in IND]) for e in IND]
    Hboth=[mm(Hup[e],gi) for e in IND]
    factor=F(-1) if mutation=='connection_factor' else F(-2)
    J={(a,b,c):factor*sum((hu[a][d]*P[d,b,c] for d in IND),F(0)) for a,b,c in TRIPLES}
    NJ={(e,a,b,c):-2*sum((Hup[e][a][d]*P[d,b,c]+hu[a][d]*NP[e,d,b,c] for d in IND),F(0)) for e,a,b,c in QUADS}
    count=0
    def equal(a,b,msg):
        nonlocal count
        need(a==b,msg);count+=1
    for key in TRIPLES:
        equal(P[key],Gamma[1][key],'FIRST_CONNECTION_VARIATION')
        equal(J[key],2*Gamma[2][key],'CONNECTION_HESSIAN')
    for e,a,b,c in QUADS:
        direct=pP[e,a,b,c]
        equal(direct,dGamma[1][e,a,b,c],'COVARIANT_CONNECTION_DERIVATIVE')
    cv=[sum((gi[b][c]*P[a,b,c] for b,c in PAIRS),F(0)) for a in IND]
    Ncv={(e,a):sum((gi[b][c]*NP[e,a,b,c] for b,c in PAIRS),F(0)) for e,a in PAIRS}
    qfactor=0 if mutation=='harmonic_quadratic_contraction' else -2
    QC=[-2*sum((hu[a][d]*cv[d] for d in IND),F(0))+qfactor*sum((hi[b][c]*P[a,b,c] for b,c in PAIRS),F(0)) for a in IND]
    NQC={(e,a):-2*sum((Hup[e][a][d]*cv[d]+hu[a][d]*Ncv[e,d] for d in IND),F(0))
          -2*sum((Hboth[e][b][c]*P[a,b,c]+hi[b][c]*NP[e,a,b,c] for b,c in PAIRS),F(0)) for e,a in PAIRS}
    for a in IND:equal(QC[a],2*Cser[2][a],'QUADRATIC_GAUGE_VECTOR')
    for e,a in PAIRS:
        equal(NQC[e,a],2*dCser[2][e,a]+sum((G[a,e,d]*QC[d] for d in IND),F(0)),'QUADRATIC_GAUGE_DERIVATIVE')
    if harmonic:
        need(all(x==0 for x in cv) and all(x==0 for x in Ncv.values()),'harmonic witness constraints')
        need(any(P.values()) and any(NP.values()),'nontrivial harmonic witness')
    pp_factor=1 if mutation=='ricci_product_factor' else 2
    stress_factor=1 if mutation=='stress_factor' else 2
    cross_factor=1 if mutation=='potential_cross_factor' else 2
    Erecipe={(a,b):sum((NJ[c,c,a,b]-NJ[b,c,a,c] for c in IND),F(0))+
       pp_factor*sum((P[c,c,d]*P[d,a,b]-P[c,b,d]*P[d,a,c] for c in IND for d in IND),F(0))-
       stress_factor*chi1[a]*chi1[b]-F(2,3)*(V[2]*chi*chi*g[a][b]+cross_factor*V[1]*chi*h[a][b]) for a,b in PAIRS}
    Hphi={(a,b):phi2[a][b]-sum((G[c,a,b]*phi1[c] for c in IND),F(0)) for a,b in PAIRS}
    Hchi={(a,b):chi2[a][b]-sum((G[c,a,b]*chi1[c] for c in IND),F(0)) for a,b in PAIRS}
    hhi=mm(hu,hi)
    scalar_factor=0 if mutation=='scalar_connection' else 2
    Srecipe=(2*sum((hhi[a][b]*Hphi[a,b]-hi[a][b]*Hchi[a,b] for a,b in PAIRS),F(0))+
        scalar_factor*sum((hi[a][b]*P[c,a,b]*phi1[c] for a,b,c in TRIPLES),F(0))+
        2*sum((hu[c][d]*cv[d]*phi1[c] for c,d in PAIRS),F(0))-
        2*sum((cv[c]*chi1[c] for c in IND),F(0))-V[3]*chi*chi)
    equal(Srecipe,Sdirect,'SCALAR_HESSIAN')
    NQdown={(a,b):sum((g[b][c]*NQC[a,c] for c in IND),F(0)) for a,b in PAIRS}
    lowering={(a,b):sum((H[a,b,c]*cv[c]+h[b][c]*Ncv[a,c] for c in IND),F(0)) for a,b in PAIRS}
    cvdown=[sum((g[a][c]*cv[c] for c in IND),F(0)) for a in IND]
    sym=F(1) if mutation=='reduction_symmetrization' else F(1,2)
    lower_factor=0 if mutation=='reduction_lowering_cross' else 1
    gaugeH={(a,b):sym*(NQdown[a,b]+NQdown[b,a])+lower_factor*(lowering[a,b]+lowering[b,a])-
                 2*sum((P[c,a,b]*cvdown[c] for c in IND),F(0)) for a,b in PAIRS}
    for a,b in PAIRS:
        equal(Erecipe[a,b],Edirect[a,b],'NATURAL_EINSTEIN_HESSIAN')
        equal(Erecipe[a,b]-gaugeH[a,b],RedDirect[a,b],'REDUCED_EINSTEIN_HESSIAN')
        if harmonic:
            equal(gaugeH[a,b],(NQdown[a,b]+NQdown[b,a])/2,'HARMONIC_REDUCTION_SIMPLIFICATION')
    return {'seed':seed,'curved_nonnormal_background':not harmonic,'first_harmonic_field':harmonic,'exact_scalar_comparisons':count}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--inject',choices=sorted(ERRORS))
    args=parser.parse_args()
    need(not args.output.exists(),'output already exists')
    reports=[audit_one(seed,harm,args.inject) for seed,harm in [(161061,False),(161062,False),(161063,True)]]
    out={'status':'PASS_INDEPENDENT_NATURAL_AND_REGISTERED_REDUCED_SOURCE_JETS',
         'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'dimension':N,
         'witnesses':reports,'comparison_count':sum(r['exact_scalar_comparisons'] for r in reports),
         'convention':'second epsilon derivative of affine metric/scalar path; Ered=E-nabla_(a C_b), symmetrization one half; fixed reference connection',
         'scope':'Exact finite curved/nonnormal and harmonic local-jet witnesses. Universal identities are proved in the companion derivation; no actual-source norm is evaluated by this algebra checker.'}
    with args.output.open('x') as f:json.dump(out,f,indent=2,sort_keys=True);f.write('\n')
    print(out['status'])


if __name__=='__main__':
    try:main()
    except (ValueError,OSError) as e:
        print(str(e),file=sys.stderr);sys.exit(2)
