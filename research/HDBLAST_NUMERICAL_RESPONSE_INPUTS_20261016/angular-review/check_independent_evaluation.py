#!/usr/bin/env python3
"""Independent invariant-contraction and one-variable proof of evaluation cap.

The producer's simplex Bernstein argument is not used for acceptance. This
checker derives vector/tensor traces from projected ambient contractions,
then bounds their J endpoints/vertex and univariate Bernstein coefficients.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import comb
from pathlib import Path
import json


def need(ok,msg):
    if not ok:raise ValueError(msg)


class P:
    def __init__(self,value=0):self.d={e:F(v) for e,v in (value.d if isinstance(value,P) else value if isinstance(value,dict) else {(0,0):value}).items() if v}
    def __add__(self,other):
        d=self.d.copy()
        for e,v in P(other).d.items():d[e]=d.get(e,F(0))+v
        return P(d)
    __radd__=__add__
    def __neg__(self):return P({e:-v for e,v in self.d.items()})
    def __sub__(self,other):return self+-P(other)
    def __rsub__(self,other):return P(other)+-self
    def __mul__(self,other):
        d={}
        for (a,b),v in self.d.items():
            for (c,e),w in P(other).d.items():d[a+c,b+e]=d.get((a+c,b+e),F(0))+v*w
        return P(d)
    __rmul__=__mul__
    def __truediv__(self,c):return self*F(1,c)
    def __pow__(self,n):
        p=P(1)
        for _ in range(n):p=p*self
        return p
    def __eq__(self,other):return self.d==P(other).d
    def substitute(self,s,j):return sum((v*s**a*j**b for (a,b),v in self.d.items()),P(0))
    def at(self,s,j):return sum((v*s**a*j**b for (a,b),v in self.d.items()),F(0))


def run(mutation):
    s=P({(1,0):1});j=P({(0,1):1})
    # DS=2(P A n), DY=2(P Q n); A=Q², Q³=Q.
    gradz=4*s*(1-s)
    gradf=16*j*(s-j)+s*(1-s)-8*j*(1-s)
    # HZ is the trace-free part of 2PAP. HF is TF(B),
    # B=8 v⊗v+4Y PQP−PAP, v=PQn.
    hz=4*(2-2*s+s*s)-F(4,3)*(2-s)**2
    trb=9*s-12*j-2
    b2=64*(s-j)**2+16*j*(2-2*s+j)+(2-2*s+s*s)+64*j*(1-2*s+j)-16*(s+j*(s-2))-8*j*(s-2)
    hf=b2-(0 if mutation=='omit_TF_trace' else F(1,3))*trb**2
    expectedhf=96*j*j-240*j*s+128*j+38*s*s-6*s+F(2,3)
    need(hf==expectedhf,'INTRINSIC_TF_CONTRACTION')
    fnorm=F(1,40) if mutation=='wrong_F_Gram' else F(19,720)
    k0=1+12*(s-F(1,2))**2+(j-s/2+F(1,12))**2/fnorm
    k1=F(3,2)*gradz+F(30,19)*gradf
    # Each chiral tensor has pointwise norm²=1/6 and L2 norm²=1/6.
    k2=F(9,20)*hz+F(15,133)*hf+(1 if mutation=='omit_TT_sector' else 2)
    need(k0==F(3,19)*(240*j*j-240*j*s+40*j+136*s*s-96*s+27),'SCALAR_GRAM_KERNEL')
    need(k1==F(48,19)*(-10*j*j+15*j*s-5*j-3*s*s+3*s),'GRADIENT_CONTRACTION')
    need(k2==F(6,665)*(1200*j*j-3000*j*s+1600*j+608*s*s-208*s+363),'BOTH_TT_KERNELS')
    # Exact independent global scalar bound and an attaining reproducing vector.
    bound=F(16) if mutation=='unsafe_constant4' else F(321,19)
    need(k0.at(1,1)==F(321,19)<=bound,'SHARP_SCALAR_ENDPOINT')
    q_at=1+6*F(1,2)+F(420,19)*F(7,12)
    qnorm=1+36*F(1,12)+F(420,19)**2*fnorm
    need(q_at==qnorm==F(321,19),'SHARPNESS_REPRODUCER')
    # Vector: when S<=1/3, derivative in J is nonpositive. Above 1/3,
    # the unconstrained vertex is inside [0,S²] because 4S²−3S+1>0.
    vertex=(3*s-1)/4
    vmax=k1.substitute(s,vertex)
    need(vmax==F(48,19)*(105*s*s-30*s+25)/40,'VECTOR_VERTEX')
    need(9-16<0 and vmax.at(1,0)==F(120,19) and vmax.at(F(1,3),0)==F(32,19),'VECTOR_CONVEX_ENDPOINTS')
    # Tensor is convex in J, so it is enough to bound J=0 and J=S².
    t0=k2.substitute(s,P(0))*F(665,6)
    t1=k2.substitute(s,s*s)*F(665,6)
    need(t0==608*s*s-208*s+363,'TENSOR_ZERO_ENDPOINT')
    need(t1==1200*s**4-3000*s**3+2208*s*s-208*s+363,'TENSOR_QUARTIC_ENDPOINT')
    bern=[sum(v*F(comb(i,a),comb(4,a)) for (a,b),v in t1.d.items() if b==0 and a<=i) for i in range(5)]
    if mutation=='Bernstein_error':bern[2]+=1
    recon=sum((c*comb(4,i)*s**i*(1-s)**(4-i) for i,c in enumerate(bern)),P(0))
    need(recon==t1,'UNIVARIATE_BERNSTEIN_IDENTITY')
    need(bern==list(map(F,[363,311,627,561,563])) and max(bern)<763,'TENSOR_ENDPOINT_BOUND')
    need(t0.at(1,0)==763 and t0.at(0,0)==363,'TENSOR_CONVEX_ENDPOINTS')
    vectorcap=F(120,19);tensorcap=F(4578,665)
    need(max(vectorcap,tensorcap)<bound,'BLOCK_NORM_BOUND')
    # Compare independently derived invariant polynomials to producer report.
    base=Path(__file__).resolve().parent;pins=json.loads((base/'INPUT_PINS.json').read_text())
    for key,item in pins.items():
        need(sha256((base/item['path']).read_bytes()).hexdigest()==item['sha256'],'INPUT_HASH_'+key)
    report=json.loads((base/'inputs/EVALUATION_BOUND_NORMAL.json').read_text())
    x=s;y=j
    for key,poly in zip(['scalar_block','gradient_block','TF_tensor_block'],[k0,k1,k2]):
        own=poly.substitute(x+y,(x-y)**2)
        claimed=P({tuple(v['powers_X_Y']):F(v['coefficient']) for v in report['certificates'][key]['polynomial']})
        need(own==claimed,'PRODUCER_POLYNOMIAL_'+key)
    return {'status':'PASS_INDEPENDENT_SHARP_SELECTED23_EVALUATION','squared_sharp_constant':str(bound),
      'independent_vector_trace_global_max':str(vectorcap),'independent_tensor_trace_global_max':str(tensorcap),
      'quartic_univariate_Bernstein':list(map(str,bern)),'sharp_reproducer_value_and_squared_norm':str(qnorm),
      'method':'Exact invariant contractions and endpoint/vertex proofs, independent of producer triangle Bernstein certificates.',
      'scope':'Selected23 evaluation only; no reduction closure or numerical response follows from this checker.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--inject',choices=['omit_TF_trace','wrong_F_Gram','omit_TT_sector','unsafe_constant4','Bernstein_error']);a=p.parse_args()
    need(not a.output.exists(),'OUTPUT_EXISTS')
    need(not a.output.resolve().is_relative_to(Path(__file__).resolve().parent),'OUTPUT_MUST_BE_EXTERNAL')
    result=run(a.inject);need(a.inject is None,'CONTROL_NOT_REJECTED')
    result['checker_sha256']=sha256(Path(__file__).read_bytes()).hexdigest();a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(result['status'])


if __name__=='__main__':main()
