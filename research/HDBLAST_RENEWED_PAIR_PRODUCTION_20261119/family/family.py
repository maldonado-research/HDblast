"""Bound the best chemical occupation of two independently adjustable late gates."""
import argparse
from fractions import Fraction as F
from itertools import product
import json
from math import comb
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


def multiply(a,b):
    out={}
    for (i,j),u in a.items():
        for (k,l),v in b.items():
            ij=(i+k,j+l)
            out[ij]=out.get(ij,F(0))+u*v
    return out


def power(p,n):
    result={(0,0):F(1)}
    for _ in range(n):result=multiply(result,p)
    return result


def polynomial(matrix,lo):
    def value(i,j):return F(float(matrix[i,j]))/F(float(matrix[0,0]))
    functions=({(0,0):F(1)},{(0,0):F(lo),(1,0):F(1-lo)},
               {(0,0):F(lo),(0,1):F(1-lo)})
    p={}
    for i,j in product(range(3), repeat=2):
        for key,coef in multiply(functions[i],functions[j]).items():
            p[key]=p.get(key,F(0))+coef*value(i,j)
    return p


def bernstein(p,degree):
    rows=[]
    for i in range(degree+1):
        row=[]
        for j in range(degree+1):
            row.append(sum(v*F(comb(i,a),comb(degree,a))*F(comb(j,b),comb(degree,b))
                           for (a,b),v in p.items() if a<=i and b<=j))
        rows.append(tuple(row))
    return tuple(rows)


def subdivide(net,axis):
    if axis:
        l,r=subdivide(tuple(zip(*net)),False)
        return tuple(zip(*l)),tuple(zip(*r))
    left,right=[],[]
    for row in net:
        a,b=[row[0]],[row[-1]]
        while len(row)>1:
            row=tuple((i+j)/2 for i,j in zip(row,row[1:]))
            a.append(row[0]);b.append(row[-1])
        left.append(tuple(a));right.append(tuple(reversed(b)))
    return tuple(left),tuple(right)


def certify(n,e,low,upper):
    # All numbers below are exact rationals for the stored binary64 Gram entries.
    npoly,epoly=polynomial(n,low),polynomial(e,low)
    limit=F(upper)**4
    residual=power(npoly,4)
    for key,v in power(epoly,3).items():residual[key]=residual.get(key,F(0))-limit*v
    boxes=[(bernstein(npoly,2),bernstein(epoly,2),bernstein(residual,8),0)]
    visited=0
    while boxes:
        bn,be,bp,depth=boxes.pop()
        visited+=1
        elo=min(map(min,be))
        if elo>0 and max(map(max,bn))**4<=limit*elo**3:continue
        if max(map(max,bp))<=0:continue
        if depth>=110:raise ArithmeticError("Upper ratio could not be certified")
        children=[subdivide(net,depth%2) for net in (bn,be,bp)]
        boxes.extend((*subbox,depth+1) for subbox in zip(*children))
    return visited


def score(n,e,params):
    v=np.r_[1.,params]
    rn,re=(v@n@v)/n[0,0],(v@e@v)/e[0,0]
    return float(rn/re**.75)


def maximum(n,e,low):
    def objective(x):
        v=np.r_[1.,x]
        rn,re=v@n@v,v@e@v
        s=score(n,e,x)
        return -s,-s*(2*(n@v)[1:]/rn-1.5*(e@v)[1:]/re)
    points=[(0.,0.)]
    for x0 in product(sorted({0.,low,1.}),repeat=2):
        sol=minimize(objective,x0,jac=True,bounds=[(low,1.)]*2,method="L-BFGS-B",
                     options={"gtol":1e-13,"ftol":1.e-16,"maxiter":200})
        points.extend([sol.x,x0])
    best=max(points,key=lambda x:score(n,e,x))
    achieved=score(n,e,best)
    upper=achieved+1.e-11
    boxes=certify(n,e,low,upper)
    return dict(coefficients=list(map(float,best)),ratio_to_A=achieved,certified_upper=upper,
                boxes=boxes)


def analyze(emission):
    d=json.loads(emission.read_bytes())
    indices=[d["matrix_basis_order"].index(g) for g in ("A","first","last")]
    results={}
    for cutoff in ("120","1920"):
        walls={}
        for arm in ("quiet","plus","minus","far"):
            src=d["arms"][arm]["cutoffs"][cutoff]
            n,e=(np.array(src[key]["real"])[np.ix_(indices,indices)] for key in
                 ("pair_number_matrix","both_sign_conformal_energy_matrix"))
            n,e=((v+v.T)*.5 for v in (n,e))
            assert min(np.linalg.eigvalsh(e))>0 and min(np.linalg.eigvalsh(n))>0
            walls[arm]={}
            for label,low in (("signed",-1),("nonnegative",0)):
                result=maximum(n,e,low)
                lookup=result if arm=="quiet" else walls["quiet"][label]
                result["quiet_coefficients_on_this_history"]=score(n,e,lookup["coefficients"])
                walls[arm][label]=result
        results[cutoff]=walls
    return results


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--emission",type=Path,required=True)
    cli.add_argument("--output",type=Path,required=True)
    args=cli.parse_args()
    args.output.write_text(json.dumps(analyze(args.emission),indent=2,allow_nan=False)+"\n")
