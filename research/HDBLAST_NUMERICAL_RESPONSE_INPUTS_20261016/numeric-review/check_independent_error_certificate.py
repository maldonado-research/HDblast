#!/usr/bin/env python3
"""Independent rational replay and manufactured-identity review of C3 errors."""
import argparse
from fractions import Fraction as F
import hashlib
from math import comb, factorial
import json
from pathlib import Path


def need(ok, text):
    if not ok:
        raise ValueError(text)


def mv(a, v):
    return [sum((x*y for x, y in zip(row, v)), F(0)) for row in a]


def vsum(*vs):
    return [sum(row, F(0)) for row in zip(*vs)]


def scale(c, v):
    return [c*x for x in v]


def det(a):
    return sum((a[0][j]*(a[1][(j+1)%3]*a[2][(j+2)%3]
                -a[1][(j+2)%3]*a[2][(j+1)%3]) for j in range(3)), F(0))


def adjugate_inverse(a):
    # Cyclic cofactor construction, independent of producer elimination.
    denominator = det(a)
    need(denominator != 0, 'nonzero adjugate denominator')
    return [[(a[(j+1)%3][(i+1)%3]*a[(j+2)%3][(i+2)%3]
              -a[(j+1)%3][(i+2)%3]*a[(j+2)%3][(i+1)%3])/denominator
             for j in range(3)] for i in range(3)]


# Bivariate exact polynomial arithmetic for independent manufactured fields.
def bp(value=0):
    return {(0, 0): F(value)} if value else {}


def plus(*args):
    out = {}
    for arg in args:
        for key, value in arg.items():
            out[key] = out.get(key, F(0))+value
    return {k: v for k, v in out.items() if v}


def times(a, b):
    out = {}
    for i, v in a.items():
        for j, w in b.items():
            key = i[0]+j[0], i[1]+j[1]
            out[key] = out.get(key, F(0))+v*w
    return {k: v for k, v in out.items() if v}


def smul(c, a):
    return {k: c*v for k, v in a.items() if c*v}


def derivative(a, x=0, r=0):
    return {(i-x,j-r): v*F(factorial(i),factorial(i-x))*F(factorial(j),factorial(j-r))
            for (i,j),v in a.items() if i>=x and j>=r}


def initial(a):
    return {key: v for key, v in a.items() if key[0] == 0}


def manufactured(injection):
    e = {(2,0):F(1,2),(2,2):F(1,2),(3,0):F(1,3),(3,2):F(1,6),(4,2):F(1,24)}
    a = {(0,0):F(1),(1,0):F(1),(1,2):F(1)}
    b = {(0,1):F(1),(2,1):F(1)}
    d = {(0,2):F(1),(2,2):F(1)}
    linear = lambda field: plus(times(a,derivative(field,1)),times(b,derivative(field,r=1)),times(d,field))
    wave = lambda field: plus(derivative(field,2),smul(-1,derivative(field,r=2)))
    residual = plus(wave(e),smul(-1,linear(e)))
    for k in (1,2):
        commutator = {}
        for j in range(1,k+1):
            for coefficient, ix, ir in ((a,1,0),(b,0,1),(d,0,0)):
                term = smul(comb(k,j),times(derivative(coefficient,j),derivative(e,k-j+ix,ir)))
                if injection == 'omit_time_commutator' and k==1 and coefficient is a:
                    term = {}
                commutator = plus(commutator,term)
        need(wave(derivative(e,k)) == plus(linear(derivative(e,k)),commutator,derivative(residual,k)),
             'manufactured differentiated coefficient commutator')
    third = plus(times(initial(a),initial(residual)),initial(derivative(residual,1)))
    if injection == 'omit_initial_acceleration_offset':
        third = initial(derivative(residual,1))
    need(initial(derivative(e,3)) == third, 'manufactured stabilized third initial jet')
    rr = derivative(residual,r=1)
    if injection == 'omit_radial_residual':
        rr = {}
    recovered = plus(derivative(e,2,1),smul(-1,derivative(linear(e),r=1)),smul(-1,rr))
    need(derivative(e,r=3) == recovered, 'manufactured third radial recovery')
    need(not initial(e) and not initial(derivative(e,1)), 'manufactured zero first Cauchy jets')
    need(all(j%2==0 for (i,j) in e), 'manufactured Neumann boundary parity')
    return 6


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--producer',type=Path,default=Path(__file__).resolve().parent.parent/'numeric')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--inject',choices=('omit_time_commutator','omit_initial_acceleration_offset','omit_radial_residual'))
    args = parser.parse_args()
    count = manufactured(args.inject)
    root = args.producer
    raw_report = (root/'APPROXIMATION_ERRORS_NORMAL.json').read_bytes()
    report = json.loads(raw_report)
    pins = json.loads((root/'INPUT_PINS.json').read_text())
    inputs = {}
    for key,pin in pins.items():
        raw = (root/pin['path']).read_bytes()
        need(hashlib.sha256(raw).hexdigest()==pin['sha256'] and len(raw)==pin['bytes'],'pinned input '+key)
        inputs[key] = json.loads(raw)
    master, co, bump = (inputs[k] for k in ('master','coefficients','bump'))
    m,bd = master['exact']['interior'],master['exact']['boundary']
    ell, x, th = F(1,10**6),F(3),F(3,2)
    qcap = list(map(F,bump['dimensionless_Q_derivative_sup']))
    velocity = [x*ell,ell,ell*(qcap[1]+ell),ell*(qcap[2]+2*ell*qcap[1]+ell**2)]
    lift = list(map(F,(bd['y'],bd['yprime'],m['ysecond'],m['ythird'])))
    kernel = [x*x*F(m['forcing']['f'])/2,F(m['first_error']),max(map(F,m['second_error'].values())),max(map(F,m['third_error'].values()))]
    ew = [sum(row,F(0)) for row in zip(velocity,lift,kernel)]
    need(ew==[F(v['exact_upper']) for v in report['master_error_budget_by_scaled_total_order']['actual_minus_comparison']], 'independent actual master C3 errors')
    c = 6*ell/130
    pcap = [5*c]+[c*(5*qcap[k]+k*qcap[k-1]) for k in range(1,4)]
    need(pcap==list(map(F,report['gauge_comparison_absolute_derivative_caps_by_total_order'])), 'prefactor-complete comparator derivative caps')
    actual = list(map(F,master['source_bound_adapter']['scaled_order_bounds']))
    indices = [(i,j) for k in range(3) for i in range(k+1) for j in (k-i,)]
    mag = lambda item:max(map(lambda s:abs(F(s)),item['exact']))
    matrices = {name:{ij:[[mag(entry[f'{ij[0]},{ij[1]}']) for entry in row] for row in co['matrix_mixed_derivative_intervals'][name]] for ij in indices} for name in ('A','B','D')}
    ps,p0 = (co['source_profile_mixed_derivative_intervals'][name] for name in ('P_s','P_0'))
    residual = {}
    for i,j in indices:
        value=[F(0)]*3;source=F(0)
        for a in range(i+1):
            for b in range(j+1):
                weight=comb(i,a)*comb(j,b);k=i+j-a-b;key=f'{a},{b}'
                delta=max(abs(F(v)-(c if a+b==0 else 0)) for v in ps[key]['exact'])
                source+=weight*(delta*qcap[k+1]+mag(ps[key])*ew[k+1]+mag(p0[key])*actual[k])
                for channel in range(3):
                    value[channel]+=weight*((matrices['A'][a,b][channel][0]+matrices['B'][a,b][channel][0])*pcap[k+1]+matrices['D'][a,b][channel][0]*pcap[k])
        value[0]+=source;residual[i,j]=value
        need(source==F(report['source_approximation_defect_mixed_scaled_derivative_caps'][f'{i},{j}']['exact_upper']), 'independent mixed source defect')
        need(value==[F(v['exact_upper']) for v in report['complete_gauge_residual_mixed_scaled_derivative_caps'][f'{i},{j}']], 'independent complete gauge residual')
    A,B,D=(matrices[k] for k in ('A','B','D'))
    linear=lambda key:[[A[key][i][j]+B[key][i][j]+th*D[key][i][j] for j in range(3)] for i in range(3)]
    M=[[x*v for v in row] for row in linear((0,0))]
    rownorm=max(map(sum,M));need(rownorm<1,'true contraction')
    inv=adjugate_inverse([[F(i==j)-M[i][j] for j in range(3)] for i in range(3)])
    need(all(v>=0 for row in inv for v in row),'positive adjugate resolvent')
    r0,rx,rr,rxx=(residual[key] for key in ((0,0),(1,0),(0,1),(2,0)))
    n0=mv(inv,scale(x,r0))
    n1=mv(inv,vsum(scale(2,r0),scale(x,vsum(mv(linear((1,0)),n0),rx))))
    h1=vsum(mv(A[0,0],r0),rx)
    free2=[max(r0[i]/th+2*h1[i],rr[i]+h1[i]) for i in range(3)]
    n2=mv(inv,vsum(free2,scale(x,vsum(scale(2,mv(linear((1,0)),n1)),mv(linear((2,0)),n0),rxx))))
    for key,v in (('n0',n0),('n1',n1),('n2',n2)):
        need(v==list(map(F,report['gauge_error_resolvent'][key])), 'independent '+key+' differentiated resolvent')
    errors={(0,0):scale(th,n0),(1,0):n0,(0,1):n0,
            (2,0):[min(n1[i],th*n2[i]) for i in range(3)],
            (1,1):n1,(3,0):n2,(2,1):n2}
    errors[0,2]=vsum(errors[2,0],mv(A[0,0],n0),mv(B[0,0],n0),mv(D[0,0],errors[0,0]),r0)
    errors[1,2]=vsum(n2,mv(A[1,0],n0),mv(A[0,0],errors[2,0]),mv(B[1,0],n0),mv(B[0,0],n1),mv(D[1,0],errors[0,0]),mv(D[0,0],n0),rx)
    errors[0,3]=vsum(n2,mv(A[0,1],n0),mv(A[0,0],n1),mv(B[0,1],n0),mv(B[0,0],errors[0,2]),mv(D[0,1],errors[0,0]),mv(D[0,0],n0),rr)
    for (i,j),value in errors.items():
        need(value==[F(v['exact_upper']) for v in report['actual_minus_gauge_comparison_mixed_scaled_derivative_caps'][f'{i},{j}']], 'independent final mixed gauge derivative errors')
    out={'status':'PASS_INDEPENDENT_C3_ERROR_REVIEW','manufactured_polynomial_checks':count,
         'all_point_analytic_review':['actual master decomposition with reflected shell terms and dynamic lift through third derivative','prefactor-complete comparison derivatives','physical homogeneous mixed boundary rows survive time differentiation','nonzero initial acceleration and third initial jet retained','radial recovery retains radial residual'],
         'exact_replay':['four master errors','all six source defect derivatives','all six three-channel residual derivatives','all three error resolvent vectors by independent adjugate inverse','all ten mixed gauge derivative errors'],
         'producer_report_sha256':hashlib.sha256(raw_report).hexdigest(),
         'producer_source_sha256':hashlib.sha256((root/'check_approximation_errors.py').read_bytes()).hexdigest(),
         'input_sha256':{k:v['sha256'] for k,v in pins.items()},
         'limits':'Raw master and coordinate gauge-amplitude approximation errors only; no second response, common-gauge source error, or background point approximation is certified by this review.'}
    need(not args.output.exists(),'refuse existing review output')
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')


if __name__=='__main__':
    main()
